"""Generic post-reduction cleanup for Goodman 1D coadds.

The 3580-3720 A "ghost" / blue-edge region on Goodman 400_M1/M2 often leaves
residual spikes after flat-fielding.  Those spikes are bright enough to flip
SNID-SAGE type assignments (e.g. Ia -> Ic) even when the rest of the spectrum is
fine.  :func:`clean_goodman_coadd_blue_end` soft-masks and replaces pathological
pixels using a continuum model anchored just redward of the ghost band.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter

logger = logging.getLogger(__name__)

# Empirically problematic band on Goodman red 400_M1/M2 (see ghost_correction tests).
# Extend below 3500 A: coadds often retain wild edge pixels from 3050-3500 A that
# are outside the classical 3580-3720 ghost but still poison SNID.
GHOST_WAVE_LO = 3000.0
GHOST_WAVE_HI = 3800.0
REF_WAVE_LO = 3800.0
REF_WAVE_HI = 4500.0


def _continuum_model(
    wave: np.ndarray,
    flux: np.ndarray,
    *,
    ref_lo: float,
    ref_hi: float,
    medfilt: int = 51,
) -> np.ndarray:
    """Build a smooth continuum and extrapolate it across the full wave grid."""
    ref = (
        np.isfinite(wave)
        & np.isfinite(flux)
        & (wave >= ref_lo)
        & (wave <= ref_hi)
        & (flux > 0)
    )
    if int(np.count_nonzero(ref)) < 20:
        # Fall back to a global positive median so callers still get a vector.
        med = float(np.nanmedian(flux[np.isfinite(flux) & (flux > 0)]))
        return np.full_like(flux, med if np.isfinite(med) and med > 0 else 1.0)

    # Median-filter the reference band, then fit log-flux vs wave (linear).
    f_ref = flux[ref].copy()
    w_ref = wave[ref]
    k = min(medfilt, max(5, (int(np.count_nonzero(ref)) // 2) * 2 - 1))
    f_smooth = median_filter(f_ref, size=k)
    good = f_smooth > 0
    if int(np.count_nonzero(good)) < 10:
        med = float(np.nanmedian(f_ref[f_ref > 0]))
        return np.full_like(flux, med if np.isfinite(med) and med > 0 else 1.0)

    coeff = np.polyfit(w_ref[good], np.log(f_smooth[good]), 1)
    return np.exp(np.polyval(coeff, wave))


def _repair_invalid_wavelengths(
    wave: np.ndarray,
    flux: np.ndarray,
    *,
    wave_grid_mid: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Interpolate non-positive / non-finite wavelengths; return (wave, mask_flags, n_fixed).

    PypeIt coadds occasionally write ``wave=0`` samples inside the grid. SNID-SAGE
    then reports a 0-lambdamax range and mis-classifies the spectrum. Prefer
    ``wave_grid_mid`` when available, else linear interpolation from neighbors.
    """
    wave = wave.copy()
    bad = ~np.isfinite(wave) | (wave <= 0)
    n_fixed = int(np.count_nonzero(bad))
    if n_fixed == 0:
        return wave, bad, 0

    if wave_grid_mid is not None:
        mid = np.asarray(wave_grid_mid, dtype=float)
        use_mid = bad & np.isfinite(mid) & (mid > 0)
        wave[use_mid] = mid[use_mid]
        bad = ~np.isfinite(wave) | (wave <= 0)

    good = np.isfinite(wave) & (wave > 0)
    if int(np.count_nonzero(good)) >= 2 and np.any(bad):
        wave[bad] = np.interp(np.flatnonzero(bad), np.flatnonzero(good), wave[good])
    return wave, bad, n_fixed


def clean_goodman_coadd_blue_end(
    coadd_fits: str | Path,
    *,
    ghost_lo: float = GHOST_WAVE_LO,
    ghost_hi: float = GHOST_WAVE_HI,
    ref_lo: float = REF_WAVE_LO,
    ref_hi: float = REF_WAVE_HI,
    nsig: float = 4.0,
    ratio_hi: float = 1.8,
    ratio_lo: float = 0.25,
    inplace: bool = True,
) -> dict[str, Any]:
    """Mask/replace pathological blue-end pixels and repair invalid wavelengths.

    Parameters
    ----------
    coadd_fits
        Path to a PypeIt ``coadd1d`` product (``SPECTRUM`` BinTable with
        ``wave`` / ``flux`` / ``ivar`` / ``sigma`` / ``mask``).
    ghost_lo, ghost_hi
        Wavelength window cleaned (A).
    ref_lo, ref_hi
        Continuum reference band used to model the blue end (A).
    nsig
        Robust sigma threshold for residual clipping inside the ghost band.
    ratio_hi, ratio_lo
        Additional flux/continuum ratio cuts (unitless).
    inplace
        If True, rewrite ``coadd_fits``; otherwise write ``*_blueclean.fits``.

    Returns
    -------
    dict
        ``path``, ``n_flagged``, ``n_ghost``, ``n_wave_fixed``, ``ghost_lo``, ``ghost_hi``.
    """
    path = Path(coadd_fits)
    with fits.open(path, mode='readonly') as hdul:
        if 'SPECTRUM' not in hdul:
            raise ValueError(f'{path} has no SPECTRUM extension')
        spec = hdul['SPECTRUM']
        wave = np.asarray(spec.data['wave'], dtype=float).copy()
        flux = np.asarray(spec.data['flux'], dtype=float).copy()
        names = set(spec.data.dtype.names or ())
        wave_grid_mid = (
            np.asarray(spec.data['wave_grid_mid'], dtype=float)
            if 'wave_grid_mid' in names
            else None
        )
        ivar = (
            np.asarray(spec.data['ivar'], dtype=float).copy()
            if 'ivar' in names
            else None
        )
        sigma = (
            np.asarray(spec.data['sigma'], dtype=float).copy()
            if 'sigma' in names
            else None
        )
        mask = (
            np.asarray(spec.data['mask'], dtype=np.int64).copy()
            if 'mask' in names
            else np.ones(wave.shape, dtype=np.int64)
        )
        primary = hdul[0].header.copy()
        spec_header = spec.header.copy()
        table_data = np.array(spec.data)

    wave, wave_bad, n_wave_fixed = _repair_invalid_wavelengths(
        wave, flux, wave_grid_mid=wave_grid_mid
    )
    if n_wave_fixed:
        mask[wave_bad] = 0
        if ivar is not None:
            ivar[wave_bad] = 0.0
        logger.warning(
            'Repaired %s invalid wavelength sample(s) in %s',
            n_wave_fixed,
            path.name,
        )

    ghost = np.isfinite(wave) & (wave >= ghost_lo) & (wave <= ghost_hi)
    n_ghost = int(np.count_nonzero(ghost))
    n_flagged = 0
    if n_ghost == 0:
        logger.info('No ghost-band pixels in %s; wavelength repair only', path.name)
    else:
        cont = _continuum_model(wave, flux, ref_lo=ref_lo, ref_hi=ref_hi)
        resid = flux - cont
        ref = np.isfinite(wave) & (wave >= ref_lo) & (wave <= ref_hi)
        ref_resid = resid[ref & np.isfinite(resid)]
        if ref_resid.size:
            mad = float(np.nanmedian(np.abs(ref_resid - np.nanmedian(ref_resid))))
            scale = 1.4826 * mad if mad > 0 else float(np.nanstd(ref_resid))
        else:
            scale = float(np.nanstd(resid[ghost]))
        if not np.isfinite(scale) or scale <= 0:
            scale = float(np.nanmedian(np.abs(flux[ref]))) * 0.05 if np.any(ref) else 1.0

        ratio = np.full_like(flux, np.nan)
        okc = cont > 0
        ratio[okc] = flux[okc] / cont[okc]
        flagged = ghost & np.isfinite(flux) & (
            (np.abs(resid) > nsig * scale)
            | (ratio > ratio_hi)
            | (ratio < ratio_lo)
            | (flux <= 0)
        )
        n_flagged = int(np.count_nonzero(flagged))
        if n_flagged:
            flux[flagged] = cont[flagged]
            mask[flagged] = 0
            if ivar is not None:
                ivar[flagged] = np.minimum(ivar[flagged], 1.0 / (scale * scale + 1e-30)) * 0.01
            if sigma is not None:
                sigma[flagged] = np.maximum(sigma[flagged], scale)

    table_data['wave'] = wave
    table_data['flux'] = flux
    if 'mask' in table_data.dtype.names:
        table_data['mask'] = mask
    if ivar is not None and 'ivar' in table_data.dtype.names:
        table_data['ivar'] = ivar
    if sigma is not None and 'sigma' in table_data.dtype.names:
        table_data['sigma'] = sigma

    out_path = path if inplace else path.with_name(path.stem + '_blueclean.fits')
    primary['PASSTABC'] = (
        True,
        'PASSTA blue-end ghost clean applied',
    )
    primary['BBCNFLG'] = (n_flagged, 'N pixels replaced in blue clean')
    primary['BBCNWAV'] = (n_wave_fixed, 'N invalid wavelengths repaired')
    primary['BBCGLO'] = (ghost_lo, 'Blue clean ghost wave lo [A]')
    primary['BBCGHI'] = (ghost_hi, 'Blue clean ghost wave hi [A]')

    hdu_list = fits.HDUList(
        [
            fits.PrimaryHDU(header=primary),
            fits.BinTableHDU(data=table_data, header=spec_header, name='SPECTRUM'),
        ]
    )
    hdu_list.writeto(out_path, overwrite=True)
    logger.info(
        'Blue-end clean %s: flagged %s/%s ghost-band pixels (%.0f-%.0f A); '
        'repaired %s bad wavelengths',
        out_path.name,
        n_flagged,
        n_ghost,
        ghost_lo,
        ghost_hi,
        n_wave_fixed,
    )
    return {
        'path': str(out_path),
        'n_flagged': n_flagged,
        'n_ghost': n_ghost,
        'n_wave_fixed': n_wave_fixed,
        'ghost_lo': ghost_lo,
        'ghost_hi': ghost_hi,
    }
