"""Tests for Goodman blue-end coadd cleanup."""

from __future__ import annotations

import numpy as np
from astropy.io import fits

from passta.spectrum_cleanup import clean_goodman_coadd_blue_end


def _write_fake_coadd(path, wave, flux):
    primary = fits.PrimaryHDU()
    table = fits.BinTableHDU.from_columns(
        [
            fits.Column(name='wave', format='D', array=wave),
            fits.Column(name='wave_grid_mid', format='D', array=wave),
            fits.Column(name='flux', format='D', array=flux),
            fits.Column(name='ivar', format='D', array=np.ones_like(flux)),
            fits.Column(name='sigma', format='D', array=np.ones_like(flux)),
            fits.Column(name='mask', format='K', array=np.ones(wave.shape, dtype=np.int64)),
        ],
        name='SPECTRUM',
    )
    fits.HDUList([primary, table]).writeto(path, overwrite=True)


def test_clean_goodman_coadd_blue_end_flags_spike(tmp_path):
    wave = np.linspace(3500.0, 5000.0, 800)
    # Continuum ~ wave/5000, plus a huge ghost spike near 3650 A.
    flux = wave / 5000.0
    spike = (wave > 3620) & (wave < 3680)
    flux = flux.copy()
    flux[spike] = 5.0

    path = tmp_path / 'toy.coadd1d.fits'
    _write_fake_coadd(path, wave, flux)

    result = clean_goodman_coadd_blue_end(path, inplace=True)
    assert result['n_flagged'] > 10

    with fits.open(path) as hdul:
        out = hdul['SPECTRUM'].data
        assert np.nanmax(out['flux'][spike]) < 1.5
        assert np.all(out['mask'][spike] == 0)
        assert hdul[0].header['PASSTABC'] is True


def test_clean_repairs_zero_wavelengths(tmp_path):
    wave = np.linspace(4000.0, 6000.0, 200)
    flux = np.ones_like(wave)
    wave = wave.copy()
    wave[50] = 0.0
    wave[51] = 0.0
    path = tmp_path / 'badwave.coadd1d.fits'
    _write_fake_coadd(path, wave, flux)

    result = clean_goodman_coadd_blue_end(path, inplace=True)
    assert result['n_wave_fixed'] == 2
    with fits.open(path) as hdul:
        out_w = np.asarray(hdul['SPECTRUM'].data['wave'], float)
        assert np.all(out_w > 0)
        assert np.min(out_w) > 3900

