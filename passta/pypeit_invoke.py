"""In-process wrappers for PypeIt script entry points (replaces shelling out to CLIs)."""

from __future__ import annotations

import os

# Must be set before matplotlib is imported (parallel reductions + QA plots).
os.environ.setdefault('MPLBACKEND', 'Agg')

import argparse
import contextlib
import copy
import logging
import shlex
from collections.abc import Callable, Iterator
from importlib import resources
from pathlib import Path
from typing import Any

from passta.pypeit_msgs_bridge import pypeit_logging_bridge
from passta import astropy_iers as _astropy_iers  # noqa: F401 - IERS patch
from passta import pypeit_files as _pypeit_files  # noqa: F401 - InputFile read patch
from passta.pypeit_files import sanitize_pypeit_file, sanitize_pypeit_files

logger = logging.getLogger(__name__)

_GOODMAN_SPEC_PREFIX = 'soar_goodman'
_pypeit_standard_utf8_patch_applied = False
_arc_fit_qa_patch_applied = False
_headless_matplotlib_applied = False


def ensure_headless_matplotlib() -> None:
    """Use the non-interactive Agg backend (safe with ``--jobs`` > 1)."""
    global _headless_matplotlib_applied
    os.environ.setdefault('MPLBACKEND', 'Agg')
    try:
        import matplotlib

        if matplotlib.get_backend().lower() != 'agg':
            matplotlib.use('Agg', force=True)
    except Exception:
        pass
    _headless_matplotlib_applied = True


def _ensure_pypeit_arc_fit_qa_patch() -> None:
    """Skip arc wavelength QA plots when matplotlib mathtext/layout fails."""
    global _arc_fit_qa_patch_applied
    if _arc_fit_qa_patch_applied:
        return
    import functools

    from pypeit.core.wavecal import autoid

    _orig_arc_fit_qa = autoid.arc_fit_qa

    @functools.wraps(_orig_arc_fit_qa)
    def _safe_arc_fit_qa(*args: Any, **kwargs: Any) -> None:
        ensure_headless_matplotlib()
        try:
            _orig_arc_fit_qa(*args, **kwargs)
        except ValueError:
            logger.warning(
                'PypeIt arc_fit_qa failed (mathtext/layout); skipping wavelength QA plot',
                exc_info=True,
            )
            try:
                import matplotlib.pyplot as plt

                plt.close('all')
            except Exception:
                pass

    autoid.arc_fit_qa = _safe_arc_fit_qa
    _arc_fit_qa_patch_applied = True


def _ensure_pypeit_standard_info_utf8_patch() -> None:
    """Normalize PypeIt standard catalog text files to ASCII (en-dash breaks Table.read)."""
    global _pypeit_standard_utf8_patch_applied
    if _pypeit_standard_utf8_patch_applied:
        return

    from pathlib import Path

    import pypeit

    standards_root = Path(pypeit.__file__).resolve().parent / 'data' / 'standards'
    if not standards_root.is_dir():
        _pypeit_standard_utf8_patch_applied = True
        return

    replacements = (
        (b'\xe2\x80\x93', b'-'),  # en-dash
        (b'\xe2\x80\x94', b'-'),  # em-dash
    )
    for path in standards_root.rglob('*_info.txt'):
        raw = path.read_bytes()
        fixed = raw
        for old, new in replacements:
            fixed = fixed.replace(old, new)
        if fixed != raw:
            path.write_bytes(fixed)
            logger.info('Normalized non-ASCII punctuation in %s', path)

    _pypeit_standard_utf8_patch_applied = True


ensure_headless_matplotlib()
_ensure_pypeit_standard_info_utf8_patch()
_ensure_pypeit_arc_fit_qa_patch()


def goodman_flatfield_param_file() -> str | None:
    """Path to packaged Goodman flat-field PypeIt parameter block, if present."""
    try:
        path = resources.files('passta').joinpath('config/goodman_flatfield.par')
        if path.is_file():
            return str(path)
    except (ModuleNotFoundError, TypeError, FileNotFoundError):
        pass
    fallback = Path(__file__).resolve().parent / 'config' / 'goodman_flatfield.par'
    return str(fallback) if fallback.is_file() else None


def goodman_telluric_param_file() -> str | None:
    """Path to packaged Goodman telluric PypeIt parameter block, if present."""
    try:
        path = resources.files('passta').joinpath('config/goodman_telluric.par')
        if path.is_file():
            return str(path)
    except (ModuleNotFoundError, TypeError, FileNotFoundError):
        pass
    fallback = Path(__file__).resolve().parent / 'config' / 'goodman_telluric.par'
    return str(fallback) if fallback.is_file() else None


def _goodman_detector_data_hdu_index(hdu: Any) -> int:
    """Index of the HDU that carries Goodman ``GAIN`` / science header (0 or 1).

    PypeIt's Goodman classes assume ``hdu[1]`` (typical fpacked layout). Single-HDU
    archive ``.fits`` only have ``hdu[0]``, which caused ``IndexError`` in
    ``get_detector_par`` during reduction (e.g. ``bpm`` -> ``get_rawimage``).
    """
    if hdu is None or len(hdu) == 0:
        return 0
    if len(hdu) == 1:
        return 0
    try:
        hdr1 = hdu[1].header
        if 'GAIN' in hdr1:
            return 1
    except (IndexError, KeyError, AttributeError):
        pass
    for i in range(len(hdu)):
        hdr = getattr(hdu[i], 'header', None)
        if hdr is not None and 'GAIN' in hdr:
            return i
    return 0


def _passta_goodman_red_get_detector_par(self, det: int, hdu: Any = None) -> Any:  # noqa: ANN401
    """``SOARGoodmanRedSpectrograph.get_detector_par`` with correct HDU index for single-ext .fits."""
    import numpy as np

    from pypeit import PypeItError
    from pypeit.images import detector_container

    if hdu is None:
        binning = '2,2'
        gain = None
        ronoise = None
        datasec = None
        oscansec = None
        edx = 1
    else:
        edx = _goodman_detector_data_hdu_index(hdu)
        hdr = hdu[edx].header
        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
        gain = np.atleast_1d(hdr['GAIN'])
        ronoise = np.atleast_1d(hdr['RDNOISE'])
        datasec = None
        oscansec = None

    detector_dict = dict(
        binning=binning,
        det=1,
        dataext=edx,
        specaxis=1,
        specflip=False,
        spatflip=False,
        platescale=0.15,
        darkcurr=0.0,
        saturation=65535.0,
        nonlinear=1.0,
        mincounts=-1e10,
        numamplifiers=1,
        gain=gain,
        ronoise=ronoise,
        datasec=datasec,
        oscansec=oscansec,
    )

    if hdu is None:
        return detector_container.DetectorContainer(**detector_dict)

    hdr = hdu[edx].header
    if binning != '2,2':
        raise PypeItError('Ask the developers to add your binning.  Or add it yourself.')
    col0 = int(hdr['TRIMSEC'][1:].split(':')[0])
    detector_dict['datasec'] = np.atleast_1d(f'[:,{col0 * 2}:]')
    detector_dict['oscansec'] = np.atleast_1d(f'[:,1:{int(col0 * 2) - 2}:]')
    return detector_container.DetectorContainer(**detector_dict)


def _passta_goodman_blue_get_detector_par(self, det: int, hdu: Any = None) -> Any:  # noqa: ANN401
    """``SOARGoodmanBlueSpectrograph.get_detector_par`` with correct HDU index for single-ext .fits."""
    import numpy as np

    from pypeit.images import detector_container

    if hdu is None:
        binning = '2,2'
        gain = None
        ronoise = None
        datasec = None
        oscansec = None
        edx = 1
    else:
        edx = _goodman_detector_data_hdu_index(hdu)
        hdr = hdu[edx].header
        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
        gain = np.atleast_1d(hdr['GAIN'])
        ronoise = np.atleast_1d(hdr['RDNOISE'])
        datasec = None
        oscansec = None

    detector_dict = dict(
        binning=binning,
        det=1,
        dataext=edx,
        specaxis=1,
        specflip=False,
        spatflip=False,
        platescale=0.15,
        darkcurr=0.0,
        saturation=65535.0,
        nonlinear=1.0,
        mincounts=-1e10,
        numamplifiers=1,
        gain=gain,
        ronoise=ronoise,
        datasec=np.asarray(['[:,20:4112]']),
        oscansec=np.asarray(['[:,2:16]']),
    )

    if hdu is None:
        return detector_container.DetectorContainer(**detector_dict)
    return detector_container.DetectorContainer(**detector_dict)


@contextlib.contextmanager
def goodman_raw_fits_runtime(spectrograph: str, extension: str | None) -> Iterator[None]:
    """Compatibility for SOAR/Goodman raw ``.fits`` for setup **and** reduction.

    Upstream :class:`~pypeit.spectrographs.soar_goodman.SOARGoodmanSpectrograph`
    sets ``allowed_extensions = ['.fz']`` only. That breaks ``pypeit_setup`` with
    ``-e .fits`` and later :class:`~pypeit.pypeit.PypeIt` (``_check_extensions`` on
    raw paths listed in the ``.pypeit`` file). PASSTA widens ``allowed_extensions``
    for ``.fits`` / ``None`` for the whole time this context is active.

    Uncompressed Goodman frames are often single-HDU while ``meta`` uses
    ``ext=1`` (matching fpacked ``.fz``). We duplicate the primary header in
    :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.get_headarr` when only
    one extension is present so setup metadata matches PypeIt expectations.

    Reduction assumes ``hdu[1]`` in ``get_detector_par``; PASSTA patches Red/Blue
    ``get_detector_par`` to use the HDU that actually contains ``GAIN`` and sets
    ``dataext`` accordingly so ``get_rawimage`` / BPM construction work on
    single-extension ``.fits``.
    """
    if not spectrograph.startswith(_GOODMAN_SPEC_PREFIX):
        yield
        return
    if extension not in ('.fits', None):
        yield
        return
    from pypeit.spectrographs.soar_goodman import (
        SOARGoodmanBlueSpectrograph,
        SOARGoodmanRedSpectrograph,
        SOARGoodmanSpectrograph,
    )

    orig_ext = list(SOARGoodmanSpectrograph.allowed_extensions)
    orig_get_headarr: Callable[..., Any] = SOARGoodmanSpectrograph.get_headarr
    orig_red_gdp: Callable[..., Any] = SOARGoodmanRedSpectrograph.get_detector_par
    orig_blue_gdp: Callable[..., Any] = SOARGoodmanBlueSpectrograph.get_detector_par

    def _patched_get_headarr(self, inp, strict=True):  # noqa: ANN001
        headers = orig_get_headarr(self, inp, strict=strict)
        if headers is not None and len(headers) == 1:
            return [headers[0], copy.deepcopy(headers[0])]
        return headers

    patched = False
    if '.fits' not in orig_ext:
        SOARGoodmanSpectrograph.allowed_extensions = list(dict.fromkeys([*orig_ext, '.fits']))
        patched = True
    SOARGoodmanSpectrograph.get_headarr = _patched_get_headarr  # type: ignore[method-assign]
    SOARGoodmanRedSpectrograph.get_detector_par = _passta_goodman_red_get_detector_par  # type: ignore[method-assign]
    SOARGoodmanBlueSpectrograph.get_detector_par = _passta_goodman_blue_get_detector_par  # type: ignore[method-assign]
    try:
        yield
    finally:
        SOARGoodmanSpectrograph.get_headarr = orig_get_headarr  # type: ignore[method-assign]
        SOARGoodmanRedSpectrograph.get_detector_par = orig_red_gdp  # type: ignore[method-assign]
        SOARGoodmanBlueSpectrograph.get_detector_par = orig_blue_gdp  # type: ignore[method-assign]
        if patched:
            SOARGoodmanSpectrograph.allowed_extensions = orig_ext


def pypeit_verbosity() -> int:
    """Map the root logger level to PypeIt's 0 / 1 / 2 verbosity scheme."""
    lvl = logging.getLogger().getEffectiveLevel()
    if lvl <= logging.DEBUG:
        return 2
    if lvl <= logging.INFO:
        return 1
    return 0


def _log_pypeit_cli(name: str, *argv: str) -> None:
    """Log the shell-equivalent PypeIt command (in-process invocation uses the same args)."""
    cmd = ' '.join([name, *(shlex.quote(a) for a in argv)])
    logger.info('PypeIt equivalent command: %s', cmd)


def run_setup(
    spectrograph: str,
    root: str,
    output_path: str,
    *,
    extension: str = '.fits',
) -> None:
    """Run :class:`pypeit.scripts.setup.Setup` (``pypeit_setup``).

    Mirrors ``pypeit_setup -s <spec> -c all -r <root> -d <output_path>`` with
    defaults aligned to PASSTA's prior shell invocation. PypeIt ``msgs`` are
    routed through :func:`passta.pypeit_msgs_bridge.pypeit_logging_bridge`.

    Parameters
    ----------
    extension : str, optional
        Filename suffix for frames under ``root`` (e.g. ``.fits`` or ``.fz``).
        For SOAR/Goodman raw ``.fits``, use :func:`goodman_raw_fits_runtime` in
        :func:`passta.pypeit_reduce.run_reduction` so setup and ``run_pypeit`` share
        the same compatibility window.
    """
    ensure_headless_matplotlib()
    from pypeit.scripts.setup import Setup

    param_block_file = None
    if spectrograph.startswith(_GOODMAN_SPEC_PREFIX):
        param_block_file = goodman_flatfield_param_file()
        if param_block_file:
            logger.info('PypeIt flat-field param block %s', param_block_file)

    args = argparse.Namespace(
        spectrograph=spectrograph,
        root=[root],
        extension=extension,
        output_path=output_path,
        overwrite=False,
        cfg_split='all',
        background=False,
        manual_extraction=False,
        flexure=False,
        keep_bad_frames=False,
        gui=False,
        version_override=None,
        date_override=None,
        param_block_file=param_block_file,
        verbosity=pypeit_verbosity(),
        log_file=None,
        log_level=None,
    )
    setup_argv = [
        '-s',
        spectrograph,
        '-c',
        'all',
        '-r',
        root,
        '-d',
        output_path,
        '-e',
        extension,
    ]
    if param_block_file:
        setup_argv.extend(['--param_block_file', param_block_file])
    _log_pypeit_cli('pypeit_setup', *setup_argv)
    with pypeit_logging_bridge(logger):
        Setup.main(args)
    n_fixed = sanitize_pypeit_files(sorted(Path(output_path).rglob('*.pypeit')))
    if n_fixed:
        logger.info('Normalized %s .pypeit file(s) to ASCII after setup', n_fixed)


def run_pypeit(
    pypeit_file: str,
    redux_path: str,
    *,
    overwrite: bool = False,
    reuse_calibs: bool = True,
) -> None:
    """Run :class:`pypeit.scripts.run_pypeit.RunPypeIt` (``run_pypeit``).

    Parameters
    ----------
    overwrite
        If True, re-reduce frames even when Science products already exist.
        PASSTA enables this automatically when some science frames lack
        ``spec1d`` products (incomplete prior reductions).
    reuse_calibs
        Passed through to PypeIt (default True).
    """
    ensure_headless_matplotlib()
    sanitize_pypeit_file(pypeit_file)
    from pypeit.scripts.run_pypeit import RunPypeIt

    args = argparse.Namespace(
        pypeit_file=pypeit_file,
        verbosity=pypeit_verbosity(),
        redux_path=redux_path,
        reuse_calibs=reuse_calibs,
        show=False,
        overwrite=overwrite,
        calib_only=False,
        log_file=None,
        log_level=None,
    )
    argv = [pypeit_file, '-r', redux_path]
    if overwrite:
        argv.append('-o')
    _log_pypeit_cli('run_pypeit', *argv)
    with pypeit_logging_bridge(logger):
        RunPypeIt.main(args)


def run_sensfunc(
    spec1d_fits: str,
    outfile: str,
    par_outfile: str,
) -> None:
    """Run :class:`pypeit.scripts.sensfunc.SensFunc` (``pypeit_sensfunc``)."""
    from pypeit.scripts.sensfunc import SensFunc

    # Match ``pypeit.scripts.sensfunc.SensFunc.get_parser`` (API drift: ``use_flat``, ``extr``, ``spec1dfiles``).
    args = argparse.Namespace(
        spec1dfiles=[spec1d_fits],
        outfile=outfile,
        par_outfile=par_outfile,
        algorithm=None,
        multi=None,
        sens_file=None,
        use_flat=False,
        extr=None,
        debug=False,
        verbosity=pypeit_verbosity(),
        log_file=None,
        log_level=None,
    )
    sens_argv = [spec1d_fits, '-o', outfile, '--par_outfile', par_outfile]
    if args.use_flat:
        sens_argv.append('--use_flat')
    _log_pypeit_cli('pypeit_sensfunc', *sens_argv)
    with pypeit_logging_bridge(logger):
        SensFunc.main(args)


def run_flux_calib(
    flux_file: str,
    par_outfile: str,
    *,
    try_old: bool = True,
) -> None:
    """Run :class:`pypeit.scripts.flux_calib.FluxCalib` (``pypeit_flux_calib``)."""
    from pypeit.scripts.flux_calib import FluxCalib

    args = argparse.Namespace(
        flux_file=flux_file,
        par_outfile=par_outfile,
        try_old=try_old,
        verbosity=pypeit_verbosity(),
        log_file=None,
        log_level=None,
    )
    flux_argv = [flux_file, '--par_outfile', par_outfile]
    if try_old:
        flux_argv.append('--try_old')
    _log_pypeit_cli('pypeit_flux_calib', *flux_argv)
    with pypeit_logging_bridge(logger):
        FluxCalib.main(args)


def run_coadd_1dspec(
    coadd1d_file: str,
    par_outfile: str,
) -> None:
    """Run :class:`pypeit.scripts.coadd_1dspec.CoAdd1DSpec` (``pypeit_coadd_1dspec``)."""
    from pypeit.scripts.coadd_1dspec import CoAdd1DSpec

    args = argparse.Namespace(
        coadd1d_file=coadd1d_file,
        par_outfile=par_outfile,
        debug=False,
        show=False,
        verbosity=pypeit_verbosity(),
        log_file=None,
        log_level=None,
    )
    _log_pypeit_cli('pypeit_coadd_1dspec', coadd1d_file, '--par_outfile', par_outfile)
    with pypeit_logging_bridge(logger):
        CoAdd1DSpec.main(args)


def run_tellfit(
    spec1d_or_coadd: str,
    *,
    objmodel: str = 'poly',
    par_outfile: str = 'telluric.par',
    tell_file: str | None = None,
    chk_version: bool = False,
) -> None:
    """Run PASSTA telluric fitting (:func:`passta.tellfit.run_tellfit`)."""
    from passta.tellfit import run_tellfit as _run_tellfit

    tell_argv = [spec1d_or_coadd, '--objmodel', objmodel, '--par_outfile', par_outfile]
    if tell_file:
        tell_argv.extend(['--tell_file', tell_file])
    if chk_version:
        tell_argv.append('--chk_version')
    _log_pypeit_cli('goodman-pypeit-tellfit', *tell_argv)
    _run_tellfit(
        spec1d_or_coadd,
        objmodel=objmodel,
        par_outfile=par_outfile,
        tell_file=tell_file,
        chk_version=chk_version,
    )
