"""In-process wrappers for PypeIt script entry points (replaces shelling out to CLIs)."""

from __future__ import annotations

import argparse
import logging

from passta.pypeit_msgs_bridge import pypeit_logging_bridge

logger = logging.getLogger(__name__)


def pypeit_verbosity() -> int:
    """Map the root logger level to PypeIt's 0 / 1 / 2 verbosity scheme."""
    lvl = logging.getLogger().getEffectiveLevel()
    if lvl <= logging.DEBUG:
        return 2
    if lvl <= logging.INFO:
        return 1
    return 0


def run_setup(
    spectrograph: str,
    root: str,
    output_path: str,
) -> None:
    """Run :class:`pypeit.scripts.setup.Setup` (``pypeit_setup``).

    Mirrors ``pypeit_setup -s <spec> -c all -r <root> -d <output_path>`` with
    defaults aligned to PASSTA's prior shell invocation. PypeIt ``msgs`` are
    routed through :func:`passta.pypeit_msgs_bridge.pypeit_logging_bridge`.
    """
    from pypeit.scripts.setup import Setup

    args = argparse.Namespace(
        spectrograph=spectrograph,
        root=[root],
        extension='.fits',
        output_path=output_path,
        overwrite=False,
        cfg_split='all',
        background=False,
        manual_extraction=False,
        keep_bad_frames=False,
        gui=False,
        version_override=None,
        date_override=None,
        verbosity=pypeit_verbosity(),
    )
    logger.info('PypeIt setup: spectrograph=%s root=%s output=%s', spectrograph, root, output_path)
    with pypeit_logging_bridge(logger):
        Setup.main(args)


def run_pypeit(
    pypeit_file: str,
    redux_path: str,
) -> None:
    """Run :class:`pypeit.scripts.run_pypeit.RunPypeIt` (``run_pypeit``)."""
    from pypeit.scripts.run_pypeit import RunPypeIt

    args = argparse.Namespace(
        pypeit_file=pypeit_file,
        verbosity=pypeit_verbosity(),
        redux_path=redux_path,
        reuse_calibs=True,
        show=False,
        overwrite=False,
        calib_only=False,
    )
    logger.info('PypeIt reduce: file=%s redux_path=%s', pypeit_file, redux_path)
    with pypeit_logging_bridge(logger):
        RunPypeIt.main(args)


def run_sensfunc(
    spec1d_fits: str,
    outfile: str,
    par_outfile: str,
) -> None:
    """Run :class:`pypeit.scripts.sensfunc.SensFunc` (``pypeit_sensfunc``)."""
    from pypeit.scripts.sensfunc import SensFunc

    args = argparse.Namespace(
        spec1dfile=spec1d_fits,
        outfile=outfile,
        par_outfile=par_outfile,
        algorithm=None,
        multi=None,
        sens_file=None,
        flatfile=None,
        debug=False,
        verbosity=pypeit_verbosity(),
    )
    logger.info('PypeIt sensfunc: in=%s out=%s', spec1d_fits, outfile)
    with pypeit_logging_bridge(logger):
        SensFunc.main(args)


def run_flux_calib(
    flux_file: str,
    par_outfile: str,
    *,
    try_old: bool = False,
) -> None:
    """Run :class:`pypeit.scripts.flux_calib.FluxCalib` (``pypeit_flux_calib``)."""
    from pypeit.scripts.flux_calib import FluxCalib

    args = argparse.Namespace(
        flux_file=flux_file,
        par_outfile=par_outfile,
        try_old=try_old,
        verbosity=pypeit_verbosity(),
    )
    logger.info('PypeIt flux_calib: flux_file=%s par_out=%s', flux_file, par_outfile)
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
    )
    logger.info('PypeIt coadd_1dspec: %s', coadd1d_file)
    with pypeit_logging_bridge(logger):
        CoAdd1DSpec.main(args)


def run_tellfit(
    spec1d_or_coadd: str,
    *,
    objmodel: str = 'poly',
    par_outfile: str = 'telluric.par',
    chk_version: bool = False,
) -> None:
    """Run :class:`pypeit.scripts.tellfit.TellFit` (``pypeit_tellfit``)."""
    from pypeit.scripts.tellfit import TellFit

    args = argparse.Namespace(
        spec1dfile=spec1d_or_coadd,
        objmodel=objmodel,
        redshift=None,
        tell_grid=None,
        pca_file=None,
        tell_file=None,
        debug=False,
        plot=False,
        par_outfile=par_outfile,
        chk_version=chk_version,
        verbosity=pypeit_verbosity(),
    )
    logger.info('PypeIt tellfit: %s objmodel=%s', spec1d_or_coadd, objmodel)
    with pypeit_logging_bridge(logger):
        TellFit.main(args)
