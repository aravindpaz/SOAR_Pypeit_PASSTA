"""Run telluric correction using PypeIt with the PASSTA-bundled atmospheric grid."""

from __future__ import annotations

import argparse
import logging
from collections.abc import Sequence
from pathlib import Path

from passta.pypeit_invoke import goodman_telluric_param_file, pypeit_verbosity
from passta.pypeit_msgs_bridge import pypeit_logging_bridge
from passta.tellfit.grid import telluric_grid_path

logger = logging.getLogger(__name__)


def _resolve_tell_grid(tell_grid: str | Path | None) -> str:
    if tell_grid is not None:
        path = Path(tell_grid).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f'Telluric grid not found: {path}')
        return str(path)
    return str(telluric_grid_path())


def run_tellfit(
    spec1d_or_coadd: str | Path,
    *,
    objmodel: str | None = None,
    redshift: float | None = None,
    tell_grid: str | Path | None = None,
    pca_file: str | None = None,
    tell_file: str | Path | None = None,
    par_outfile: str = 'telluric.par',
    debug: bool = False,
    plot: bool = False,
    chk_version: bool = False,
) -> None:
    """Telluric-correct a spec1d or coadd1d FITS file (PypeIt ``TellFit`` API).

    Unless ``tell_grid`` is set, uses the PCA grid packaged in ``passta/data/telluric/``.
    """
    from pypeit.scripts.tellfit import TellFit

    grid_path = _resolve_tell_grid(tell_grid)
    spec_path = str(Path(spec1d_or_coadd).expanduser().resolve())
    tell_path = None
    if tell_file is not None:
        tell_path = str(Path(tell_file).expanduser().resolve())

    args = argparse.Namespace(
        spec1dfile=spec_path,
        objmodel=objmodel,
        redshift=redshift,
        tell_grid=grid_path,
        pca_file=pca_file,
        tell_file=tell_path,
        debug=debug,
        plot=plot,
        par_outfile=par_outfile,
        chk_version=chk_version,
        verbosity=pypeit_verbosity(),
        log_file=None,
        log_level=None,
    )
    logger.info('Telluric grid: %s', grid_path)
    if tell_path:
        logger.info('Telluric parameter file: %s', tell_path)
    with pypeit_logging_bridge(logger):
        TellFit.main(args)


def default_goodman_tell_file() -> str | None:
    """Packaged Goodman telluric par block, if present."""
    return goodman_telluric_param_file()


def build_cli_parser():
    """Argument parser matching ``pypeit_tellfit`` plus a PASSTA ``--version`` flag."""
    from pypeit.scripts.tellfit import TellFit

    from passta.versioning import add_version_argument

    parser = TellFit.get_parser()
    add_version_argument(parser)
    parser.description = (
        'Telluric-correct a coadd1d or spec1d spectrum using PypeIt telluric fitting '
        'and the PASSTA-bundled atmospheric PCA grid.'
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point for ``goodman-pypeit-tellfit``."""
    parser = build_cli_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)

    if args.tell_grid is None:
        args.tell_grid = str(telluric_grid_path())

    try:
        from pypeit.scripts.tellfit import TellFit

        with pypeit_logging_bridge(logger):
            TellFit.main(args)
    except Exception:
        logger.exception('Telluric fitting failed')
        return 1
    return 0
