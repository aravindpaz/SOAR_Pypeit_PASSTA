"""Batch PypeIt reduction for every PASSTA target under ``OUTDIR/ToOs/``."""

from __future__ import annotations

import argparse
import logging
import os
import sys
from typing import Sequence

from passta.logging_config import add_logging_arguments, configure_from_parsed_args
from passta.pypeit_reduce import (
    default_caldb_dir,
    discover_to_os_reduction_jobs,
    run_reduction,
)
from passta.versioning import add_version_argument

logger = logging.getLogger(__name__)


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse CLI for :func:`main`."""
    parser = argparse.ArgumentParser(
        description=(
            'Run the same PASSTA Goodman / PypeIt reduction as goodman-pypeit-reduce for '
            'each target under OUTDIR/ToOs/ (flat ``ToOs/<name>/rawdata`` from gdrive, '
            'or legacy ``ToOs/<name>/<YYYYMMDD>/rawdata``).'
        ),
    )
    add_version_argument(parser)
    parser.add_argument(
        'outdir',
        type=str,
        help='Base directory containing ``ToOs/`` (same layout as goodman-spec-download gdrive).',
    )
    parser.add_argument(
        '--frame-extension',
        default='.fits',
        help=(
            'Suffix passed to ``pypeit_setup`` for frames linked into each object directory '
            '(default ``.fits`` for uncompressed Drive data; try ``.fits.fz`` for fpacked-only).'
        ),
    )
    parser.add_argument(
        '--snid',
        default=False,
        action='store_true',
        help='Run SNID SAGE classification on telluric-corrected spectra.',
    )
    parser.add_argument(
        '--slack',
        default=False,
        action='store_true',
        help='Post SNID results to Slack (requires ``--snid`` and ``SLACK_TOKEN_FILE_PASSTA``).',
    )
    parser.add_argument(
        '--slack-dry-run',
        dest='slack_dry_run',
        default=False,
        action='store_true',
        help=(
            'Build the same Slack message as ``--slack`` (write ``*.slack.debug`` under '
            '``spectra/slack/``) and log it; run SNID; do not call the Slack API or read tokens.'
        ),
    )
    parser.add_argument(
        '--caldb',
        type=str,
        default=None,
        help='Directory with caldb.txt and sensitivity functions (default: packaged caldb).',
    )
    parser.add_argument(
        '--jobs',
        type=int,
        default=1,
        metavar='N',
        help='Parallel PypeIt setup/reduce/SNID workers (default 1). Use 2--4 on multi-core hosts.',
    )
    parser.add_argument(
        '--slit-pos',
        type=float,
        default=505.0,
        metavar='PIX',
        dest='slit_pos',
        help=(
            'Spatial pixel position (spat_pixpos) for choosing the science trace in '
            'parse_pypeit_output (default 505). Pixel tolerance is still widened in '
            'steps of 15, 30, then 60.'
        ),
    )
    parser.add_argument(
        '--no-try-old',
        dest='try_old',
        action='store_false',
        help=(
            'Enforce PypeIt ``SpecObj`` datamodel version checks during ``pypeit_flux_calib`` '
            '(default: pass ``--try_old``, skipping version checks).'
        ),
    )
    parser.set_defaults(try_old=True)
    add_logging_arguments(parser)
    args = parser.parse_args(argv)
    args.outdir = os.path.abspath(os.path.expanduser(args.outdir))
    if not str(args.frame_extension).startswith('.'):
        parser.error('--frame-extension must start with a dot (e.g. .fits)')
    if args.slack_dry_run:
        args.snid = True
    return args


def main(argv: Sequence[str] | None = None) -> None:
    """Console entry point for ``goodman-pypeit-reduce-toos``."""
    args = parse_arguments(argv)
    configure_from_parsed_args(args)
    caldb_dir = args.caldb if args.caldb is not None else default_caldb_dir()

    jobs = discover_to_os_reduction_jobs(args.outdir)
    if not jobs:
        logger.info(
            'No ToOs targets with rawdata or workspace under %s',
            os.path.join(args.outdir, 'ToOs'),
        )
        return

    for date_label, od in jobs:
        logger.info('ToOs reduction: date=%s outdir=%s', date_label, od)
        run_reduction(
            date_label,
            od,
            caldb_dir=caldb_dir,
            snid=args.snid,
            slack=args.slack,
            slack_dry_run=args.slack_dry_run,
            parallel_jobs=args.jobs,
            frame_extension=args.frame_extension,
            slit_pos=args.slit_pos,
            try_old=args.try_old,
        )


if __name__ == '__main__':
    main()
