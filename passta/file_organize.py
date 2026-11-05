"""Prepare local Goodman FITS for PASSTA layout (rawdata, workspace, fpack)."""

from __future__ import annotations

import argparse
import glob
import logging
import os
import subprocess
from typing import Sequence

from passta.logging_config import add_logging_arguments, configure_from_parsed_args
from passta.versioning import add_version_argument

from astropy.io import fits
from astropy.time import Time

logger = logging.getLogger(__name__)


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments for :func:`organize_files`.

    Parameters
    ----------
    argv : sequence of str or None
        Arguments to parse; defaults to ``sys.argv[1:]`` when None.

    Returns
    -------
    argparse.Namespace
        Parsed ``date``, ``input_dir``, and ``outdir``.
    """
    parser = argparse.ArgumentParser(
        description='Organize Goodman FITS into dated rawdata tree and fpack.',
    )
    add_version_argument(parser)
    parser.add_argument(
        'date',
        type=str,
        help='Night/date label (parsed by Astropy Time, e.g. YYYY-MM-DD).',
    )
    parser.add_argument(
        'input_dir',
        type=str,
        help='Directory containing input ``.fits`` / ``.fz`` files.',
    )
    parser.add_argument(
        '--outdir',
        type=str,
        default='.',
        help='Base output directory (default: current directory).',
    )
    add_logging_arguments(parser)
    return parser.parse_args(argv)


def organize_files(date: str, input_dir: str, outdir: str = '.') -> None:
    """Copy FITS into ``<outdir>/<YYYYMMDD>/rawdata``, fix headers, and fpack.

    Object frames have ``OBSTYPE`` rewritten to ``SPECTRUM`` when the input
    value is ``OBJECT`` (case-insensitive). Unpacked ``.fits`` are removed after
    successful ``fpack`` when a ``.fz`` is produced.

    Parameters
    ----------
    date : str
        Date string understood by :class:`astropy.time.Time`.
    input_dir : str
        Source directory of FITS data.
    outdir : str, optional
        Base directory for outputs (default ``.``).

    Raises
    ------
    subprocess.CalledProcessError
        If ``fpack`` exits non-zero.
    """
    datedir = Time(date).datetime.strftime('%Y%m%d')
    fulloutdir = os.path.join(outdir, datedir, 'rawdata')
    fullworkdir = os.path.join(outdir, datedir, 'workspace')

    os.makedirs(fulloutdir, exist_ok=True)
    os.makedirs(fullworkdir, exist_ok=True)

    input_files = glob.glob(os.path.join(input_dir, '*.fits'))
    input_files += glob.glob(os.path.join(input_dir, '*.fz'))

    for path in input_files:
        with fits.open(path) as hdu:
            obstype = hdu[0].header['OBSTYPE']
            if obstype.lower() == 'object':
                hdu[0].header['OBSTYPE'] = 'SPECTRUM'

            outfile = os.path.join(fulloutdir, os.path.basename(path))
            hdu.writeto(outfile, overwrite=True, output_verify='silentfix')

        if outfile.endswith('.fits'):
            packfile = outfile + '.fz'
            if not os.path.exists(packfile):
                cmd = ['fpack', outfile]
                logger.info('Running: %s', ' '.join(cmd))
                subprocess.run(cmd, check=True)
            os.remove(outfile)


def main(argv: Sequence[str] | None = None) -> None:
    """Console entry point."""
    args = parse_arguments(argv)
    configure_from_parsed_args(args)
    organize_files(args.date, args.input_dir, outdir=args.outdir)


if __name__ == '__main__':
    main()
