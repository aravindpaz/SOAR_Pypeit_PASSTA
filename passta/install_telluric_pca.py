"""Download PypeIt telluric PCA from Google Drive and register it with ``pypeit_install_telluric``."""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

logger = logging.getLogger(__name__)

# https://drive.google.com/file/d/1DvIGYiFJG9oWkYDKRP8x8ronD6fkW9_c/view?usp=drive_link
DEFAULT_GDRIVE_FILE_ID = '1DvIGYiFJG9oWkYDKRP8x8ronD6fkW9_c'
DEFAULT_OUT_NAME = 'TellPCA_3000_26000_R15000.fits'
# Reject HTML / virus-scan stub pages mistaken for the FITS.
MIN_VALID_BYTES = 100_000


def _gdown_download_uc(file_id: str, dest: Path, *, quiet: bool) -> None:
    import gdown

    url = f'https://drive.google.com/uc?id={file_id}'
    dest.parent.mkdir(parents=True, exist_ok=True)
    gdown.download(url, str(dest), quiet=quiet)


def download_telluric_pca_from_gdrive(
    dest_path: str | Path,
    *,
    file_id: str = DEFAULT_GDRIVE_FILE_ID,
    quiet: bool = False,
) -> Path:
    """Download ``TellPCA_3000_26000_R15000.fits`` from the given Google Drive file id.

    Parameters
    ----------
    dest_path
        Full path for the downloaded FITS (parent directories are created).
    file_id
        Google Drive document id from the share link.
    quiet
        Passed through to ``gdown.download``.

    Returns
    -------
    pathlib.Path
        Path to the downloaded file.
    """
    out = Path(dest_path).expanduser().resolve()
    logger.info('Downloading telluric PCA (id=%s) -> %s', file_id, out)
    _gdown_download_uc(file_id, out, quiet=quiet)
    if not out.is_file():
        raise FileNotFoundError(f'Download did not create {out}')
    size = out.stat().st_size
    if size < MIN_VALID_BYTES:
        out.unlink(missing_ok=True)
        raise RuntimeError(
            f'Downloaded file is only {size} bytes; expected a large FITS. '
            'Check Drive permissions or network (login / quota pages save as tiny HTML).'
        )
    logger.info('Downloaded %s (%d bytes)', out.name, size)
    return out


def install_telluric_pca_via_pypeit(local_fits: str | Path) -> None:
    """Run ``pypeit_install_telluric --local_file`` so PypeIt caches under ``telluric/atm_grids``."""
    fits = Path(local_fits).resolve()
    if not fits.is_file():
        raise FileNotFoundError(fits)

    exe = shutil.which('pypeit_install_telluric')
    if exe is None:
        raise RuntimeError(
            'pypeit_install_telluric not found on PATH. Activate the conda env where PypeIt is installed.'
        )

    cmd = [exe, '--local_file', str(fits)]
    logger.info('Running: %s', ' '.join(cmd))
    subprocess.run(cmd, check=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            'Download TellPCA_3000_26000_R15000.fits from Google Drive and install it with '
            'pypeit_install_telluric --local_file (avoids PypeIt S3 when DNS blocks s3-west.nrp-nautilus.io).'
        ),
    )
    parser.add_argument(
        '--file-id',
        default=DEFAULT_GDRIVE_FILE_ID,
        help='Google Drive file id (default: PypeIt TellPCA grid from the official share link).',
    )
    parser.add_argument(
        '--workdir',
        type=Path,
        default=None,
        help='Directory for the temporary download (default: system temp directory).',
    )
    parser.add_argument(
        '--keep-download',
        action='store_true',
        help='Do not delete the downloaded FITS after installing to the PypeIt cache.',
    )
    parser.add_argument(
        '-q',
        '--quiet',
        action='store_true',
        help='Less gdown console output.',
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        help='Log at DEBUG.',
    )
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(levelname)s %(name)s: %(message)s',
    )

    work = args.workdir
    if work is not None:
        work = work.expanduser().resolve()
        work.mkdir(parents=True, exist_ok=True)

    try:
        if work is not None:
            dest = work / DEFAULT_OUT_NAME
            download_telluric_pca_from_gdrive(dest, file_id=args.file_id, quiet=args.quiet)
            install_telluric_pca_via_pypeit(dest)
            if not args.keep_download:
                dest.unlink(missing_ok=True)
                logger.info('Removed download copy: %s', dest)
        else:
            with tempfile.TemporaryDirectory(prefix='passta_tellpca_') as td:
                dest = Path(td) / DEFAULT_OUT_NAME
                download_telluric_pca_from_gdrive(dest, file_id=args.file_id, quiet=args.quiet)
                install_telluric_pca_via_pypeit(dest)
    except subprocess.CalledProcessError as exc:
        logger.error('pypeit_install_telluric failed with exit code %s', exc.returncode)
        return exc.returncode or 1
    except Exception:
        logger.exception('Telluric PCA install failed')
        return 1

    logger.info('Done. TellPCA grid should be in the PypeIt cache (telluric/atm_grids).')
    return 0


if __name__ == '__main__':
    sys.exit(main())
