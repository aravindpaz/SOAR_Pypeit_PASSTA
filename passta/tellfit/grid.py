"""Packaged PypeIt telluric PCA atmospheric grid."""

from __future__ import annotations

from importlib import resources
from pathlib import Path

DEFAULT_TELLURIC_GRID_NAME = 'TellPCA_3000_26000_R15000.fits'


def telluric_grid_path(name: str = DEFAULT_TELLURIC_GRID_NAME) -> Path:
    """Return the absolute path to a packaged telluric grid FITS file.

    The grid is shipped under ``passta/data/telluric/`` in the same float64
    multi-HDU format used by PypeIt (``read_telluric_pca``).

    Raises
    ------
    FileNotFoundError
        If the requested grid is not present in the installed package.
    """
    try:
        pkg_path = resources.files('passta').joinpath('data/telluric', name)
        if pkg_path.is_file():
            return Path(pkg_path)
    except (ModuleNotFoundError, TypeError, FileNotFoundError, OSError):
        pass

    fallback = Path(__file__).resolve().parents[1] / 'data' / 'telluric' / name
    if fallback.is_file():
        return fallback

    raise FileNotFoundError(
        f'Telluric grid {name!r} not found in passta package data '
        f'(expected passta/data/telluric/{name}).'
    )
