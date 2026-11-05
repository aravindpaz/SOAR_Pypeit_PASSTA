"""Package version resolution (setuptools-scm file, then distribution metadata)."""

from __future__ import annotations

import argparse
import importlib.metadata


def get_version() -> str:
    """Return the installed or development package version string.

    Tries ``passta._version`` (from setuptools-scm on install/build), then
    :func:`importlib.metadata.version` for the ``passta`` distribution, then
    ``0.0.0+unknown`` (e.g. a raw checkout without metadata).
    """
    try:
        from passta._version import __version__ as v

        return v
    except ImportError:
        pass
    try:
        return importlib.metadata.version('passta')
    except importlib.metadata.PackageNotFoundError:
        return '0.0.0+unknown'


def add_version_argument(parser: argparse.ArgumentParser) -> None:
    """Register ``--version`` on a CLI parser (uses :func:`get_version`)."""
    parser.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s {get_version()}',
    )
