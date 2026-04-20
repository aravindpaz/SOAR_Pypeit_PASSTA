"""CLI logging setup for PASSTA console scripts."""

from __future__ import annotations

import argparse
import logging
import sys
from typing import TextIO


def add_logging_arguments(parser: argparse.ArgumentParser) -> None:
    """Register ``--log-level`` and ``-v`` / ``--verbose`` on an ``ArgumentParser``."""
    parser.add_argument(
        '--log-level',
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        help='Minimum message level for the root logger (default: INFO).',
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='count',
        default=0,
        help='Increase verbosity (-v sets DEBUG and overrides --log-level).',
    )


def log_level_from_args(verbose: int, log_level: str) -> int:
    """Map CLI ``verbose`` count and ``log_level`` string to a ``logging`` level."""
    if verbose >= 1:
        return logging.DEBUG
    return getattr(logging, str(log_level).upper(), logging.INFO)


def configure_from_parsed_args(args: argparse.Namespace) -> None:
    """Apply :func:`configure_logging` using ``verbose`` and ``log_level`` on ``args``."""
    configure_logging(log_level_from_args(args.verbose, args.log_level))


def configure_logging(
    level: int = logging.INFO,
    *,
    stream: TextIO | None = None,
    force: bool = True,
) -> None:
    """Configure the root logger (timestamps, level, stderr).

    Safe to call from each CLI entry point; uses ``force=True`` so repeated
    invocations in tests replace the previous handler configuration.

    Parameters
    ----------
    level : int
        A ``logging`` level constant (e.g. ``logging.INFO``).
    stream : text stream, optional
        Log sink; defaults to ``sys.stderr``.
    force : bool, optional
        Passed through to :func:`logging.basicConfig` (Python 3.8+).
    """
    fmt = '%(asctime)s %(levelname)s [%(name)s] %(message)s'
    datefmt = '%Y-%m-%dT%H:%M:%S'
    logging.basicConfig(
        level=level,
        format=fmt,
        datefmt=datefmt,
        stream=stream or sys.stderr,
        force=force,
    )
