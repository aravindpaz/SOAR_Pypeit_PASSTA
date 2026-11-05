"""PASSTA telluric fitting (PypeIt-compatible, bundled atmospheric grid)."""

from passta.tellfit.core import run_tellfit
from passta.tellfit.grid import DEFAULT_TELLURIC_GRID_NAME, telluric_grid_path

# Pipeline integration reads this flag; standalone ``goodman-pypeit-tellfit`` is unaffected.
TELLURIC_CORRECTION_ENABLED = False

__all__ = [
    'DEFAULT_TELLURIC_GRID_NAME',
    'TELLURIC_CORRECTION_ENABLED',
    'run_tellfit',
    'telluric_grid_path',
]
