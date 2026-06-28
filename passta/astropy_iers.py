"""Astropy IERS helpers for threaded PypeIt reductions."""

from __future__ import annotations

import functools
import logging
import threading
import urllib.error

import astropy.units as u
import numpy as np

logger = logging.getLogger(__name__)

_iers_patch_applied = False
_iers_init_lock = threading.Lock()
_iers_ready = False


def _scalar_mjd_value(mjd: object) -> float:
    value = mjd.value if hasattr(mjd, 'value') else mjd
    return float(value)


def _iers_a_with_day_mjd(iers_a):
    """Ensure the IERS-A ``MJD`` column has day units for Astropy merges."""
    mjd_col = iers_a['MJD']
    if hasattr(mjd_col, 'unit') and mjd_col.unit is not None:
        return iers_a
    iers_a = iers_a.copy()
    mjd_val = mjd_col.value if hasattr(mjd_col, 'value') else mjd_col
    iers_a['MJD'] = np.asarray(mjd_val) * u.day
    return iers_a


def _ensure_iers_combine_patch() -> None:
    """Work around astropy 7.2.0 assuming ``MJD[i]`` is always a Quantity."""
    global _iers_patch_applied
    if _iers_patch_applied:
        return

    from astropy.utils.iers import IERS_A

    _orig_combine = IERS_A._combine_a_b_columns.__func__

    @classmethod
    @functools.wraps(_orig_combine)
    def _combine_a_b_columns_patched(cls, iers_a):
        iers_a = _iers_a_with_day_mjd(iers_a)
        try:
            return _orig_combine(cls, iers_a)
        except AttributeError as exc:
            if 'value' not in str(exc):
                raise
            table = iers_a[np.isfinite(iers_a['UT1_UTC_A']) & (iers_a['PolPMFlag_A'] != '')]
            table = cls._substitute_iers_b(table)
            p_index = min(
                np.searchsorted(table['UT1Flag_A'], 'P'),
                np.searchsorted(table['PolPMFlag_A'], 'P'),
            )
            table.meta['predictive_index'] = p_index
            table.meta['predictive_mjd'] = _scalar_mjd_value(table['MJD'][p_index])
            return table

    IERS_A._combine_a_b_columns = _combine_a_b_columns_patched
    _iers_patch_applied = True


def _load_iers_table():
    from astropy.utils import iers

    iers.IERS_Auto.iers_table = None
    table = iers.IERS_Auto.open()
    if isinstance(table, iers.IERS):
        return table

    logger.warning(
        'IERS_Auto.open() returned %s; retrying with bundled table only.',
        type(table).__name__,
    )
    iers.IERS_Auto.iers_table = None
    iers.conf.auto_download = False
    table = iers.IERS_Auto.open()
    if not isinstance(table, iers.IERS):
        raise ValueError(
            f'IERS_Auto.open() returned {type(table).__name__}, expected IERS table'
        )
    return table


def ensure_astropy_iers_ready() -> None:
    """Load IERS tables once before parallel reductions use heliocentric correction."""
    global _iers_ready
    if _iers_ready:
        return

    _ensure_iers_combine_patch()

    from astropy.utils import iers

    with _iers_init_lock:
        if _iers_ready:
            return
        tries = 0
        while tries < 5:
            try:
                table = _load_iers_table()
                iers.conf.auto_download = False
                iers.earth_orientation_table.get()
                _iers_ready = True
                logger.debug(
                    'Astropy IERS table ready (%s rows, predictive_mjd=%s).',
                    len(table),
                    table.meta.get('predictive_mjd'),
                )
                return
            except (urllib.error.URLError, OSError, ValueError, AttributeError) as exc:
                tries += 1
                logger.warning(
                    'Astropy IERS initialization failed (attempt %s/5): %s',
                    tries,
                    exc,
                )
                iers.IERS_Auto.iers_table = None
        raise RuntimeError('Failed to initialize Astropy IERS tables after 5 attempts')


_ensure_iers_combine_patch()
