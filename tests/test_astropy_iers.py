"""Tests for Astropy IERS initialization used during PypeIt heliocentric correction."""

from __future__ import annotations

import numpy as np
from astropy import units as u
from astropy.table import QTable

from passta.astropy_iers import _scalar_mjd_value, ensure_astropy_iers_ready


def test_scalar_mjd_value_accepts_quantity_and_float():
    assert _scalar_mjd_value(61210.0 * u.day) == 61210.0
    assert _scalar_mjd_value(np.float64(61210.0)) == 61210.0


def test_iers_combine_patch_handles_scalar_mjd_column():
    from astropy.utils.iers import IERS_A

    table = QTable()
    table['UT1_UTC_A'] = np.array([0.0, 0.0])
    table['PolPMFlag_A'] = np.array(['I', 'P'])
    table['UT1_UTC_B'] = np.array([np.nan, np.nan])
    table['UT1Flag_A'] = np.array(['I', 'P'])
    table['PM_X_B'] = np.array([np.nan, np.nan])
    table['PM_Y_B'] = np.array([np.nan, np.nan])
    table['PM_x_A'] = np.array([0.0, 0.0])
    table['PM_y_A'] = np.array([0.0, 0.0])
    table['dX_2000A_A'] = np.array([0.0, 0.0])
    table['dY_2000A_A'] = np.array([0.0, 0.0])
    table['dX_2000A_B'] = np.array([np.nan, np.nan])
    table['dY_2000A_B'] = np.array([np.nan, np.nan])
    table['NutFlag_A'] = np.array(['I', 'P'])
    table['MJD'] = np.array([61209.0, 61210.0])  # plain float, not Quantity

    out = IERS_A._combine_a_b_columns(table)
    assert out.meta['predictive_mjd'] == 61210.0


def test_ensure_astropy_iers_ready_loads_table():
    ensure_astropy_iers_ready()
    from astropy.utils import iers

    table = iers.earth_orientation_table.get()
    assert len(table) > 0
    assert 'predictive_mjd' in table.meta


def test_geomotion_velocity_after_iers_ready():
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time
    from pypeit.core.wave import geomotion_velocity

    ensure_astropy_iers_ready()
    loc = EarthLocation.from_geodetic(-70.8068 * u.deg, -30.2403 * u.deg, height=2743 * u.m)
    obstime = Time('2026-06-24T00:07:28.548', location=loc)
    coord = SkyCoord(141.84 * u.deg, -32.01 * u.deg, frame='icrs')
    vel = geomotion_velocity(obstime, coord, frame='heliocentric')
    assert np.isfinite(vel)
