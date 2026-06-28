"""Tests for SNID SAGE ``.output`` parsing used by ``run_snid_sage`` / ``--snid``."""

from __future__ import annotations

import textwrap

from passta.pypeit_reduce import (
    _parse_snid_sage_output_legacy,
    _parse_snid_sage_output_pypi,
    _parse_snid_sage_summary,
)


SAMPLE_PYPY = textwrap.dedent(
    """\
    ================================
    SNID-SAGE CLASSIFICATION RESULTS
    ================================

    ---------------------------------------------------------------------------
    | Type   | Subtype      |        z |      z_err |  age [d] |  age_err [d] |
    ---------------------------------------------------------------------------
    | II     | II-87A       | 0.005938 |   0.000872 |    190.6 |        278.1 |

    TEMPLATE MATCHES (from Best Cluster - Auto-Selected):
    ------------------------------------------------------------------------------
      # Template         Type   Subtype   HσLAP-CCC    Redshift      ±Error    Age
    ------------------------------------------------------------------------------
      1 PTF10qts         Ic     Ic-broad       2.67    0.062096    0.005672   -1.2
      2 SN2007gr         Ic     Ic-narrow      2.10    0.061000    0.006000    0.0
    """
)


def test_parse_snid_sage_summary():
    summary = _parse_snid_sage_summary(SAMPLE_PYPY)
    assert summary is not None
    assert summary['type'] == 'II'
    assert summary['subtype'] == 'II-87A'
    assert summary['redshift'] == '0.005938'


def test_parse_snid_sage_output_pypi_table():
    tbl = _parse_snid_sage_output_pypi(SAMPLE_PYPY)
    assert tbl is not None
    assert len(tbl) == 2
    assert tbl['name'][0] == 'PTF10qts'
    assert tbl['type'][0] == 'Ic'
    assert tbl['subtype'][0] == 'Ic-broad'
    assert tbl['redshift'][0] == '0.062096'
    assert tbl['age'][1] == '0.0'


def test_parse_snid_sage_output_legacy_table():
    legacy = textwrap.dedent(
        """\
        header
        --------------------
        1  PTF10qts  Ic  Ic-broad  2.67  0.062  0.005  -1.2
        2  SN2007gr  Ic  Ic-narrow  2.1  0.061  0.006  0.0
        """
    )
    tbl = _parse_snid_sage_output_legacy(legacy)
    assert tbl is not None
    assert len(tbl) == 2
    assert tbl['name'][0] == 'PTF10qts'


def test_goodman_reduce_toos_uses_same_run_snid_sage():
    """``goodman-pypeit-reduce-toos`` calls ``run_reduction`` → ``run_snid_sage`` (same as reduce)."""
    import inspect

    from passta import pypeit_reduce
    from passta import pypeit_reduce_toos

    src_reduce = inspect.getsource(pypeit_reduce.run_reduction)
    assert 'run_snid_sage' in src_reduce
    src_toos = inspect.getsource(pypeit_reduce_toos.main)
    assert 'run_reduction' in src_toos
    assert 'snid=args.snid' in src_toos
