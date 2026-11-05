"""Tests for LCO archive RLEVEL handling in ``goodman-spec-download``."""

from __future__ import annotations

from passta.lcogt import (
    LCO_RLEVEL_CFZST,
    LCO_RLEVEL_ECFZST,
    LCO_RLEVEL_LEGACY_REDUCED,
    LCO_RLEVEL_RAW,
    LCO_RLEVEL_SPECTROSCOPY,
    LCO_RLEVEL_WECFZST,
)
from passta.spec_download import (
    _filter_bad_basenames,
    _is_pipeline_reduced_frame,
)


def test_lco_rlevel_spectroscopy_includes_legacy_and_2026_codes():
    assert LCO_RLEVEL_RAW in LCO_RLEVEL_SPECTROSCOPY
    for level in LCO_RLEVEL_LEGACY_REDUCED:
        assert level in LCO_RLEVEL_SPECTROSCOPY
    assert LCO_RLEVEL_CFZST in LCO_RLEVEL_SPECTROSCOPY
    assert LCO_RLEVEL_ECFZST in LCO_RLEVEL_SPECTROSCOPY
    assert LCO_RLEVEL_WECFZST in LCO_RLEVEL_SPECTROSCOPY


def test_is_pipeline_reduced_frame_detects_all_pipeline_prefixes():
    for prefix in ('cfzst', 'ecfzst', 'wecfzst'):
        assert _is_pipeline_reduced_frame({'basename': f'{prefix}_20260418_001'})
        assert _is_pipeline_reduced_frame({'filename': f'{prefix}_20260418_001.fits.fz'})


def test_is_pipeline_reduced_frame_accepts_raw_frames():
    raw = {'basename': '20260418_001', 'filename': '20260418_001.fits.fz'}
    assert not _is_pipeline_reduced_frame(raw)


def test_filter_bad_basenames_keeps_raw_only():
    frames = [
        {'id': 1, 'basename': '20260418_001'},
        {'id': 2, 'basename': 'cfzst_20260418_002'},
        {'id': 3, 'basename': 'ecfzst_20260418_003'},
        {'id': 4, 'basename': 'wecfzst_20260418_004'},
    ]
    kept = _filter_bad_basenames(frames)
    assert [f['id'] for f in kept] == [1]
