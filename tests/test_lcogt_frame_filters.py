"""Tests for LCO archive row filtering in :class:`passta.lcogt.LCOGT`."""

from __future__ import annotations

from passta.lcogt import LCOGT, LCO_RLEVEL_SPECTROSCOPY


def _row(*, telid='4m0a', configuration_type='SPECTRUM', rlevel=0):
    return {'TELID': telid, 'configuration_type': configuration_type, 'RLEVEL': rlevel}


def test_frame_row_matches_tuple_rlevel():
    row = _row()
    assert LCOGT._frame_row_matches(row, '4m0a', 'SPECTRUM', LCO_RLEVEL_SPECTROSCOPY)
    assert LCOGT._frame_row_matches(row, '4m0a', 'SPECTRUM', (0, 62, 63))
    assert not LCOGT._frame_row_matches(row, '4m0a', 'SPECTRUM', (62, 63, 64))


def test_frame_row_matches_list_rlevel():
    row = _row(rlevel=62)
    assert LCOGT._frame_row_matches(row, '4m0a', 'SPECTRUM', list(LCO_RLEVEL_SPECTROSCOPY))
    assert not LCOGT._frame_row_matches(row, '4m0a', 'SPECTRUM', [0])


def test_frame_row_matches_telid_is_exact_not_substring():
    row = _row()
    assert LCOGT._frame_row_matches(row, '4m0a', None, None)
    assert not LCOGT._frame_row_matches(row, '4m0b', None, None)
    assert not LCOGT._frame_row_matches(_row(telid='4'), '4m0a', None, None)
