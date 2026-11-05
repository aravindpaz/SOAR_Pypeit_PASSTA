"""Tests for Goodman-only filtering in ``goodman-spec-download``."""

from __future__ import annotations

from passta.spec_download import (
    _filter_goodman_spectrograph_frames,
    _is_goodman_spectrograph_frame,
)


def test_is_goodman_spectrograph_frame_accepts_ghts_red_blue():
    assert _is_goodman_spectrograph_frame({'INSTRUME': 'GHTS_RED'})
    assert _is_goodman_spectrograph_frame({'INSTRUME': 'ghts_blue'})


def test_is_goodman_spectrograph_frame_rejects_tspec_and_imager():
    assert not _is_goodman_spectrograph_frame({'INSTRUME': 'triplespec'})
    assert not _is_goodman_spectrograph_frame({'INSTRUME': 'TSpec'})
    assert not _is_goodman_spectrograph_frame({'INSTRUME': 'ghts_red_imager'})
    assert not _is_goodman_spectrograph_frame({'INSTRUME': 'GHTS_RED_IMAGER'})
    assert not _is_goodman_spectrograph_frame({'INSTRUME': 'ghts_blue_imager'})


def test_filter_goodman_spectrograph_frames():
    frames = [
        {'id': 1, 'INSTRUME': 'ghts_red'},
        {'id': 2, 'INSTRUME': 'triplespec'},
        {'id': 3, 'INSTRUME': 'ghts_red_imager'},
    ]
    kept = _filter_goodman_spectrograph_frames(frames)
    assert [f['id'] for f in kept] == [1]
