"""Tests for PypeIt frame extension inference."""

from __future__ import annotations

from passta.pypeit_reduce import _infer_frame_extension


def test_infer_frame_extension_fpacked_lco(tmp_path):
    (tmp_path / 'target.ut260512.94509043.fits.fz').write_bytes(b'')
    assert _infer_frame_extension(str(tmp_path)) == '.fz'


def test_infer_frame_extension_uncompressed_fits(tmp_path):
    (tmp_path / 'target.ut260512.94509043.fits').write_bytes(b'')
    assert _infer_frame_extension(str(tmp_path)) == '.fits'


def test_infer_frame_extension_empty_dir(tmp_path):
    assert _infer_frame_extension(str(tmp_path), preferred='.fits') == '.fits'
