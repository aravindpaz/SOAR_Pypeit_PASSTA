"""Tests for Goodman WAVMODE embedding in download filenames."""

from __future__ import annotations

from astropy.io import fits

from passta.lcogt import LCOGT


def test_goodman_wavmode_from_fits_reads_hdu1(tmp_path):
    path = tmp_path / 'test.fits.fz'
    phdu = fits.PrimaryHDU()
    shdu = fits.ImageHDU()
    shdu.header['WAVMODE'] = '400_M2'
    fits.HDUList([phdu, shdu]).writeto(path, overwrite=True)

    assert LCOGT._goodman_wavmode_from_fits(str(path)) == '400_M2'
