"""Tests for Goodman slit-aware reduction directory grouping."""

from __future__ import annotations

import os

from astropy.io import fits

from passta.pypeit_reduce import (
    _build_reduction_dirs,
    _goodman_slit_from_header,
    _link_reduction_directory,
    _science_object_slit_groups,
)


def test_goodman_slit_from_header_prefers_slit():
    hdr = fits.Header()
    hdr['SLIT'] = '1.0_LONG_SLIT'
    hdr['DECKER'] = 'other'
    assert _goodman_slit_from_header(hdr) == '1.0_LONG_SLIT'


def test_goodman_slit_from_header_falls_back_to_decker():
    hdr = fits.Header()
    hdr['DECKER'] = '1.5 LONG SLIT'
    assert _goodman_slit_from_header(hdr) == '1.5_LONG_SLIT'


def test_science_object_slit_groups():
    wavedata = [
        {'object': 'hr5501', 'slit': '1.0_LONG_SLIT', 'obstype': 'SPECTRUM'},
        {'object': 'hr5501', 'slit': '1.5_LONG_SLIT', 'obstype': 'SPECTRUM'},
        {'object': '2026pub', 'slit': '1.0_LONG_SLIT', 'obstype': 'SPECTRUM'},
        {'object': 'dflat', 'slit': '1.0_LONG_SLIT', 'obstype': 'LAMPFLAT'},
    ]
    assert _science_object_slit_groups(wavedata) == [
        ('2026pub', '1.0_LONG_SLIT'),
        ('hr5501', '1.0_LONG_SLIT'),
        ('hr5501', '1.5_LONG_SLIT'),
    ]


def test_link_reduction_directory_filters_by_slit(tmp_path):
    raw = tmp_path / 'raw'
    raw.mkdir()
    flat10 = raw / 'dflat.10.fits.fz'
    flat15 = raw / 'dflat.15.fits.fz'
    arc10 = raw / 'arc.10.fits.fz'
    sci10 = raw / 'hr5501.10.fits.fz'
    sci15 = raw / 'hr5501.15.fits.fz'
    for p in (flat10, flat15, arc10, sci10, sci15):
        p.write_bytes(b'')

    wavedata = [
        {
            'filename': str(flat10),
            'object': 'dflat',
            'slit': '1.0_LONG_SLIT',
            'obstype': 'LAMPFLAT',
            'mjd': 60000.0,
        },
        {
            'filename': str(flat15),
            'object': 'dflat',
            'slit': '1.5_LONG_SLIT',
            'obstype': 'LAMPFLAT',
            'mjd': 60000.1,
        },
        {
            'filename': str(arc10),
            'object': 'arc',
            'slit': '1.0_LONG_SLIT',
            'obstype': 'ARC',
            'mjd': 60000.2,
        },
        {
            'filename': str(sci10),
            'object': 'hr5501',
            'slit': '1.0_LONG_SLIT',
            'obstype': 'SPECTRUM',
            'mjd': 60000.3,
        },
        {
            'filename': str(sci15),
            'object': 'hr5501',
            'slit': '1.5_LONG_SLIT',
            'obstype': 'SPECTRUM',
            'mjd': 60000.4,
        },
    ]

    objdir10 = tmp_path / 'work' / '400_M1' / 'hr5501' / '1.0_LONG_SLIT'
    objdir15 = tmp_path / 'work' / '400_M1' / 'hr5501' / '1.5_LONG_SLIT'
    _link_reduction_directory(wavedata, str(objdir10), 'hr5501', '1.0_LONG_SLIT')
    _link_reduction_directory(wavedata, str(objdir15), 'hr5501', '1.5_LONG_SLIT')

    linked10 = {os.path.basename(os.readlink(p)) for p in objdir10.iterdir()}
    linked15 = {os.path.basename(os.readlink(p)) for p in objdir15.iterdir()}

    assert linked10 == {'dflat.10.fits.fz', 'hr5501.10.fits.fz', 'arc.10.fits.fz'}
    assert linked15 == {'dflat.15.fits.fz', 'hr5501.15.fits.fz'}
    assert 'arc.10.fits.fz' not in linked15


def test_build_reduction_dirs_layout(tmp_path):
    raw = tmp_path / 'raw'
    raw.mkdir()
    paths = {}
    for name in ('sci.fits.fz', 'flat.fits.fz'):
        p = raw / name
        p.write_bytes(b'')
        paths[name] = str(p)

    wavedata = [
        {
            'filename': paths['flat.fits.fz'],
            'object': 'dflat',
            'slit': '1.0_LONG_SLIT',
            'obstype': 'LAMPFLAT',
            'mjd': 1.0,
        },
        {
            'filename': paths['sci.fits.fz'],
            'object': 'target',
            'slit': '1.0_LONG_SLIT',
            'obstype': 'SPECTRUM',
            'mjd': 2.0,
        },
    ]
    work = tmp_path / 'workspace'
    dirs = _build_reduction_dirs(str(work), {'400_M1': wavedata})
    assert dirs == [str(work / '400_M1' / 'target' / '1.0_LONG_SLIT')]
    assert (work / '400_M1' / 'target' / '1.0_LONG_SLIT' / 'sci.fits.fz').is_symlink()
