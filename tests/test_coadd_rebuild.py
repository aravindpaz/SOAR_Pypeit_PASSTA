"""Regression tests for recovering incomplete Science products into the coadd."""

from __future__ import annotations

from pathlib import Path

from passta.pypeit_reduce import _coadd_list_matches_objdata


def _write_coadd1d(path: Path, rows: list[tuple[str, str]]) -> None:
    lines = [
        '[coadd1d]',
        '  coaddfile = /tmp/out.coadd1d.fits',
        '  wave_method = linear',
        '',
        'coadd1d read',
        '  path /tmp',
        'filename | obj_id',
    ]
    for fitsname, obj_id in rows:
        lines.append(f'{fitsname} | {obj_id}')
    lines.append('coadd1d end')
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def test_coadd_list_matches_objdata_same_rows(tmp_path: Path):
    rows = [
        ('/tmp/spec1d_a.fits', 'SPAT0518-SLIT0493-DET01'),
        ('/tmp/spec1d_b.fits', 'SPAT0518-SLIT0493-DET01'),
        ('/tmp/spec1d_c.fits', 'SPAT0518-SLIT0493-DET01'),
    ]
    cfg = tmp_path / 'target.coadd1d'
    _write_coadd1d(cfg, rows)
    # objdata tuples: (spec1d_txt, objname, ...)
    objdata = [(r[0].replace('.fits', '.txt'), r[1], 't', 0.0, '400_M1', '1.0') for r in rows]
    assert _coadd_list_matches_objdata(str(cfg), objdata)


def test_coadd_list_matches_objdata_detects_missing_exposure(tmp_path: Path):
    cfg = tmp_path / 'target.coadd1d'
    _write_coadd1d(
        cfg,
        [
            ('/tmp/spec1d_a.fits', 'SPAT0518-SLIT0493-DET01'),
            ('/tmp/spec1d_b.fits', 'SPAT0518-SLIT0493-DET01'),
        ],
    )
    objdata = [
        ('/tmp/spec1d_a.txt', 'SPAT0518-SLIT0493-DET01', 't', 0.0, '400_M1', '1.0'),
        ('/tmp/spec1d_b.txt', 'SPAT0518-SLIT0493-DET01', 't', 0.0, '400_M1', '1.0'),
        ('/tmp/spec1d_c.txt', 'SPAT0518-SLIT0493-DET01', 't', 0.0, '400_M1', '1.0'),
    ]
    assert not _coadd_list_matches_objdata(str(cfg), objdata)
