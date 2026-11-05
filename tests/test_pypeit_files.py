"""Tests that PASSTA-generated PypeIt input files remain ASCII and ingestible."""

from __future__ import annotations

from pathlib import Path

import pytest

from passta.pypeit_files import (
    normalize_pypeit_text,
    packaged_config_paths,
    packaged_pypeit_fixture_paths,
    read_pypeit_lines,
    sanitize_pypeit_file,
    validate_ascii_file,
    validate_packaged_pypeit_inputs,
    validate_pypeit_ingest,
)


@pytest.fixture
def flatfield_param_text() -> str:
    path = Path(__file__).resolve().parents[1] / 'passta' / 'config' / 'goodman_flatfield.par'
    return path.read_text(encoding='utf-8')


def test_packaged_pypeit_config_files_are_ascii():
    errors: list[str] = []
    for path in packaged_config_paths():
        errors.extend(validate_ascii_file(path))
    assert not errors, '\n'.join(errors)


def test_packaged_pypeit_inputs_validate():
    errors = validate_packaged_pypeit_inputs()
    assert not errors, '\n'.join(errors)


@pytest.mark.parametrize('fixture_path', packaged_pypeit_fixture_paths(), ids=lambda p: p.name)
def test_fixture_pypeit_files_are_ascii_and_ingestible(fixture_path: Path):
    assert not validate_ascii_file(fixture_path)
    validate_pypeit_ingest(fixture_path)


def test_goodman_flatfield_param_block_uses_ascii_dashes(flatfield_param_text: str):
    assert '\u2013' not in flatfield_param_text
    assert '3580-3720' in flatfield_param_text


def test_normalize_pypeit_text_replaces_en_dash():
    assert normalize_pypeit_text('3580\u20133720 A') == '3580-3720 A'


def test_sanitize_pypeit_file_rewrites_unicode(tmp_path: Path):
    path = tmp_path / 'bad.pypeit'
    path.write_bytes(
        b'# Auto-generated\n[rdx]\n    spectrograph = soar_goodman_red\n'
        b'# Blue-end ghost mitigation (3580\xe2\x80\x933720 A)\n'
        b'setup read\nSetup A:\n  binning: 2,2\nsetup end\n'
        b'data read\n path /tmp\n filename | frametype\n'
        b'a.fits | arc\n'
        b'data end\n'
    )
    assert sanitize_pypeit_file(path)
    assert not validate_ascii_file(path)
    validate_pypeit_ingest(path)


def test_read_pypeit_lines_normalizes_unicode(tmp_path: Path):
    path = tmp_path / 'unicode.pypeit'
    path.write_text('# 3580\u20133720\n', encoding='utf-8')
    lines = read_pypeit_lines(path)
    assert lines[0] == '# 3580-3720'


def test_pypeit_inputfiles_patch_reads_unicode_without_ascii_decode_error(tmp_path: Path):
    from pypeit import inputfiles

    path = tmp_path / 'unicode.pypeit'
    path.write_bytes(
        b'# Auto-generated\n[rdx]\n    spectrograph = soar_goodman_red\n'
        b'# ghost (3580\xe2\x80\x933720 A)\n'
        b'setup read\nSetup A:\n  binning: 2,2\nsetup end\n'
        b'data read\n path /tmp\n filename | frametype\n'
        b'a.fits | arc\n'
        b'data end\n'
    )
    table = inputfiles.PypeItFile.from_file(str(path), vet=False)
    assert table.setup is not None
    assert table.data is not None


def test_param_block_merged_into_pypeit_stays_ascii(tmp_path: Path, flatfield_param_text: str):
    """Regression: param-block comments must not break PypeIt ASCII reads."""
    pypeit_path = tmp_path / 'merged.pypeit'
    pypeit_path.write_text(
        '# Auto-generated\n[rdx]\n    spectrograph = soar_goodman_red\n'
        + '\n'.join(flatfield_param_text.strip().splitlines()[2:])
        + '\n'
        'setup read\nSetup A:\n  binning: 2,2\nsetup end\n'
        'data read\n path /tmp\n filename | frametype\n'
        'a.fits | arc\n'
        'data end\n',
        encoding='utf-8',
    )
    assert not validate_ascii_file(pypeit_path)
    validate_pypeit_ingest(pypeit_path)
