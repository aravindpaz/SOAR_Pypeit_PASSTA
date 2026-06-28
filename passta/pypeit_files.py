"""Validate and normalize PypeIt input files for PASSTA."""

from __future__ import annotations

import logging
import os
from importlib import resources
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)

_UNICODE_REPLACEMENTS = str.maketrans(
    {
        '\u2013': '-',  # en-dash
        '\u2014': '-',  # em-dash
        '\u2212': '-',  # minus sign
        '\u00a0': ' ',  # nbsp
        '\u2018': "'",
        '\u2019': "'",
        '\u201c': '"',
        '\u201d': '"',
    }
)

_pypeit_inputfiles_patch_applied = False


def normalize_pypeit_text(text: str) -> str:
    """Replace common Unicode punctuation with ASCII for PypeIt input files."""
    return text.translate(_UNICODE_REPLACEMENTS)


def non_ascii_byte_offsets(raw: bytes) -> list[int]:
    return [index for index, byte in enumerate(raw) if byte > 127]


def validate_ascii_file(path: Path) -> list[str]:
    """Return human-readable errors when ``path`` is not strict ASCII."""
    raw = path.read_bytes()
    bad = non_ascii_byte_offsets(raw)
    if not bad:
        return []
    messages = [
        f'{path}: non-ASCII byte 0x{raw[index]:02x} at offset {index}'
        for index in bad[:5]
    ]
    if len(bad) > 5:
        messages.append(f'{path}: ... and {len(bad) - 5} more non-ASCII bytes')
    return messages


def packaged_config_paths() -> list[Path]:
    """Paths to packaged ``passta/config/*.par`` blocks merged into ``.pypeit`` files."""
    paths: list[Path] = []
    try:
        root = resources.files('passta').joinpath('config')
        for name in ('goodman_flatfield.par', 'goodman_telluric.par'):
            candidate = root.joinpath(name)
            if candidate.is_file():
                paths.append(Path(str(candidate)))
    except (ModuleNotFoundError, TypeError, FileNotFoundError):
        fallback = Path(__file__).resolve().parent / 'config'
        paths.extend(sorted(fallback.glob('*.par')))
    return paths


def packaged_pypeit_fixture_paths() -> list[Path]:
    """Paths to checked-in ``.pypeit`` fixtures under ``tests/fixtures``."""
    fixtures = Path(__file__).resolve().parent.parent / 'tests' / 'fixtures'
    if not fixtures.is_dir():
        return []
    return sorted(fixtures.glob('*.pypeit'))


def sanitize_pypeit_file(path: str | Path) -> bool:
    """Rewrite a ``.pypeit`` file to ASCII when needed.

    Returns True when the file on disk was modified.
    """
    p = Path(path)
    raw = p.read_bytes()
    try:
        raw.decode('ascii')
        return False
    except UnicodeDecodeError:
        text = normalize_pypeit_text(raw.decode('utf-8'))
        p.write_bytes(text.encode('ascii'))
        logger.warning('Normalized non-ASCII characters in %s for PypeIt', p)
        return True


def sanitize_pypeit_files(paths: list[str | Path]) -> int:
    """Sanitize each path; return the number of files rewritten."""
    return sum(1 for path in paths if sanitize_pypeit_file(path))


def read_pypeit_lines(path: str | Path) -> np.ndarray:
    """Read a PypeIt input file the way PASSTA expects PypeIt to ingest it."""
    lines = []
    with open(path, 'r', encoding='utf-8') as handle:
        for line in handle:
            lines.append(normalize_pypeit_text(line.replace('\t', ' ').rstrip()))
    return np.array(lines)


def validate_pypeit_ingest(path: str | Path) -> None:
    """Raise when PypeIt cannot parse ``path`` as a reduction file."""
    from pypeit import inputfiles

    inputfiles.PypeItFile.from_file(str(path), vet=False)


def ensure_pypeit_inputfiles_read_patch() -> None:
    """Patch PypeIt ``InputFile.readlines`` for UTF-8 and ASCII normalization."""
    global _pypeit_inputfiles_patch_applied
    if _pypeit_inputfiles_patch_applied:
        return

    from pypeit import inputfiles
    from pypeit import PypeItError

    _orig_readlines = inputfiles.InputFile.readlines

    def _readlines_patched(ifile: str) -> np.ndarray:
        if not os.path.isfile(ifile):
            raise PypeItError('The filename does not exist -\n' + ifile)
        return read_pypeit_lines(ifile)

    inputfiles.InputFile.readlines = staticmethod(_readlines_patched)  # type: ignore[method-assign]
    _pypeit_inputfiles_patch_applied = True


def validate_packaged_pypeit_inputs() -> list[str]:
    """Validate packaged config and fixture files used by the reduction pipeline."""
    errors: list[str] = []
    for path in packaged_config_paths():
        errors.extend(validate_ascii_file(path))
    for path in packaged_pypeit_fixture_paths():
        errors.extend(validate_ascii_file(path))
        try:
            validate_pypeit_ingest(path)
        except Exception as exc:
            errors.append(f'{path}: PypeIt ingest failed: {exc}')
    return errors


def main() -> None:
    """Console entry point: ``passta-check-pypeit-files``."""
    errors = validate_packaged_pypeit_inputs()
    if errors:
        for message in errors:
            print(message)
        raise SystemExit(1)
    print('OK: packaged PypeIt config and fixtures are ASCII and ingestible.')


ensure_pypeit_inputfiles_read_patch()
