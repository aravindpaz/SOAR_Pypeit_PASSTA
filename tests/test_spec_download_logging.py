"""Tests for goodman-spec-download progress logging."""

from __future__ import annotations

import logging
from unittest.mock import MagicMock, patch

import pytest

from passta.logging_config import configure_logging
from passta.spec_download import download_night


@pytest.fixture(autouse=True)
def _logging_configured():
    configure_logging(logging.INFO, force=True)


def test_download_night_warns_on_empty_propid(caplog: pytest.LogCaptureFixture):
    with caplog.at_level(logging.INFO):
        download_night('2026-05-12', '  ,  ', outdir='/tmp/passta-test-out')
    assert 'No proposal IDs' in caplog.text


@patch('passta.spec_download.LCOGT')
def test_download_night_logs_when_no_science(mock_lco_class, caplog: pytest.LogCaptureFixture, tmp_path):
    lco = mock_lco_class.return_value
    lco.get_username_password.return_value = ('user@example.com', 'secret')
    lco.get_token_header.return_value = {'Authorization': 'Token abc'}
    lco.params = {'date_format': '%Y-%m-%dT%H:%M:%S'}
    lco.get_obslist.return_value = []
    lco.download_obslist.return_value = None

    with caplog.at_level(logging.INFO):
        download_night('2026-05-12', 'SOAR2025B-004', outdir=str(tmp_path))

    text = caplog.text
    assert 'Starting LCO download' in text
    assert 'Authenticated to LCO archive' in text
    assert 'Querying LCO archive' in text
    assert '[SOAR2025B-004] science SPECTRUM' in text
    assert 'No Goodman spectrograph science spectra found' in text
    lco.download_obslist.assert_not_called()


@patch('passta.spec_download.LCOGT')
def test_download_night_exits_on_missing_credentials(mock_lco_class):
    lco = mock_lco_class.return_value
    lco.get_username_password.return_value = (None, None)

    with pytest.raises(SystemExit) as exc:
        download_night('2026-05-12', 'SOAR2025B-004', outdir='/tmp/x')
    assert exc.value.code == 1
