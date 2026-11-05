"""Download SOAR Goodman PASSTA spectra and calibrations from the LCO archive or Google Drive."""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import glob
import inspect
import logging
import os
import re
import shutil
import sys
import tempfile
from typing import Sequence

from astropy import units as u
from astropy.time import Time, TimeDelta

from passta.logging_config import (
    add_logging_arguments,
    configure_from_parsed_args,
)
from passta.lcogt import LCOGT, LCO_RLEVEL_RAW, LCO_RLEVEL_SPECTROSCOPY
from passta.versioning import add_version_argument

logger = logging.getLogger(__name__)


def _call_gdown(fn, **kwargs):
    """Call a gdown function with only kwargs supported by the installed gdown version."""
    allowed = set(inspect.signature(fn).parameters)
    return fn(**{k: v for k, v in kwargs.items() if k in allowed})


def _fetch_gdrive_folder_manifest(url: str, output_td: str, *, use_cookies: bool):
    """Return gdown's file list for a folder (``skip_download=True``), or ``None`` on hard failure."""
    import gdown

    return _call_gdown(
        gdown.download_folder,
        url=url,
        output=output_td,
        skip_download=True,
        quiet=True,
        use_cookies=use_cookies,
        remaining_ok=True,
    )


_GDOWN_FAQ = 'https://github.com/wkentaro/gdown?tab=readme-ov-file#faq'
# gdown always loads Google session cookies from this path when use_cookies=True.
_GDOWN_COOKIE_CACHE = os.path.join(os.path.expanduser('~'), '.cache', 'gdown', 'cookies.txt')


def _resolve_gdrive_cookies_source(cli_path: str | None) -> str | None:
    """Prefer CLI path, then ``PASSTA_GDRIVE_COOKIES``, then ``GDRIVE_COOKIES``."""
    if cli_path and str(cli_path).strip():
        return os.path.abspath(os.path.expanduser(str(cli_path).strip()))
    for key in ('PASSTA_GDRIVE_COOKIES', 'GDRIVE_COOKIES'):
        v = os.environ.get(key)
        if v and v.strip():
            return os.path.abspath(os.path.expanduser(v.strip()))
    return None


def _install_gdrive_cookies_for_gdown(source: str) -> None:
    """Copy Netscape-format cookies into gdown's cache path (gdown has no custom-path API)."""
    if not os.path.isfile(source):
        raise SystemExit(f'goodman-spec-download gdrive: cookies file not found: {source}')
    dest_dir = os.path.dirname(_GDOWN_COOKIE_CACHE)
    os.makedirs(dest_dir, exist_ok=True)
    shutil.copy2(source, _GDOWN_COOKIE_CACHE)
    try:
        os.chmod(_GDOWN_COOKIE_CACHE, 0o600)
    except OSError:
        pass
    logger.info('Installed Google cookies for gdown at %s (from %s)', _GDOWN_COOKIE_CACHE, source)


def _is_gdown_file_url_error(exc: BaseException) -> bool:
    """True for gdown's public-link / confirmation-page failure (class name stable across versions)."""
    return exc.__class__.__name__ == 'FileURLRetrievalError' and exc.__class__.__module__.startswith(
        'gdown',
    )


def _gdrive_access_help(folder_id: str, file_id: str, root: str, detail: str) -> str:
    return (
        'Google Drive per-file download failed (listing the folder can still succeed).\n'
        f"  Target: {root!r}  file id: {file_id}\n"
        '  Typical causes: files or parent folders are not shared as '
        '"Anyone with the link" (viewer), Google rate limits after many requests, '
        'or Google returned a confirmation page gdown cannot parse.\n'
        '  Team / restricted shares: export Netscape cookies while logged into Google '
        f'(see FAQ), then pass --gdrive-cookies FILE or set PASSTA_GDRIVE_COOKIES '
        f'(installed to {_GDOWN_COOKIE_CACHE}).\n'
        f'  Folder: https://drive.google.com/drive/folders/{folder_id}\n'
        f'  Try in browser: https://drive.google.com/uc?id={file_id}\n'
        f'  gdown FAQ: {_GDOWN_FAQ}\n'
        f'  Detail: {detail.strip()}'
    )


def _night_tree_only_empty_or_incomplete(rawdir: str, workdir: str) -> bool:
    """Whether it is safe to remove a partial night tree (no completed rawdata files)."""
    if _dir_nonempty(workdir):
        return False
    if not os.path.isdir(rawdir):
        return True
    for name in os.listdir(rawdir):
        p = os.path.join(rawdir, name)
        if os.path.isfile(p) and os.path.getsize(p) > 0:
            return False
    return True


_BAD_BASE = ('cfzst', 'ecfzst', 'wecfzst')

# Goodman spectrograph modes handled by PASSTA / PypeIt (not TripleSpec/TSpec or imager).
_GOODMAN_SPEC_INSTRUME = frozenset({'ghts_red', 'ghts_blue'})


def _is_goodman_spectrograph_frame(frame: dict) -> bool:
    """True only for ``ghts_red`` / ``ghts_blue`` spectrograph rows.

  Excludes TripleSpec/TSpec, ``*_imager`` (e.g. ``ghts_red_imager``), and any other
    ``INSTRUME`` not in :data:`_GOODMAN_SPEC_INSTRUME`.
    """
    instr = str(frame.get('INSTRUME', '')).strip().lower()
    if not instr:
        return False
    if 'triple' in instr or instr.startswith('tspec'):
        return False
    if instr.endswith('_imager') or instr.endswith('imager') or 'imager' in instr:
        return False
    return instr in _GOODMAN_SPEC_INSTRUME


def _filter_goodman_spectrograph_frames(
    frames: list[dict],
    *,
    label: str = '',
) -> list[dict]:
    """Keep Goodman spectrograph frames only; log how many were skipped."""
    kept = [f for f in frames if _is_goodman_spectrograph_frame(f)]
    skipped = len(frames) - len(kept)
    if skipped:
        suffix = f' ({label})' if label else ''
        logger.info(
            'Skipped %s frame(s) not ghts_red/ghts_blue spectrograph%s '
            '(TripleSpec/TSpec, *_imager, or other INSTRUME)',
            skipped,
            suffix,
        )
    return kept


def _archive_name_has_pipeline_prefix(name: str) -> bool:
    """True if an LCO archive basename/filename is a pipeline-reduced product."""
    stem = name.lower()
    for ext in ('.fits.fz', '.fits', '.fz'):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
            break
    return stem.startswith(_BAD_BASE)


def _is_pipeline_reduced_frame(frame: dict) -> bool:
    """True if an archive row refers to cfzst/ecfzst/wecfzst data (not raw FITS)."""
    for key in ('basename', 'filename'):
        val = frame.get(key)
        if val and _archive_name_has_pipeline_prefix(str(val)):
            return True
    return False


def _dedupe_frames_by_id(frames: list[dict]) -> list[dict]:
    """Return frames in original order, dropping duplicate archive ``id`` values."""
    seen: set[int] = set()
    out = []
    for f in frames:
        fid = f['id']
        if fid not in seen:
            seen.add(fid)
            out.append(f)
    return out


def _bad_archive_name(basename: str) -> bool:
    return _archive_name_has_pipeline_prefix(basename)


def _filter_bad_basenames(frames: list[dict]) -> list[dict]:
    return [f for f in frames if not _is_pipeline_reduced_frame(f)]


def _goodman_science_frames(frames: list[dict]) -> list[dict]:
    """Pipeline-safe Goodman spectrograph science rows."""
    return _filter_goodman_spectrograph_frames(
        [f for f in frames if not _is_pipeline_reduced_frame(f)],
        label='after pipeline filter',
    )


def _rmtree_date_if_unused(fulloutdir: str, fullworkdir: str, outdir: str, datedir: str) -> bool:
    """Remove an empty night tree; return True if ``outdir/datedir`` was removed."""
    date_path = os.path.join(outdir, datedir)
    if not glob.glob(os.path.join(fulloutdir, '*')) and not glob.glob(os.path.join(fullworkdir, '*')):
        shutil.rmtree(date_path, ignore_errors=True)
        return True
    return False


def _ensure_lco_archive_auth(lco: LCOGT) -> None:
    """Verify LCO credentials and archive token; exit with a clear message on failure."""
    username, password = lco.get_username_password()
    if not username or not password:
        logger.error(
            'LCO archive credentials missing; set LCOGT_USERNAME and LCOGT_PASSWORD '
            'in the environment before running goodman-spec-download.',
        )
        raise SystemExit(1)
    headers = lco.get_token_header(username, password)
    if headers is None:
        logger.error('LCO archive authentication failed for %s', username)
        raise SystemExit(1)
    logger.info('Authenticated to LCO archive as %s', username)


def _parse_gdrive_folder_id(url_or_id: str) -> str:
    """Return a Drive folder id from a full ``folders/...`` URL or from a bare id string."""
    s = url_or_id.strip()
    m = re.search(r'/drive/(?:u/\d+/)?folders/([^/?&#]+)', s)
    if m:
        return m.group(1)
    if re.fullmatch(r'[A-Za-z0-9_-]+', s):
        return s
    raise ValueError(f'Could not parse Google Drive folder id from: {url_or_id!r}')


def _dir_nonempty(path: str) -> bool:
    return os.path.isdir(path) and bool(os.listdir(path))


def _rel_under_drive_root(gd_path: str, root: str) -> str:
    """Strip ``root/`` from a gdown path (POSIX-style)."""
    p = gd_path.replace('\\', '/')
    prefix = root + '/'
    if not p.startswith(prefix):
        raise ValueError(f'Path {gd_path!r} does not start with {root!r}/')
    return p[len(prefix) :]


def _infer_datedir_from_relpaths(relpaths: list[str]) -> str | None:
    """If every nested path shares one leading ``YYYYMMDD/`` segment, return that id."""
    dated: list[str] = []
    for rp in relpaths:
        parts = rp.replace('\\', '/').split('/')
        if len(parts) >= 2 and re.fullmatch(r'\d{8}', parts[0]):
            dated.append(parts[0])
    if not dated:
        return None
    uniq = set(dated)
    if len(uniq) > 1:
        raise ValueError(
            f'Multiple YYYYMMDD folder prefixes under the same target: {sorted(uniq)}; '
            'split the Drive or use a narrower --folder-url.'
        )
    return dated[0]


def _calendar_ymd_tokens_from_strings(paths: list[str]) -> list[str]:
    """Parse calendar dates embedded in path strings into ``YYYYMMDD`` (validated)."""
    out: list[str] = []
    for raw in paths:
        s = raw.replace('\\', '/')
        for m in re.finditer(r'(?<!\d)(\d{2})-(\d{2})-(\d{4})(?!\d)', s):
            d, mo, y = int(m.group(1)), int(m.group(2)), int(m.group(3))
            if 1 <= mo <= 12 and 1 <= d <= 31:
                try:
                    datetime(y, mo, d)
                    out.append(f'{y:04d}{mo:02d}{d:02d}')
                except ValueError:
                    pass
        for m in re.finditer(r'(?<!\d)(\d{4})-(\d{2})-(\d{2})(?!\d)', s):
            y, mo, d = int(m.group(1)), int(m.group(2)), int(m.group(3))
            if 1 <= mo <= 12 and 1 <= d <= 31:
                try:
                    datetime(y, mo, d)
                    out.append(f'{y:04d}{mo:02d}{d:02d}')
                except ValueError:
                    pass
    return out


def _pick_datedir_from_filenames(relpaths: list[str]) -> str | None:
    """Pick one ``YYYYMMDD`` from filenames (most common; ties go to latest calendar day)."""
    tokens = _calendar_ymd_tokens_from_strings(relpaths)
    if not tokens:
        return None
    c = Counter(tokens)
    best = max(c.values())
    top = sorted([d for d, n in c.items() if n == best])
    return top[-1]


def _rawdata_dest_name(rel: str, dup_basenames: set[str]) -> str:
    """Basename in rawdata, or ``relpath`` with slashes replaced if basename is duplicated."""
    r = rel.replace('\\', '/')
    base = os.path.basename(r)
    if base in dup_basenames:
        return r.replace('/', '_')
    return base


# Written next to ToOs/<target>/ for goodman-pypeit-reduce --all-to-os (flat layout).
GDRIVE_NIGHT_LABEL_FILE = '.passta_night'


def parse_gdrive_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments for :func:`download_gdrive_toos`."""
    parser = argparse.ArgumentParser(
        description=(
            'Download each immediate subfolder of a shared Google Drive folder into '
            '<outdir>/ToOs/<name>/rawdata (and an empty workspace sibling). A logical '
            'UT night is recorded in ToOs/<name>/.passta_night for batch reduction.'
        ),
    )
    add_version_argument(parser)
    parser.add_argument(
        '--folder-url',
        required=True,
        metavar='URL_OR_ID',
        help='Google Drive folder URL or folder id (required; no default).',
    )
    parser.add_argument(
        '--date',
        type=str,
        default=None,
        metavar='UT_DATE',
        help=(
            'Night label (Astropy Time) for <YYYYMMDD> when no leading YYYYMMDD/ folder '
            'and no parsable calendar dates appear in filenames (otherwise optional; '
            'used before falling back to UTC today).'
        ),
    )
    parser.add_argument('--outdir', type=str, default='.', help='Base output directory (ToOs is created under here).')
    parser.add_argument(
        '--no-cookies',
        action='store_true',
        help='Pass use_cookies=False to gdown (default follows gdown and uses cookies).',
    )
    parser.add_argument(
        '--gdrive-cookies',
        type=str,
        default=None,
        metavar='FILE',
        help=(
            'Netscape-format cookies.txt while logged into Google (gdown reads only '
            '~/.cache/gdown/cookies.txt; this file is copied there). Overrides '
            'PASSTA_GDRIVE_COOKIES / GDRIVE_COOKIES.'
        ),
    )
    add_logging_arguments(parser)
    return parser.parse_args(argv)


def download_gdrive_toos(
    folder_url_or_id: str,
    outdir: str = '.',
    *,
    date: str | None = None,
    use_cookies: bool = True,
    gdrive_cookies: str | None = None,
) -> None:
    """Sync each top-level Drive subdirectory into ``outdir/ToOs/<name>/rawdata``.

    Creates ``workspace`` next to ``rawdata`` (no ``YYYYMMDD`` directory in the output
    path). A single logical night ``YYYY-MM-DD`` is still inferred (same rules as
    before) and written to ``ToOs/<name>/.passta_night`` for
    ``goodman-pypeit-reduce --all-to-os``.

    Files placed directly in the shared folder (not under a subfolder) are skipped
    with a log message. If ``.../rawdata`` for a target already exists and is
    non-empty, that target is skipped.
    """
    import gdown

    cookie_src = _resolve_gdrive_cookies_source(gdrive_cookies)
    if cookie_src:
        if not use_cookies:
            logger.info('Using --gdrive-cookies; enabling gdown cookie session (--no-cookies ignored).')
        _install_gdrive_cookies_for_gdown(cookie_src)
        effective_cookies = True
    else:
        effective_cookies = use_cookies
        if effective_cookies and not os.path.isfile(_GDOWN_COOKIE_CACHE):
            logger.warning(
                'No Google session cookies at %s. Restricted or team Drive shares '
                'usually fail until you export Netscape cookies while logged in and pass '
                '--gdrive-cookies FILE or set PASSTA_GDRIVE_COOKIES (see gdown FAQ).',
                _GDOWN_COOKIE_CACHE,
            )

    folder_id = _parse_gdrive_folder_id(folder_url_or_id)
    url = f'https://drive.google.com/drive/folders/{folder_id}'
    toos = os.path.join(outdir, 'ToOs')
    os.makedirs(toos, exist_ok=True)

    with tempfile.TemporaryDirectory() as td:
        manifest = _fetch_gdrive_folder_manifest(url, td, use_cookies=effective_cookies)
        # Session cookies can make link-shared folders list as empty (wrong Google identity
        # vs anonymous "anyone with the link"). Retry listing without cookies; file downloads
        # still use effective_cookies below.
        if not manifest and effective_cookies:
            logger.warning(
                'Drive folder listing returned no files with session cookies; '
                'retrying listing without cookies (common for link-shared team folders).'
            )
            manifest = _fetch_gdrive_folder_manifest(url, td, use_cookies=False)

    if manifest is None:
        logger.warning(
            'gdown could not retrieve folder contents for %s (HTTP/parse failure). '
            'Confirm the folder id and sharing; see %s',
            folder_id,
            _GDOWN_FAQ,
        )
        return
    if len(manifest) == 0:
        logger.warning(
            'No files found under folder %s (empty listing). If the Drive UI shows data, '
            'check sharing ("Anyone with the link") or try refreshing exported cookies.',
            folder_id,
        )
        return

    by_root: dict[str, list] = {}
    loose: list[str] = []
    for f in manifest:
        p = f.path.replace('\\', '/')
        if '/' not in p:
            loose.append(p)
            continue
        root = p.split('/', 1)[0]
        by_root.setdefault(root, []).append(f)

    if loose:
        logger.info(
            'Skipping %s item(s) at the Drive folder root (not under a subdirectory): %s',
            len(loose),
            ', '.join(loose[:12]) + (' ...' if len(loose) > 12 else ''),
        )

    for root in sorted(by_root):
        relpaths = [_rel_under_drive_root(f.path, root) for f in by_root[root]]
        try:
            inferred = _infer_datedir_from_relpaths(relpaths)
        except ValueError as exc:
            raise SystemExit(f'goodman-spec-download gdrive: {exc}') from exc

        if inferred is not None:
            datedir = inferred
        else:
            from_names = _pick_datedir_from_filenames(relpaths)
            if from_names is not None:
                datedir = from_names
                uniq = set(_calendar_ymd_tokens_from_strings(relpaths))
                if len(uniq) > 1:
                    logger.warning(
                        'Drive folder %r: multiple calendar dates in filenames %s; using %s.',
                        root,
                        ', '.join(sorted(uniq)),
                        datedir,
                    )
                else:
                    logger.info(
                        'Drive folder %r: using YYYYMMDD=%s from filenames.',
                        root,
                        datedir,
                    )
            elif date is not None:
                datedir = Time(date).datetime.strftime('%Y%m%d')
            else:
                datedir = datetime.now(timezone.utc).strftime('%Y%m%d')
                logger.warning(
                    'Drive folder %r: could not infer YYYYMMDD from paths or filenames; '
                    'using UTC today %s. Pass --date UT_DATE to set the night label.',
                    root,
                    datedir,
                )

        target_root = os.path.join(toos, root)
        fulloutdir = os.path.join(target_root, 'rawdata')
        fullworkdir = os.path.join(target_root, 'workspace')
        night_label_iso = f'{datedir[:4]}-{datedir[4:6]}-{datedir[6:8]}'

        n_manifest = len(by_root[root])
        logger.info(
            "Drive folder %r: %s file(s) in manifest -> %s (logical night %s)",
            root,
            n_manifest,
            target_root,
            night_label_iso,
        )

        if _dir_nonempty(fulloutdir):
            n_existing = len(os.listdir(fulloutdir))
            logger.info(
                'Skip %s (%s item(s); delete or empty rawdata to re-fetch)',
                fulloutdir,
                n_existing,
            )
            continue

        os.makedirs(fulloutdir, exist_ok=True)
        os.makedirs(fullworkdir, exist_ok=True)
        sidecar = os.path.join(target_root, GDRIVE_NIGHT_LABEL_FILE)
        with open(sidecar, 'w', encoding='ascii') as sf:
            sf.write(night_label_iso + '\n')

        basenames = [os.path.basename(r.replace('\\', '/')) for r in relpaths]
        dup_basenames = {b for b, n in Counter(basenames).items() if n > 1}

        for f, rel in zip(by_root[root], relpaths, strict=True):
            name = _rawdata_dest_name(rel, dup_basenames)
            out_path = os.path.join(fulloutdir, name)
            if os.path.isfile(out_path) and os.path.getsize(out_path) > 0:
                continue
            logger.info('Download %s', out_path)
            try:
                _call_gdown(
                    gdown.download,
                    id=f.id,
                    output=out_path,
                    quiet=True,
                    use_cookies=effective_cookies,
                    resume=True,
                )
            except Exception as exc:
                if not _is_gdown_file_url_error(exc):
                    raise
                logger.error('%s', _gdrive_access_help(folder_id, f.id, root, str(exc)))
                if _night_tree_only_empty_or_incomplete(fulloutdir, fullworkdir):
                    shutil.rmtree(fulloutdir, ignore_errors=True)
                    shutil.rmtree(fullworkdir, ignore_errors=True)
                    try:
                        os.remove(os.path.join(target_root, GDRIVE_NIGHT_LABEL_FILE))
                    except OSError:
                        pass
                    logger.info('Removed incomplete ToOs tree under %s', target_root)
                else:
                    logger.warning('Leaving partial downloads under %s', target_root)
                raise SystemExit(2) from None


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments for :func:`download_night`."""
    parser = argparse.ArgumentParser(
        description='Download PASSTA Goodman data from the LCO archive for one night.',
    )
    add_version_argument(parser)
    parser.add_argument('date', type=str, help='UT date for the night (parsed by Astropy Time).')
    parser.add_argument(
        '--propid',
        type=str,
        default='SOAR2025B-004',
        help='Comma-separated proposal IDs (default: SOAR2025B-004).',
    )
    parser.add_argument('--outdir', type=str, default='.', help='Base output directory.')
    parser.add_argument(
        '--download-workers',
        type=int,
        default=max(4, min(8, (os.cpu_count() or 4))),
        metavar='N',
        help='Parallel LCO frame downloads per batch (default: min(8, max(4, cpu_count))).',
    )
    add_logging_arguments(parser)
    return parser.parse_args(argv)


def download_night(
    date: str,
    propid: str,
    outdir: str = '.',
    *,
    download_workers: int = 1,
) -> None:
    """Download science, standards, flats, and arcs for each proposal ID."""
    pids = [x.strip() for x in propid.split(',') if x.strip()]
    if not pids:
        logger.warning('No proposal IDs in --propid; nothing to download.')
        return

    outdir = os.path.abspath(os.path.expanduser(outdir))
    datedir = Time(date).datetime.strftime('%Y%m%d')
    fulloutdir = os.path.join(outdir, datedir, 'rawdata')
    logger.info(
        'Starting LCO download for night %s (proposal ID(s): %s) -> %s/%s/rawdata',
        date,
        ', '.join(pids),
        outdir,
        datedir,
    )

    lco = LCOGT()
    _ensure_lco_archive_auth(lco)

    sdate = Time(date) + TimeDelta(14 * 3600 * u.s)
    edate = sdate + TimeDelta(24 * 3600 * u.s)
    fullworkdir = os.path.join(outdir, datedir, 'workspace')
    os.makedirs(fulloutdir, exist_ok=True)
    os.makedirs(fullworkdir, exist_ok=True)

    telid = '4m0a'
    q = {'sdate': sdate, 'edate': edate, 'telid': telid}
    fmt = lco.params['date_format']
    logger.info(
        'Querying LCO archive (telid=%s, window %s .. %s UTC)',
        telid,
        sdate.datetime.strftime(fmt),
        edate.datetime.strftime(fmt),
    )

    callist_raw = lco.get_obslist(
        **q, propid=['calibrate'], rlevel=LCO_RLEVEL_SPECTROSCOPY, obstype='SPECTRUM'
    )
    callist = _filter_goodman_spectrograph_frames(
        _filter_bad_basenames(callist_raw),
        label='calibrate standards',
    )
    logger.info(
        'Calibrate standards: %s archive row(s), %s Goodman after filters',
        len(callist_raw),
        len(callist),
    )

    flatlist_raw = lco.get_obslist(
        **q, propid=['calibrate'], rlevel=LCO_RLEVEL_SPECTROSCOPY, obstype='LAMPFLAT'
    )
    flatlist = _filter_goodman_spectrograph_frames(flatlist_raw, label='calibrate lamp flats')
    logger.info(
        'Calibrate lamp flats: %s archive row(s), %s Goodman',
        len(flatlist_raw),
        len(flatlist),
    )

    arcs_cal_raw = lco.get_obslist(
        **q, propid=['calibrate'], rlevel=[LCO_RLEVEL_RAW], obstype='ARC'
    )
    arcs_cal = _filter_goodman_spectrograph_frames(
        _filter_bad_basenames(arcs_cal_raw),
        label='calibrate arcs',
    )
    logger.info(
        'Calibrate arcs: %s archive row(s), %s Goodman after filters',
        len(arcs_cal_raw),
        len(arcs_cal),
    )

    obs_by_pid: dict[str, list[dict]] = {}
    for pid in pids:
        obslist = lco.get_obslist(**q, propid=[pid], rlevel=LCO_RLEVEL_SPECTROSCOPY, obstype='SPECTRUM')
        good = _goodman_science_frames(obslist)
        logger.info(
            '[%s] science SPECTRUM: %s archive row(s), %s Goodman raw (non-pipeline)',
            pid,
            len(obslist),
            len(good),
        )
        if good:
            obs_by_pid[pid] = list(good)

    if not obs_by_pid:
        logger.warning(
            'No Goodman spectrograph science spectra found for night %s with proposal ID(s) %s '
            'in the archive window above; nothing downloaded.',
            date,
            ', '.join(pids),
        )
        if _rmtree_date_if_unused(fulloutdir, fullworkdir, outdir, datedir):
            logger.info('Removed empty night directory %s', os.path.join(outdir, datedir))
        return

    dw = max(1, download_workers)
    dl_kw = dict(
        outrootdir=fulloutdir,
        skip_header=True,
        funpack=False,
        max_download_workers=dw,
        embed_wavmode=True,
    )

    logger.info(
        'Downloading shared calibrators (%s standard(s), %s flat(s)) with %s worker(s)',
        len(callist),
        len(flatlist),
        dw,
    )
    for lst in (callist, flatlist):
        lco.download_obslist(lst, **dl_kw)

    total_sci = 0
    total_arcs = 0
    for pid, obslist in obs_by_pid.items():
        arcs_prop = _filter_goodman_spectrograph_frames(
            [
                f for f in lco.get_obslist(**q, propid=[pid], rlevel=[LCO_RLEVEL_RAW], obstype='ARC')
                if not _is_pipeline_reduced_frame(f)
            ],
            label=f'{pid} arcs',
        )
        arclist = _dedupe_frames_by_id(arcs_prop + arcs_cal)
        total_sci += len(obslist)
        total_arcs += len(arclist)
        logger.info(
            '[%s] Downloading %s science frame(s) and %s arc(s) with %s worker(s)',
            pid,
            len(obslist),
            len(arclist),
            dw,
        )
        for lst in (obslist, arclist):
            lco.download_obslist(lst, **dl_kw)

    n_files = len(glob.glob(os.path.join(fulloutdir, '*')))
    logger.info(
        'Download complete for %s: %s science frame(s) and %s arc(s) across %s proposal(s); '
        '%s file(s) under %s',
        datedir,
        total_sci,
        total_arcs,
        len(obs_by_pid),
        n_files,
        fulloutdir,
    )

    if _rmtree_date_if_unused(fulloutdir, fullworkdir, outdir, datedir):
        logger.warning(
            'Night directory %s was removed because rawdata and workspace are still empty '
            'after download (check archive authorization or frame URLs).',
            os.path.join(outdir, datedir),
        )


def main(argv: Sequence[str] | None = None) -> None:
    """Console entry point."""
    argv_list = list(sys.argv[1:] if argv is None else argv)
    if argv_list and argv_list[0] == 'gdrive':
        args = parse_gdrive_arguments(argv_list[1:])
        configure_from_parsed_args(args)
        download_gdrive_toos(
            args.folder_url,
            outdir=args.outdir,
            date=args.date,
            use_cookies=not args.no_cookies,
            gdrive_cookies=args.gdrive_cookies,
        )
        return

    args = parse_arguments(argv_list)
    configure_from_parsed_args(args)
    download_night(
        args.date,
        args.propid,
        outdir=args.outdir,
        download_workers=args.download_workers,
    )


if __name__ == '__main__':
    main()
