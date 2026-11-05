"""Run PypeIt on archived PASSTA Goodman data: setup, reduce, flux, and coadd."""

from __future__ import annotations

import os

os.environ.setdefault('MPLBACKEND', 'Agg')

import argparse
import astropy
import contextlib
import glob
import logging
import os
import re
import shutil
import sys
import time
import urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed
from importlib import resources
from typing import Any, Sequence

import numpy as np
from astropy.io import ascii, fits
from astropy.table import Table, unique
from astropy.time import Time

from passta import pypeit_files as pypeit_file_io  # applies PypeIt InputFile read patch
from passta import pypeit_invoke
from pypeit import inputfiles

from passta import slack_utils
from passta.logging_config import add_logging_arguments, configure_from_parsed_args
from passta.lcogt import LCOGT
from passta.spec_download import GDRIVE_NIGHT_LABEL_FILE, _pick_datedir_from_filenames
from passta.versioning import add_version_argument

logger = logging.getLogger(__name__)

_YMD_DIR_RE = re.compile(r'^\d{8}$')


def _ymd_dir_to_date_arg(ymd: str) -> str:
    """``YYYYMMDD`` directory name -> ISO date string for :class:`~astropy.time.Time`."""
    return f'{ymd[:4]}-{ymd[4:6]}-{ymd[6:8]}'


def discover_to_os_reduction_jobs(base_dir: str) -> list[tuple[str, str]]:
    """Discover ``(date, outdir)`` pairs for :func:`run_reduction` under ``ToOs``.

    Supports:

    - **Flat** (current gdrive layout): ``ToOs/<target>/rawdata`` with optional
      ``.passta_night`` (written by :func:`passta.spec_download.download_gdrive_toos`).
    - **Legacy nested**: ``ToOs/<target>/<YYYYMMDD>/rawdata``.

    ``outdir`` for each job is ``base_dir/ToOs/<target>``.

    Returns
    -------
    list[tuple[str, str]]
        Sorted by target name then night label; each item is ``(date_arg, outdir)``.
    """
    toos = os.path.join(base_dir, 'ToOs')
    if not os.path.isdir(toos):
        raise SystemExit(f'goodman-pypeit-reduce --all-to-os: missing directory: {toos}')
    rows: list[tuple[str, str, str]] = []
    for target in sorted(os.listdir(toos)):
        tpath = os.path.join(toos, target)
        if not os.path.isdir(tpath) or target.startswith('.'):
            continue

        nested_nights: list[str] = []
        for ymd in sorted(os.listdir(tpath)):
            if not _YMD_DIR_RE.match(ymd):
                continue
            night = os.path.join(tpath, ymd)
            if not os.path.isdir(night):
                continue
            raw = os.path.join(night, 'rawdata')
            work = os.path.join(night, 'workspace')
            if os.path.isdir(raw) or os.path.isdir(work):
                nested_nights.append(ymd)

        if nested_nights:
            for ymd in sorted(nested_nights):
                rows.append((target, _ymd_dir_to_date_arg(ymd), tpath))
            continue

        flat_raw = os.path.join(tpath, 'rawdata')
        flat_work = os.path.join(tpath, 'workspace')
        if not os.path.isdir(flat_raw) and not os.path.isdir(flat_work):
            continue

        meta = os.path.join(tpath, GDRIVE_NIGHT_LABEL_FILE)
        date_label: str | None = None
        if os.path.isfile(meta):
            with open(meta, encoding='ascii', errors='replace') as sf:
                date_label = sf.readline().strip()
        if not date_label:
            basenames = [
                os.path.basename(p)
                for p in glob.glob(os.path.join(flat_raw, '*.fits'))
                + glob.glob(os.path.join(flat_raw, '*.fz'))
            ]
            ymd_pick = _pick_datedir_from_filenames(basenames)
            if ymd_pick is not None:
                date_label = _ymd_dir_to_date_arg(ymd_pick)
        if not date_label:
            logger.warning(
                'Skipping %s: flat ToOs tree needs %s or parsable dates in frame names',
                tpath,
                GDRIVE_NIGHT_LABEL_FILE,
            )
            continue
        rows.append((target, date_label, tpath))

    rows.sort(key=lambda r: (r[0], r[1]))
    return [(date_label, tpath) for target, date_label, tpath in rows]


@contextlib.contextmanager
def _working_directory(path: str):
    prev = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev)


def _run_shell(cmd: str) -> None:
    logger.info('Running: %s', cmd)
    os.system(cmd)


def _fits_stem(filename: str) -> str:
    return filename.replace('.fits', '').replace('.fz', '')


def _list_rawdata_frames(raw_dir: str) -> list[str]:
    """Return Goodman FITS paths under ``raw_dir`` (``.fz`` preferred over same-stem ``.fits``)."""
    fz = sorted(glob.glob(os.path.join(raw_dir, '*.fz')))
    fits = sorted(glob.glob(os.path.join(raw_dir, '*.fits')))
    fz_stems = {os.path.splitext(os.path.basename(p))[0] for p in fz}
    fits_only = [
        p for p in fits if os.path.splitext(os.path.basename(p))[0] not in fz_stems
    ]
    out = fz + fits_only
    out.sort(key=os.path.basename)
    return out


def _infer_frame_extension(raw_dir: str, preferred: str = '.fits') -> str:
    """Pick the PypeIt ``-e`` suffix from filenames under ``raw_dir``.

    LCO archive downloads are typically ``*.fits.fz``; PypeIt matches them with
    ``-e .fz``. Uncompressed Drive / local trees use ``-e .fits``.
    """
    files = _list_rawdata_frames(raw_dir)
    if not files:
        return preferred
    lower = [p.lower() for p in files]
    if any(p.endswith('.fits.fz') or (p.endswith('.fz') and not p.endswith('.fits')) for p in lower):
        return '.fz'
    if any(p.endswith('.fits') for p in lower):
        return '.fits'
    return preferred


def _goodman_frame_header_hdu_index(hdulist: fits.HDUList, path: str) -> int:
    """HDU index whose header carries Goodman science keys (``WAVMODE``, ``OBSTYPE``, ...).

    SOAR Goodman **uncompressed** ``.fits`` from the archive or Drive are often
    single-extension (metadata on HDU 0). **fpacked** ``.fz`` products are typically
    multi-extension with the science header on HDU 1. We prefer HDU 1 then 0 for
    ``.fz`` / ``.fits.fz``, and HDU 0 then 1 otherwise, falling back to the first
    extension that contains ``WAVMODE``.
    """
    pl = path.lower()
    prefer = (1, 0) if pl.endswith('.fz') or pl.endswith('.fits.fz') else (0, 1)

    def has_goodman_header(i: int) -> bool:
        if i < 0 or i >= len(hdulist):
            return False
        h = hdulist[i]
        hdr = getattr(h, 'header', None)
        return hdr is not None and 'WAVMODE' in hdr

    for i in prefer:
        if has_goodman_header(i):
            return i
    for i in range(len(hdulist)):
        if has_goodman_header(i):
            return i
    raise ValueError(f'No HDU with WAVMODE in {path!r} ({len(hdulist)} extension(s))')


def _goodman_slit_from_header(hdr: fits.Header) -> str | None:
    """Return slit/decker string matching PypeIt metadata ``decker`` (Goodman ``SLIT``)."""
    for key in ('SLIT', 'DECKER'):
        if key in hdr:
            val = str(hdr[key]).strip()
            if val:
                return val.replace(' ', '_')
    return None


def _symlink_if_missing(src: str, dest: str) -> None:
    if not os.path.exists(dest):
        logger.info('Linking %s -> %s', src, dest)
        os.symlink(src, dest)


def _science_object_slit_groups(
    wavedata: Sequence[dict[str, Any]],
) -> list[tuple[str, str]]:
    """Unique ``(object, slit)`` pairs that have on-sky science frames."""
    keys: set[tuple[str, str]] = set()
    for d in wavedata:
        if d['obstype'] in {'SPECTRUM', 'OBJECT'} and d['object'].strip():
            keys.add((d['object'], d['slit']))
    return sorted(keys)


def _link_reduction_directory(
    wavedata: Sequence[dict[str, Any]],
    objdir: str,
    objname: str,
    slit: str,
) -> None:
    """Symlink science and calibrators for one ``(object, slit)`` reduction tree."""
    os.makedirs(objdir, exist_ok=True)
    mjds: list[float] = []
    for data in wavedata:
        if data['slit'] != slit:
            continue
        filename = data['filename']
        basefile = os.path.basename(filename)
        outfile = os.path.join(objdir, basefile)
        obstype = data['obstype']

        if obstype in {'SPECTRUM', 'OBJECT'}:
            if data['object'] != objname:
                continue
            mjds.append(data['mjd'])
            _symlink_if_missing(filename, outfile)
            continue
        if obstype in {'ARC', 'COMP'}:
            continue
        _symlink_if_missing(filename, outfile)

    if not mjds:
        return

    arcdata = [
        d
        for d in wavedata
        if d['slit'] == slit and d['obstype'] in {'ARC', 'COMP'}
    ]
    if not arcdata:
        return

    avg_mjd = float(np.average(mjds))
    arcfile = min(arcdata, key=lambda x: abs(x['mjd'] - avg_mjd))['filename']
    _symlink_if_missing(arcfile, os.path.join(objdir, os.path.basename(arcfile)))


def _build_reduction_dirs(
    fullworkdir: str,
    all_data: dict[str, list[dict[str, Any]]],
) -> list[str]:
    """Create per-``(wavemode, object, slit)`` workspace trees and return their paths."""
    reduction_dirs: list[str] = []
    for wavemode, wavedata in all_data.items():
        wavedir = os.path.join(fullworkdir, wavemode)
        os.makedirs(wavedir, exist_ok=True)
        wavedata_sorted = sorted(wavedata, key=lambda x: x['mjd'])
        for objname, slit in _science_object_slit_groups(wavedata_sorted):
            objdir = os.path.join(wavedir, objname, slit)
            _link_reduction_directory(wavedata_sorted, objdir, objname, slit)
            reduction_dirs.append(objdir)
    return reduction_dirs


def _pypeit_files_in_objdir(spec: str, objdir: str) -> list[str]:
    return sorted(glob.glob(os.path.join(objdir, f'{spec}_*', f'{spec}_*.pypeit')))


def _pypeit_science_products_missing(science_dir: str, pypeit_table: Any) -> bool:
    """True when any science frame lacks a ``spec1d_*.txt`` product.

    This catches the common Goodman failure mode where PypeIt writes a ``spec2d``
    (and may even set ``FINDOBJ``) but optimal extraction yields no 1D spectrum.
    Incomplete Science trees must force a re-reduce **and** a coadd rebuild so a
    later-recovered exposure is not left out of a stale 2-of-N coadd.
    """
    if not os.path.isdir(science_dir):
        return True
    for row in pypeit_table.data:
        if row['frametype'] != 'science':
            continue
        framebase = _fits_stem(row['filename'])
        if not glob.glob(os.path.join(science_dir, f'spec1d_{framebase}*.txt')):
            return True
    return False


def _count_pypeit_science_frames(pypeit_table: Any) -> int:
    return sum(1 for row in pypeit_table.data if row['frametype'] == 'science')


def _coadd_list_matches_objdata(coadd_filename: str, objdata: list[tuple[Any, ...]]) -> bool:
    """Return True when an existing ``.coadd1d`` already lists the same spec1d/obj_id rows."""
    if not os.path.exists(coadd_filename):
        return False
    expected_norm = [
        ' '.join(f'{_spec1d_as_fits(obj[0])} | {obj[1]}'.split()) for obj in objdata
    ]
    try:
        with open(coadd_filename, encoding='utf-8', errors='replace') as handle:
            lines = [ln.strip() for ln in handle]
    except OSError:
        return False
    data_rows = []
    for ln in lines:
        if not ln or ln.lower().startswith('filename') or ln.lower().startswith('coadd1d'):
            continue
        if '.fits' in ln and '|' in ln:
            data_rows.append(' '.join(ln.split()))
    return data_rows == expected_norm


# Goodman lamp cards used by PypeIt ``SOARGoodmanSpectrograph.lamps`` / ``check_frame_type``.
_GOODMAN_LAMP_ARC_CARDS = ('LAMP_HGA', 'LAMP_NE', 'LAMP_AR', 'LAMP_FE', 'LAMP_CU')
_GOODMAN_LAMP_DOME_CARDS = ('LAMP_QUA', 'LAMP_BUL', 'LAMP_DOM')

# ``OBSTYPE`` values that identify calibration or non-science frames (leave unchanged).
_GOODMAN_CALIB_OBSTYPES = frozenset(
    {
        'ZERO',
        'BIAS',
        'DARK',
        'FLAT',
        'DFLAT',
        'DOMEFLAT',
        'PIXFLAT',
        'SKYFLAT',
        'ARC',
        'COMP',
        'COMPARISON',
        'HAL',
        'HALOGEN',
        'CAL',
        'LAMP',
        'WAVE',
        'CONTINUUM',
        'PINHOLE',
        'TRACE',
        'ALIGN',
        'DOME',
        'INTERNAL',
    }
)

# On-sky but not science extractions PypeIt should own via ``OBSTYPE=SPECTRUM``.
_GOODMAN_OBSTYPES_LEAVE_UNCHANGED = frozenset({'SKY', 'ACQUISITION', 'ACQ', 'FOCUS', 'SETUP'})

# Goodman spectroscopy headers use ``GHTS_RED`` / ``GHTS_BLUE``; imager+spec readout uses *_IMAGER.
# If we skip normalization, ``OBSTYPE=OBJECT`` never becomes ``SPECTRUM`` and PypeIt cannot assign
# ``science`` / ``standard`` (``soar_goodman.check_frame_type`` requires ``idname == 'SPECTRUM'``).
_GOODMAN_INSTRUME_FOR_SPECTRUM_OBSTYPE = frozenset(
    {'GHTS_RED', 'GHTS_BLUE', 'GHTS_RED_IMAGER', 'GHTS_BLUE_IMAGER'}
)


def _goodman_header_lamp_truthy(value: Any) -> bool:
    """Return True when a Goodman ``LAMP_*`` card reads as *on* for frametyping."""
    if value is None:
        return False
    s = str(value).strip().lower()
    if s in ('', 'false', 'f', 'no', '0', 'none', 'off', 'n'):
        return False
    if s in ('true', 't', 'yes', '1', 'on', 'y'):
        return True
    try:
        return float(s) != 0.0
    except (TypeError, ValueError):
        return False


def _goodman_header_has_arc_lamps(hdr: fits.Header) -> bool:
    return any(_goodman_header_lamp_truthy(hdr.get(k)) for k in _GOODMAN_LAMP_ARC_CARDS)


def _goodman_header_has_dome_lamps(hdr: fits.Header) -> bool:
    """Any dome / quartz lamp used for Goodman flat-field style exposures."""
    return any(_goodman_header_lamp_truthy(hdr.get(k)) for k in _GOODMAN_LAMP_DOME_CARDS)


def _goodman_header_is_raw_calibration(hdr: fits.Header) -> bool:
    """True when this frame should not receive ``OBSTYPE=SPECTRUM`` normalization."""
    ot = str(hdr.get('OBSTYPE', '')).strip().upper()
    if ot in _GOODMAN_CALIB_OBSTYPES:
        return True
    if ot in _GOODMAN_OBSTYPES_LEAVE_UNCHANGED:
        return True
    if _goodman_header_has_arc_lamps(hdr):
        return True
    if _goodman_header_has_dome_lamps(hdr):
        return True
    try:
        exptime = float(hdr.get('EXPTIME', 0.0) or 0.0)
    except (TypeError, ValueError):
        exptime = 0.0
    if exptime < 0.01 and ot in {'', 'ZERO', 'BIAS'}:
        return True
    return False


def _goodman_normalize_header_for_pypeit_frametyping(hdr: fits.Header) -> str | None:
    """Set ``OBSTYPE=SPECTRUM`` for Goodman on-sky frames so PypeIt can type science/standard.

    PypeIt's Goodman logic requires ``OBSTYPE`` (mapped to ``idname``) to be the
    literal ``SPECTRUM`` for both science and flux standards; ``OBJECT`` and
    other archive values never match. We rewrite only when the header does not
    already describe a bias, flat, arc, or similar calibration (by ``OBSTYPE``
    or lamp status), matching the user's exclusion list for calibration types.

    Returns
    -------
    str or None
        The previous ``OBSTYPE`` string if it was changed, else ``None``.
    """
    instr = str(hdr.get('INSTRUME', '')).strip().upper()
    if instr not in _GOODMAN_INSTRUME_FOR_SPECTRUM_OBSTYPE:
        return None

    prev = str(hdr.get('OBSTYPE', '')).strip()
    prev_u = prev.upper()
    if prev_u == 'SPECTRUM':
        return None
    if _goodman_header_is_raw_calibration(hdr):
        return None

    hdr['OBSTYPE'] = 'SPECTRUM'
    return prev if prev else '(empty)'


def _spec1d_as_fits(filename: str) -> str:
    return filename[:-4] + '.fits' if filename.endswith('.txt') else filename


def _objdata_dir_and_target(objdata: list[tuple[Any, ...]]) -> tuple[str, str]:
    first = objdata[0]
    return os.path.dirname(first[0]), first[2]


def _parallel_each(jobs: int, items: Sequence, fn) -> None:
    n = max(1, jobs)
    if n <= 1 or len(items) <= 1:
        for x in items:
            fn(x)
        return
    w = min(n, len(items))

    def _run(item):
        pypeit_invoke.ensure_headless_matplotlib()
        return fn(item)

    with ThreadPoolExecutor(max_workers=w) as pool:
        futs = [pool.submit(_run, x) for x in items]
        for fut in as_completed(futs):
            fut.result()


def _parallel_map_unordered(jobs: int, items: Sequence, fn):
    n = max(1, jobs)
    if n <= 1 or len(items) <= 1:
        return [fn(x) for x in items]
    w = min(n, len(items))

    def _run(item):
        pypeit_invoke.ensure_headless_matplotlib()
        return fn(item)

    with ThreadPoolExecutor(max_workers=w) as pool:
        futs = [pool.submit(_run, x) for x in items]
        return [fut.result() for fut in as_completed(futs)]


def default_caldb_dir() -> str:
    """Return the filesystem path to the packaged sensitivity-function database."""
    return str(resources.files('passta').joinpath('caldb'))


def refresh_astropy_cache() -> None:
    """Force-download Astropy observatory/site data and pre-load IERS tables."""
    from passta.astropy_iers import ensure_astropy_iers_ready

    tries = 0
    while tries < 5:
        try:
            astropy.coordinates.EarthLocation.get_site_names(refresh_cache=True)
            time.sleep(0.5)
            astropy.coordinates.EarthLocation.get_site_names()
            logger.debug('Astropy EarthLocation site cache refreshed.')
            break
        except urllib.error.URLError as exc:
            tries += 1
            logger.warning(
                'Astropy site cache refresh failed (attempt %s/5): %s',
                tries,
                exc,
            )
            time.sleep(5)
    ensure_astropy_iers_ready()


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments for :func:`run_reduction`.

    Parameters
    ----------
    argv : sequence of str or None
        If None, ``argparse`` reads ``sys.argv``.

    Returns
    -------
    argparse.Namespace
        Parsed ``date``, ``outdir``, ``snid``, ``slack``, ``slack_dry_run``, ``jobs``, etc.
    """
    parser = argparse.ArgumentParser(
        description='Reduce PASSTA Goodman spectra with PypeIt and optional SNID/Slack.',
    )
    add_version_argument(parser)
    parser.add_argument(
        '--all-to-os',
        dest='all_to_os',
        metavar='BASEDIR',
        default=None,
        help=(
            'Run reduction for every ``BASEDIR/ToOs/<name>/`` target (flat ``rawdata`` '
            'or legacy ``<YYYYMMDD>/rawdata``) that has data. Omit DATE and OUTDIR when '
            'using this option.'
        ),
    )
    parser.add_argument(
        'date',
        nargs='?',
        default=None,
        help='Night label matching the dated reduction tree under ``outdir``.',
    )
    parser.add_argument(
        'outdir',
        nargs='?',
        default=None,
        help='Base directory containing ``<YYYYMMDD>/rawdata`` and ``workspace``.',
    )
    parser.add_argument(
        '--snid',
        default=False,
        action='store_true',
        help='Run SNID SAGE classification on telluric-corrected spectra.',
    )
    parser.add_argument(
        '--slack',
        default=False,
        action='store_true',
        help='Post SNID results to Slack (requires ``--snid`` and ``SLACK_TOKEN_FILE_PASSTA``).',
    )
    parser.add_argument(
        '--slack-dry-run',
        dest='slack_dry_run',
        default=False,
        action='store_true',
        help=(
            'Build the same Slack message as ``--slack`` (write ``*.slack.debug`` under '
            '``spectra/slack/``) and log it; run SNID; do not call the Slack API or read tokens.'
        ),
    )
    parser.add_argument(
        '--caldb',
        type=str,
        default=None,
        help='Directory with caldb.txt and sensitivity functions (default: packaged caldb).',
    )
    parser.add_argument(
        '--jobs',
        type=int,
        default=1,
        metavar='N',
        help='Parallel PypeIt setup/reduce/SNID workers (default 1). Use 2--4 on multi-core hosts.',
    )
    parser.add_argument(
        '--frame-extension',
        default=None,
        help=(
            'Suffix passed to ``pypeit_setup`` for frames linked into each object '
            'directory (default: infer ``.fz`` for ``*.fits.fz`` LCO downloads, '
            'else ``.fits`` for uncompressed data).'
        ),
    )
    parser.add_argument(
        '--slit-pos',
        type=float,
        default=505.0,
        metavar='PIX',
        dest='slit_pos',
        help=(
            'Spatial pixel position (spat_pixpos) for choosing the science trace in '
            'parse_pypeit_output (default 505). Pixel tolerance is still widened in '
            'steps of 15, 30, then 60.'
        ),
    )
    parser.add_argument(
        '--no-try-old',
        dest='try_old',
        action='store_false',
        help=(
            'Enforce PypeIt ``SpecObj`` datamodel version checks during ``pypeit_flux_calib`` '
            '(default: pass ``--try_old``, skipping version checks).'
        ),
    )
    parser.add_argument(
        '--overwrite',
        default=False,
        action='store_true',
        help=(
            'Force PypeIt re-setup/re-reduce and regenerate flux/coadd products even when '
            'Science outputs already exist. Incomplete reductions (spec2d without spec1d) '
            'are always re-run with overwrite.'
        ),
    )
    parser.set_defaults(try_old=True)
    add_logging_arguments(parser)
    args = parser.parse_args(argv)
    if args.all_to_os is not None:
        args.all_to_os = os.path.abspath(os.path.expanduser(args.all_to_os))
        if args.date is not None or args.outdir is not None:
            parser.error('With --all-to-os BASEDIR, omit DATE and OUTDIR positionals.')
    else:
        if args.date is None or args.outdir is None:
            parser.error('Expected DATE OUTDIR (or use --all-to-os BASEDIR alone).')
        args.outdir = os.path.abspath(os.path.expanduser(args.outdir))
    if args.slack_dry_run:
        args.snid = True
    if args.frame_extension is not None and not str(args.frame_extension).startswith('.'):
        parser.error('--frame-extension must start with a dot (e.g. .fits or .fz)')
    return args


def _pypeit_setup_objdir(
    spec: str,
    objdir: str,
    frame_extension: str = '.fits',
    *,
    overwrite: bool = False,
) -> None:
    """Run PypeIt setup for one reduction directory via :mod:`passta.pypeit_invoke`.

    When ``overwrite`` is True, remove existing ``{spec}_*`` setup directories so the
    packaged Goodman parameter block (flat-field / findobj / extraction) is re-applied.
    """
    if overwrite:
        for setup_dir in sorted(glob.glob(os.path.join(objdir, f'{spec}_*'))):
            if os.path.isdir(setup_dir):
                logger.info('Removing setup directory for overwrite: %s', setup_dir)
                shutil.rmtree(setup_dir)
    pypeit_invoke.run_setup(spec, objdir, objdir, extension=frame_extension)


def _pypeit_reduce_objdir(
    spec: str,
    objdir: str,
    caldb_dir: str,
    *,
    slit_pos: float = 505.0,
    overwrite: bool = False,
) -> tuple[list[tuple[Any, ...]], bool]:
    """Run PypeIt reduction if needed and return science rows from :func:`parse_pypeit_output`.

    Returns
    -------
    objdata, force_rebuild_coadd
        ``force_rebuild_coadd`` is True when Science products were incomplete and a
        re-reduce recovered frames that a stale coadd might still omit.
    """
    pypeit_files = _pypeit_files_in_objdir(spec, objdir)
    if not pypeit_files:
        logger.error('pypeit_setup failed or missing .pypeit for %s', objdir)
        return [], False

    science_dir = os.path.join(objdir, 'Science')
    objdata: list[tuple[Any, ...]] = []
    force_rebuild_coadd = False
    for pypeit_file in pypeit_files:
        pypeit_file_io.sanitize_pypeit_file(pypeit_file)
        pypeit_table = inputfiles.PypeItFile.from_file(pypeit_file)
        setup = os.path.basename(os.path.dirname(pypeit_file))
        n_science = _count_pypeit_science_frames(pypeit_table)
        missing = _pypeit_science_products_missing(science_dir, pypeit_table)
        if overwrite or missing:
            if missing and not overwrite:
                logger.warning(
                    'Incomplete Science products for %s (%s); re-running to recover '
                    'missing extractions (will force coadd rebuild)',
                    objdir,
                    setup,
                )
                force_rebuild_coadd = True
            # Full overwrite must rebuild flats too: packaged Goodman flatfield /
            # findobj / extraction params live in the .pypeit file and only take
            # effect when Calibrations are regenerated.
            if overwrite:
                force_rebuild_coadd = True
                calib_dir = os.path.join(objdir, 'Calibrations')
                if os.path.isdir(calib_dir):
                    logger.info('Removing Calibrations for overwrite: %s', calib_dir)
                    shutil.rmtree(calib_dir)
                qa_dir = os.path.join(objdir, 'QA')
                if os.path.isdir(qa_dir):
                    shutil.rmtree(qa_dir)
                if os.path.isdir(science_dir):
                    logger.info('Removing Science for overwrite: %s', science_dir)
                    shutil.rmtree(science_dir)
            pypeit_invoke.run_pypeit(
                pypeit_file,
                objdir,
                overwrite=True,
                reuse_calibs=not overwrite,
            )
        else:
            logger.info('Reduction already done for %s (%s)', objdir, setup)
        recovered = parse_pypeit_output(science_dir, pypeit_table, caldb_dir, slit_pos=slit_pos)
        if n_science and len(recovered) < n_science:
            logger.error(
                'Only %s/%s science frames have usable 1D traces for %s (%s); '
                'coadd will omit incomplete extractions',
                len(recovered),
                n_science,
                objdir,
                setup,
            )
        objdata.extend(recovered)
    return objdata, force_rebuild_coadd


def parse_pypeit_output(
    scidir: str,
    pypeit_table: Any,
    caldb: str,
    slit_pos: float = 505.0,
    pix_tolerance: Sequence[float] = (15.0, 30.0, 60.0),
) -> list[tuple[Any, ...]]:
    """Inspect PypeIt ``Science`` products and build sensitivity functions.

    For each science row, selects the highest-S/N trace near ``slit_pos``. For
    standards, runs sensitivity-function generation via PypeIt and appends metadata to ``caldb.txt``.

    Parameters
    ----------
    scidir : str
        Path to the ``Science`` subdirectory for this reduction.
    pypeit_table : `pypeit.inputfiles.PypeItFile`
        Parsed PypeIt file object (``.pypeit``).
    caldb : str
        Directory holding ``caldb.txt`` and sensitivity FITS files.
    slit_pos : float, optional
        Expected spatial pixel position of the target trace (default 505).
    pix_tolerance : sequence of float, optional
        Increasing spatial tolerances tried when matching traces.

    Returns
    -------
    list of tuple
        Per-object tuples
        ``(spec1d_txt, objname, target, mjd, dispname, slitname)`` for science
        extractions; may be empty if no valid traces.
    """
    objdata = []
    for row in pypeit_table.data:
        if row['frametype'] == 'science':
            filename = row['filename']
            target = row['target'].lower()
            ra = row['ra']
            dec = row['dec']
            mjd = row['mjd']
            airmass = row['airmass']
            dispname = row['mode']
            slitname = row['decker']

            framebase = _fits_stem(filename)
            outfiles = glob.glob(os.path.join(scidir, f'spec1d_{framebase}*.txt'))

            if len(outfiles) == 0:
                logger.warning('No traces detected for %s', framebase)
                continue

            outfile = outfiles[0]
            spectable = ascii.read(outfile, delimiter='|', header_start=0)

            subtable = None
            for pix_tol in pix_tolerance:
                mask = np.abs(spectable['spat_pixpos'] - slit_pos) < pix_tol
                subtable = spectable[mask]
                logger.debug(
                    'Slit match table for slit_pos=%s pix_tol=%s:\n%s',
                    slit_pos,
                    pix_tol,
                    subtable,
                )
                if len(subtable) > 0:
                    break

            if subtable is None or len(subtable) == 0:
                logger.warning('Could not find a good trace in: %s', filename)
                continue

            idx = int(np.argmax(subtable['s2n'].data))
            objname = subtable[idx]['name']
            objdata.append((outfile, objname, target, mjd, dispname, slitname))

        elif row['frametype'] == 'standard':
            framebase = _fits_stem(row['filename'])
            spec1dfiles = glob.glob(os.path.join(scidir, f'spec1d_{framebase}*.fits'))
            if len(spec1dfiles) == 0:
                continue

            outfile = spec1dfiles[0]
            target = row['target'].lower()
            ra = row['ra']
            dec = row['dec']
            mjd = row['mjd']
            airmass = row['airmass']
            datestr = Time(mjd, format='mjd').datetime.strftime('%Y%m%d')
            dispname = row['mode']
            slitname = row['decker']

            sensfilename = f'{target}.{dispname}.{slitname}.{datestr}_sensfunc.fits'
            fullsensfilename = os.path.join(caldb, sensfilename)
            par_outfile = os.path.join(scidir, 'sensfunc.par')

            os.makedirs(caldb, exist_ok=True)
            os.makedirs(os.path.join(caldb, 'plots'), exist_ok=True)

            pypeit_invoke.run_sensfunc(outfile, fullsensfilename, par_outfile)

            for file in glob.glob(os.path.join(caldb, '*.pdf')):
                newfile = os.path.join(caldb, 'plots', os.path.basename(file))
                if not os.path.exists(newfile):
                    shutil.move(file, os.path.join(caldb, 'plots'))
                else:
                    os.remove(file)

            if os.path.exists(fullsensfilename):
                caldbfile = os.path.join(caldb, 'caldb.txt')
                names = (
                    'filename', 'target', 'ra', 'dec', 'mjd',
                    'dispname', 'slitname', 'airmass',
                )
                dtype = (str, str, np.float64, np.float64, np.float64, str, str, np.float64)
                row_vals = (
                    os.path.basename(sensfilename),
                    target,
                    ra,
                    dec,
                    mjd,
                    dispname,
                    slitname,
                    airmass,
                )
                if not os.path.exists(caldbfile):
                    caldb_table = Table(names=names, dtype=dtype)
                else:
                    caldb_table = ascii.read(caldbfile)
                caldb_table.add_row(row_vals)
                caldb_table = unique(caldb_table)
                caldb_table.write(caldbfile, format='ascii.ecsv', overwrite=True)

    return objdata


def make_flux_file(objdata: list[tuple[Any, ...]], fullsensfile: str) -> str | None:
    """Write a ``pypeit_flux_calib`` input file for the listed spec1d products.

    Parameters
    ----------
    objdata : list of tuple
        Rows as produced by :func:`parse_pypeit_output`.
    fullsensfile : str
        Absolute path to the sensitivity function FITS file.

    Returns
    -------
    str or None
        Path to the written ``.flux`` file, or None if nothing was written.
    """
    dirname, target = _objdata_dir_and_target(objdata)
    flux_filename = os.path.join(dirname, f'{target}.flux')

    with open(flux_filename, 'w') as f:
        # ``extrap_sens`` allows science wavelengths outside the std sensfunc grid (e.g. Goodman
        # coverage bluer than the standard); otherwise ``pypeit_flux_calib`` raises PypeItError.
        # ``extinct_correct`` must match the UVIS sensfunc default (True): the packaged
        # sensitivity functions are built with extinction-corrected standards, so science
        # frames need the same correction or the blue continuum is systematically tilted.
        f.write('[fluxcalib] \n')
        f.write('  extinct_correct = True \n')
        f.write('  extinct_file = closest \n')
        f.write('  extrap_sens = True \n\n')
        f.write('flux read \n')
        f.write(f'  path {dirname} \n')
        f.write('filename | sensfile \n')
        for obj in objdata:
            f.write(f'{_spec1d_as_fits(obj[0])} | {fullsensfile} \n')
        f.write('flux end \n')

    return flux_filename if os.path.exists(flux_filename) else None


def make_coadd_file(objdata: list[tuple[Any, ...]], fullout_1dfilename: str) -> str | None:
    """Write a ``pypeit_coadd_1dspec`` file list for one target.

    Parameters
    ----------
    objdata : list of tuple
        Rows as produced by :func:`parse_pypeit_output`.
    fullout_1dfilename : str
        Desired output coadd path passed into the config block.

    Returns
    -------
    str or None
        Path to the ``.coadd1d`` file if written.
    """
    dirname, target = _objdata_dir_and_target(objdata)
    coadd_filename = os.path.join(dirname, f'{target}.coadd1d')

    with open(coadd_filename, 'w') as f:
        f.write('[coadd1d] \n')
        f.write(f'  coaddfile = {fullout_1dfilename} \n')
        f.write('  wave_method = linear \n\n')
        f.write('coadd1d read \n')
        f.write(f'  path {dirname} \n')
        f.write('filename | obj_id \n')
        for obj in objdata:
            f.write(f'{_spec1d_as_fits(obj[0])} | {obj[1]} \n')
        f.write('coadd1d end \n')

    return coadd_filename if os.path.exists(coadd_filename) else None


def pypeit_post_process(
    objdata: list[tuple[Any, ...]],
    specdir: str,
    caldb_dir: str,
    *,
    try_old: bool = True,
    overwrite: bool = False,
    force_rebuild_coadd: bool = False,
) -> str | None:
    """Flux-calibrate, coadd, clean the blue end, and symlink the final 1D spectrum.

    Telluric correction is disabled in the pipeline (see ``passta.tellfit``);
    run ``goodman-pypeit-tellfit`` on the coadd when needed.

    Parameters
    ----------
    objdata : list of tuple
        Science rows for a single object from :func:`parse_pypeit_output`.
    specdir : str
        Directory where the coadd FITS symlink is created.
    caldb_dir : str
        Location of ``caldb.txt`` and sensitivity functions.
    try_old : bool, optional
        Passed to ``pypeit_flux_calib``.
    overwrite : bool, optional
        If True, regenerate flux/coadd products even when the coadd FITS exists.
    force_rebuild_coadd : bool, optional
        If True, rebuild the coadd even without ``overwrite`` (e.g. after an
        incomplete Science tree was repaired and newly recovered frames must
        enter the coadd).

    Returns
    -------
    str or None
        Path to the flux-calibrated coadd spectrum if produced, else None.

    Raises
    ------
    FileNotFoundError
        If ``caldb/caldb.txt`` is missing.
    """
    caldb_file = os.path.join(caldb_dir, 'caldb.txt')
    if not os.path.exists(caldb_file):
        raise FileNotFoundError(f'Could not locate calibration database: {caldb_file}')

    caldb = ascii.read(caldb_file)

    firstobj = objdata[0]
    outfile, _objname, target, mjd, dispname, slitname = firstobj
    dirname, _ = os.path.split(outfile)
    datestr = Time(mjd, format='mjd').datetime.strftime('%Y%m%d')

    mask = (caldb['dispname'] == dispname) & (caldb['slitname'] == slitname)
    goodcals = caldb[mask]
    idx = int(np.argmin(np.abs(goodcals['mjd'] - mjd)))
    fullsensfile = os.path.join(caldb_dir, goodcals[idx]['filename'])

    flux_filename = make_flux_file(objdata, fullsensfile)
    flux_par_outfile = os.path.join(dirname, f'{target}.fluxing.par')
    if flux_filename is None:
        return None

    pypeit_invoke.run_flux_calib(flux_filename, flux_par_outfile, try_old=try_old)

    out_1dfilename = f'{target}.{datestr}.{dispname}.coadd1d.fits'
    fullout_1dfilename = os.path.join(dirname, out_1dfilename)
    full_tellcorrfilename = fullout_1dfilename.replace('.fits', '_tellcorr.fits')

    # Detect stale coadds written before a later-recovered exposure existed.
    existing_coadd_cfg = os.path.join(dirname, f'{target}.coadd1d')
    stale_coadd = os.path.exists(fullout_1dfilename) and not _coadd_list_matches_objdata(
        existing_coadd_cfg, objdata
    )
    rebuild_coadd = bool(overwrite or force_rebuild_coadd or stale_coadd)
    if stale_coadd and not (overwrite or force_rebuild_coadd):
        logger.warning(
            'Existing coadd for %s lists a different exposure set than current '
            'extractions (%s frames); rebuilding coadd',
            target,
            len(objdata),
        )

    coadd_filename = make_coadd_file(objdata, fullout_1dfilename)
    par_outfile = os.path.join(dirname, 'coadd1d.par')
    if coadd_filename is None:
        return None

    if rebuild_coadd and os.path.exists(fullout_1dfilename):
        logger.info('Removing existing coadd for rebuild: %s', fullout_1dfilename)
        os.remove(fullout_1dfilename)
    if rebuild_coadd and os.path.exists(full_tellcorrfilename):
        os.remove(full_tellcorrfilename)

    if not os.path.exists(fullout_1dfilename):
        pypeit_invoke.run_coadd_1dspec(coadd_filename, par_outfile)

    if os.path.exists(fullout_1dfilename):
        from passta.spectrum_cleanup import clean_goodman_coadd_blue_end

        try:
            clean_goodman_coadd_blue_end(fullout_1dfilename, inplace=True)
        except Exception:
            logger.exception('Blue-end cleanup failed for %s', fullout_1dfilename)

    from passta.tellfit import TELLURIC_CORRECTION_ENABLED

    final_path = full_tellcorrfilename if TELLURIC_CORRECTION_ENABLED else fullout_1dfilename
    if TELLURIC_CORRECTION_ENABLED and os.path.exists(fullout_1dfilename):
        tell_par_outfile = os.path.join(dirname, 'tellfit.par')
        tell_file = pypeit_invoke.goodman_telluric_param_file()
        if tell_file:
            logger.info('Telluric parameter file %s', tell_file)
        if not os.path.exists(full_tellcorrfilename):
            with _working_directory(dirname):
                pypeit_invoke.run_tellfit(
                    fullout_1dfilename,
                    objmodel='poly',
                    par_outfile=tell_par_outfile,
                    tell_file=tell_file,
                )
    elif os.path.exists(fullout_1dfilename):
        logger.info(
            'Telluric correction skipped (pipeline disabled); use goodman-pypeit-tellfit on %s',
            fullout_1dfilename,
        )

    if os.path.exists(final_path):
        outlink = os.path.join(specdir, os.path.basename(final_path))
        if os.path.islink(outlink) or os.path.exists(outlink):
            if rebuild_coadd:
                os.remove(outlink)
            elif os.path.exists(outlink):
                logger.info('Linked %s -> %s', final_path, outlink)
                return outlink
        if not os.path.exists(outlink):
            os.symlink(final_path, outlink)
        if os.path.exists(outlink):
            logger.info('Linked %s -> %s', final_path, outlink)
            return outlink
    return None


def post_slack_results(
    file_path: str,
    results: dict[str, Any],
    slackdir: str,
    channel: str = 'passta-classifications',
    *,
    dry_run: bool = False,
) -> None:
    """Format SNID results and optionally upload the spectrum (and plot) to Slack.

    Parameters
    ----------
    file_path : str
        Path to the telluric-corrected FITS spectrum.
    results : dict
        Output structure from :func:`run_snid_sage`.
    slackdir : str
        Directory for the marker / debug text file.
    channel : str, optional
        Slack channel name (without leading ``#``).
    dry_run : bool, optional
        If True, write ``*.slack.debug`` and log the message; do not use the Slack API.

    Raises
    ------
    RuntimeError
        If ``SLACK_TOKEN_FILE_PASSTA`` is not set (live post only).
    """
    basename = os.path.basename(file_path)
    suffix = '.slack.debug' if dry_run else '.slack'
    slack_file = os.path.join(slackdir, basename.replace('.fits', suffix))

    if not dry_run and os.path.exists(slack_file):
        logger.info('Already posted to Slack; skipping %s', slack_file)
        return

    with fits.open(file_path) as hdu:
        datestr = Time(hdu['PRIMARY'].header['MJD'], format='mjd').datetime.strftime(
            '%Y-%m-%dT%H:%M:%S'
        )
        mode = hdu['PRIMARY'].header['MODE']
        objname = hdu['PRIMARY'].header['TARGET']

    res_table = results.get('results')
    summary = results.get('summary')
    if summary:
        best_class = f"{summary.get('type', '?')} {summary.get('subtype', '')}".strip()
        conf = summary.get('type_confidence')
        if conf:
            best_class = f'{best_class} ({conf})'
    elif res_table is not None and len(res_table) > 0:
        best_class = res_table[0]['subtype']
    else:
        best_class = None

    message = f'*PASSTA Spectrum*: {objname} \n'
    message += f'*Observation Time*: {datestr} \n'
    message += f'*Mode*: {mode} \n'
    message += f'*Best classification*: {best_class} \n'

    if res_table is not None and len(res_table) > 0:
        message += '*SNID SAGE Classifications:* \n'
        form = '{name: <9} {typ: <4} {subtyp: <8} {redshift: <10} {age: <8}'
        message += '```'
        message += form.format(
            name='name', typ='type', subtyp='subtype', redshift='z', age='age'
        )
        message += '\n'
        nresults = int(np.min([len(res_table), 3]))
        for i in range(nresults):
            row = res_table[i]
            message += form.format(
                name=row['name'],
                typ=row['type'],
                subtyp=row['subtype'],
                redshift=row['redshift'],
                age=row['age'],
            )
            message += '\n'
        message += '```'

    with open(slack_file, 'w') as f:
        f.write(message)

    if dry_run:
        logger.info('Slack dry-run (no API); message written to %s\n%s', slack_file, message)
        return

    token_file = os.environ.get('SLACK_TOKEN_FILE_PASSTA')
    if not token_file:
        raise RuntimeError('SLACK_TOKEN_FILE_PASSTA needs to be an environment variable')

    client = slack_utils.set_up_slack(token_file)
    file_ids = [slack_utils.get_file_upload_url(client, file_path, channel)]
    if 'plot' in results and results['plot']:
        file_ids.append(slack_utils.get_file_upload_url(client, results['plot'], channel))

    slack_utils.post_slack_files(client, file_ids, message, channel)


def _parse_snid_sage_summary(text: str) -> dict[str, str] | None:
    """Parse the SNID-SAGE classification summary block at the top of ``.output``."""
    summary: dict[str, str] | None = None
    for line in text.splitlines():
        ls = line.strip()
        if ls.startswith('|') and 'Type' not in ls and not set(ls.replace('|', '').strip()) <= {'-'}:
            parts = [p.strip() for p in ls.strip('|').split('|')]
            if len(parts) >= 5:
                typ, subtype, z, z_err, age = parts[:5]
                if typ.lower() != 'type':
                    summary = {
                        'type': typ,
                        'subtype': subtype,
                        'redshift': z,
                        'redshift_err': z_err,
                        'age_days': age,
                    }
        # e.g. "Type       Ia          High            Medium (+49.3%, vs Ic)"
        if summary is not None and ls.lower().startswith('type') and 'vs' in ls.lower():
            parts = ls.split()
            # Label MatchQuality Confidence...
            if len(parts) >= 3:
                summary['type_match_quality'] = parts[2]
            if len(parts) >= 4:
                summary['type_confidence'] = parts[3]
    return summary


def _log_snid_sage_results(
    filepath: str,
    result_file: str,
    logpath: str,
    out: dict[str, Any],
    *,
    ran_identify: bool,
) -> None:
    """Emit a concise SNID summary to the pipeline log (stdout is redirected during identify)."""
    basename = os.path.basename(filepath)
    if ran_identify:
        logger.info('SNID SAGE finished for %s (log: %s)', basename, logpath)
    else:
        logger.info('SNID SAGE reusing cached results for %s', basename)

    summary: dict[str, str] | None = None
    if os.path.exists(result_file):
        with open(result_file, encoding='utf-8', errors='replace') as handle:
            summary = _parse_snid_sage_summary(handle.read())

    if summary:
        logger.info(
            'SNID classification for %s: %s %s  z=%s+/-%s  age=%s d',
            basename,
            summary['type'],
            summary['subtype'],
            summary['redshift'],
            summary['redshift_err'],
            summary['age_days'],
        )
    else:
        logger.warning('SNID produced no parseable summary for %s (%s)', basename, result_file)

    tbl = out.get('results')
    if tbl is not None and len(tbl) > 0:
        top = tbl[0]
        logger.info(
            'SNID top template for %s: %s (%s %s)  z=%s  RLAP=%s  age=%s d',
            basename,
            top['name'],
            top['type'],
            top['subtype'],
            top['redshift'],
            top['rlap-ccc'],
            top['age'],
        )

    plot = out.get('plot')
    if plot:
        logger.info('SNID plot for %s: %s', basename, plot)


def _parse_snid_sage_output_pypi(text: str) -> Table | None:
    """Parse PyPI ``snid-sage`` text ``.output`` (TEMPLATE MATCHES table)."""
    if 'TEMPLATE MATCHES' not in text:
        return None
    rows: list[list[str]] = []
    capture = False
    for line in text.splitlines():
        ls = line.strip()
        if '# Template' in ls and 'Subtype' in ls:
            capture = True
            continue
        if not capture:
            continue
        if not ls:
            if rows:
                break
            continue
        if set(ls) <= {'-', '|'}:
            continue
        segs = [s for s in re.split(r'\s{2,}', line.strip()) if s]
        if len(segs) < 7:
            continue
        try:
            idx_name = segs[0].split(None, 1)
            if len(idx_name) != 2:
                continue
            idx_s, name = idx_name
            int(idx_s)
        except (ValueError, IndexError):
            continue
        typ, subtype, rlap, z, ze, age = segs[1], segs[2], segs[3], segs[4], segs[5], segs[6]
        rows.append([idx_s, name, typ, subtype, rlap, z, ze, age])
    if not rows:
        return None
    columns = [list(row) for row in zip(*rows)]
    return Table(
        columns,
        names=('idx', 'name', 'type', 'subtype', 'rlap-ccc', 'redshift', 'error', 'age'),
    )


def _parse_snid_sage_output_legacy(text: str) -> Table | None:
    """Parse legacy SNID SAGE ``.output`` (table after a ``--------------------`` rule)."""
    data_start = False
    all_data: list[list[str]] = []
    for line in text.splitlines():
        if line.startswith('--------------------'):
            data_start = True
        elif not data_start:
            continue
        else:
            all_data.append(line.split())
    if not all_data or not all(len(r) > 0 for r in all_data):
        return None
    columns = [list(row) for row in zip(*all_data)]
    return Table(
        columns,
        names=('idx', 'name', 'type', 'subtype', 'rlap-ccc', 'redshift', 'error', 'age'),
    )


def _parse_snid_sage_output_file(result_file: str) -> Table | None:
    """Load ``.output`` from disk; prefer PyPI ``snid-sage`` format when detected."""
    with open(result_file, encoding='utf-8', errors='replace') as outf:
        text = outf.read()
    if 'TEMPLATE MATCHES' in text:
        tbl = _parse_snid_sage_output_pypi(text)
        if tbl is not None:
            return tbl
    return _parse_snid_sage_output_legacy(text)


def _run_snid_sage_identify_in_process(filepath: str, results_dir: str, logpath: str) -> int:
    """Call **snid-sage**'s CLI entry point in-process (same as ``sage identify ...``).

    Uses ``snid_sage.interfaces.cli.main.main`` so behavior matches the installed
    ``sage`` console script without a subprocess. Stdout/stderr from the run are
    appended to ``logpath``.
    """
    try:
        pypeit_invoke.ensure_headless_matplotlib()
        import matplotlib

        matplotlib.use('Agg')
    except Exception:
        pass

    from snid_sage.interfaces.cli.main import main as snid_sage_cli_main

    # Goodman 400_M1/M2: residual blue-end ghost (~3800 A) can dominate the
    # cross-correlation and flip Ia<->Ic at very low confidence. Mask that band
    # and the O2 A-band by default; the coadd is also cleaned in-place first.
    argv = [
        '--quiet',
        'identify',
        filepath,
        '-o',
        results_dir,
        '--no-progress',
        '--aband-remove',
        '--wavelength-masks',
        '3000:3800',
        '7550:7700',
    ]
    logger.info('Running SNID SAGE in-process: sage %s', ' '.join(argv))

    rc = 1
    with open(logpath, 'w', encoding='utf-8', errors='replace') as logf:
        logf.write(f'# argv: sage {" ".join(argv)}\n')
        try:
            with contextlib.redirect_stdout(logf), contextlib.redirect_stderr(logf):
                out_rc = snid_sage_cli_main(argv)
            rc = 0 if out_rc is None else int(out_rc)
        except SystemExit as exc:
            code = exc.code
            if code is None:
                rc = 0
            elif isinstance(code, int):
                rc = code
            else:
                rc = 1
        except Exception:
            rc = 1
            logger.exception('SNID SAGE identify failed (log: %s)', logpath)
    if rc != 0:
        logger.error('SNID SAGE identify returned code %s (log: %s)', rc, logpath)
    return rc


def run_snid_sage(filepath: str, *, overwrite: bool = False) -> dict[str, Any]:
    """Run **snid-sage** identification on a 1D spectrum and parse outputs.

    Calls the same code path as the ``sage identify ... -o`` CLI via
    :func:`_run_snid_sage_identify_in_process` (no subprocess). Legacy ``.output``
    formats are still parsed if present.

    Parameters
    ----------
    filepath : str
        Flux-calibrated (optionally telluric-corrected) FITS spectrum.
    overwrite : bool, optional
        If True, delete a cached ``.output`` and re-run identify.

    Returns
    -------
    dict
        Keys ``results`` (`~astropy.table.Table` or None), ``summary`` (dict or
        None), and ``plot`` (path string or None).
    """
    dirpath, basename = os.path.split(filepath)
    results_dir = os.path.join(dirpath, 'results')
    result_file = os.path.join(results_dir, basename.replace('.fits', '.output'))
    spec_plot = os.path.join(
        results_dir, basename.replace('.fits', '_flattened_spectrum.png')
    )

    os.makedirs(results_dir, exist_ok=True)
    logpath = os.path.join(results_dir, basename.replace('.fits', '.log'))

    if overwrite and os.path.exists(result_file):
        os.remove(result_file)

    ran_identify = False
    if not os.path.exists(result_file):
        try:
            ran_identify = _run_snid_sage_identify_in_process(filepath, results_dir, logpath) == 0
        except ImportError as exc:
            logger.error('snid-sage is required for --snid: %s', exc)
            return {'results': None, 'summary': None, 'plot': None}

    out: dict[str, Any] = {'results': None, 'summary': None, 'plot': None}

    if os.path.exists(result_file):
        with open(result_file, encoding='utf-8', errors='replace') as handle:
            text = handle.read()
        out['summary'] = _parse_snid_sage_summary(text)
        if 'TEMPLATE MATCHES' in text:
            out['results'] = _parse_snid_sage_output_pypi(text)
            if out['results'] is None:
                out['results'] = _parse_snid_sage_output_legacy(text)
        else:
            out['results'] = _parse_snid_sage_output_legacy(text)

    if os.path.exists(spec_plot):
        out['plot'] = spec_plot

    _log_snid_sage_results(filepath, result_file, logpath, out, ran_identify=ran_identify)

    return out


def run_reduction(
    date: str,
    outdir: str,
    caldb_dir: str | None = None,
    snid: bool = False,
    slack: bool = False,
    slack_dry_run: bool = False,
    parallel_jobs: int = 1,
    *,
    frame_extension: str | None = None,
    slit_pos: float = 505.0,
    try_old: bool = True,
    overwrite: bool = False,
) -> None:
    """Main driver: symlink raw data, run PypeIt setup and reduction in-process, post-process.

    Parameters
    ----------
    date : str
        Night label; must match the dated tree created by the download step.
    outdir : str
        Base reduction directory (absolute path recommended).
    caldb_dir : str or None
        Sensitivity-function database directory; defaults to the packaged
        :func:`default_caldb_dir`.
    snid : bool, optional
        If True, run :func:`run_snid_sage` on outputs.
    slack : bool, optional
        If True with ``snid``, push SNID summaries via Slack (requires token env).
    slack_dry_run : bool, optional
        If True with ``snid``, build the same Slack message to ``*.slack.debug`` and log
        it; do not call the Slack API (``--slack`` is ignored when this is True).
    parallel_jobs : int, optional
        If greater than 1, run independent PypeIt setup, reduction, and SNID
        steps with a thread pool (default 1). Keep small (2--4) on NFS hosts.
    frame_extension : str or None, optional
        Filename suffix for ``pypeit_setup`` (``-e``). If omitted, infer ``.fz`` when
        ``rawdata`` holds ``*.fits.fz`` LCO frames, otherwise ``.fits``.
    slit_pos : float, optional
        Passed to :func:`parse_pypeit_output` to select the science trace nearest this
        ``spat_pixpos`` (default 505).
    try_old : bool, optional
        If True (default), run ``pypeit_flux_calib`` with ``--try_old``. Pass ``--no-try-old``
        on the CLI to enforce ``SpecObj`` datamodel version checks.
    overwrite : bool, optional
        Force re-setup, re-reduce, re-flux/coadd, and re-SNID even when products exist.

    Notes
    -----
    Returns immediately if neither ``rawdata`` nor ``workspace`` exists for the date.
    """
    if caldb_dir is None:
        caldb_dir = default_caldb_dir()

    if slack and slack_dry_run:
        logger.warning('Both slack and slack_dry_run set; only Slack dry-run will run.')
        slack = False

    refresh_astropy_cache()

    lco = LCOGT()

    datedir = Time(date).datetime.strftime('%Y%m%d')
    nested_raw = os.path.join(outdir, datedir, 'rawdata')
    nested_work = os.path.join(outdir, datedir, 'workspace')
    flat_raw = os.path.join(outdir, 'rawdata')
    flat_work = os.path.join(outdir, 'workspace')

    nested_ok = os.path.isdir(nested_raw) or os.path.isdir(nested_work)
    flat_ok = os.path.isdir(flat_raw) or os.path.isdir(flat_work)

    if nested_ok and flat_ok:
        logger.warning(
            'Both %s and %s exist; using nested dated layout for %s',
            flat_raw,
            nested_raw,
            datedir,
        )
    if nested_ok and not flat_ok:
        fulloutdir, fullworkdir = nested_raw, nested_work
    elif flat_ok and not nested_ok:
        fulloutdir, fullworkdir = flat_raw, flat_work
    elif nested_ok and flat_ok:
        fulloutdir, fullworkdir = nested_raw, nested_work
    else:
        fulloutdir, fullworkdir = nested_raw, nested_work

    fullspecdir = os.path.join(outdir, 'spectra')
    fullslackdir = os.path.join(outdir, 'spectra', 'slack')

    if not os.path.exists(fulloutdir) and not os.path.exists(fullworkdir):
        logger.info('No data to process for %s', date)
        return

    os.makedirs(fullspecdir, exist_ok=True)
    os.makedirs(fullslackdir, exist_ok=True)

    files = _list_rawdata_frames(fulloutdir)
    if not files:
        logger.info('No *.fz or *.fits frames in %s for %s', fulloutdir, date)
        return

    if frame_extension is None:
        frame_extension = _infer_frame_extension(fulloutdir)
        logger.info('PypeIt frame extension %s (inferred from rawdata filenames)', frame_extension)
    else:
        logger.info('PypeIt frame extension %s (from --frame-extension)', frame_extension)

    all_data: dict[str, list[dict[str, Any]]] = {}
    spec = 'soar_goodman_red'

    for path in files:
        with fits.open(path, mode='update') as hdu:
            hi = _goodman_frame_header_hdu_index(hdu, path)
            hdr = hdu[hi].header
            prev_obstype = _goodman_normalize_header_for_pypeit_frametyping(hdr)
            if prev_obstype is not None:
                logger.info(
                    'Goodman PypeIt frametyping: OBSTYPE %s -> SPECTRUM in %s',
                    prev_obstype,
                    os.path.basename(path),
                )
            wavemode = hdr['WAVMODE']
            obstype = str(hdr.get('OBSTYPE', '')).strip()
            objname = lco.sanitize_objname(hdr['OBJECT'])
            hdr['OBJECT'] = objname
            mjd = Time(hdr['DATE-OBS']).mjd

            if hdr['INSTRUME'].upper() == 'GHTS_BLUE':
                spec = 'soar_goodman_blue'

            slit = _goodman_slit_from_header(hdr)
            if slit is None:
                logger.warning(
                    'No SLIT/DECKER in %s; skipping frame',
                    os.path.basename(path),
                )
                continue

            data = {
                'wavemod': wavemode,
                'obstype': obstype,
                'filename': path,
                'object': objname,
                'mjd': mjd,
                'slit': slit,
            }

            if wavemode not in {'400_M1', '400_M2'}:
                continue

            all_data.setdefault(wavemode, []).append(data)

    reduction_dirs = _build_reduction_dirs(fullworkdir, all_data)

    jobs = max(1, parallel_jobs)

    # Goodman upstream only allows ``.fz``; keep PASSTA shim active for setup and
    # ``run_pypeit`` (which re-checks raw paths from the .pypeit file).
    with pypeit_invoke.goodman_raw_fits_runtime(spec, frame_extension):
        _parallel_each(
            jobs,
            reduction_dirs,
            lambda od, sp=spec, ext=frame_extension, ow=overwrite: _pypeit_setup_objdir(
                sp, od, ext, overwrite=ow
            ),
        )

        raw_obj = _parallel_map_unordered(
            jobs,
            reduction_dirs,
            lambda od, sp=spec, cdb=caldb_dir, spos=slit_pos, ow=overwrite: _pypeit_reduce_objdir(
                sp, od, cdb, slit_pos=spos, overwrite=ow
            ),
        )
    # Each reduce returns (objdata, force_rebuild_coadd).
    all_objdata: list[tuple[list[tuple[Any, ...]], bool]] = [
        (od, force) for od, force in raw_obj if len(od) > 0
    ]

    outlinks: list[str | None] = []
    for objdata, force_rebuild in all_objdata:
        outlinks.append(
            pypeit_post_process(
                objdata,
                fullspecdir,
                caldb_dir,
                try_old=try_old,
                overwrite=overwrite,
                force_rebuild_coadd=force_rebuild,
            )
        )

    links_done = [L for L in outlinks if L and os.path.exists(L)]
    if snid and links_done:
        pairs = _parallel_map_unordered(
            jobs,
            links_done,
            lambda link, ow=overwrite: (link, run_snid_sage(link, overwrite=ow)),
        )
        if slack_dry_run or slack:
            by_link = dict(pairs)
            dry = bool(slack_dry_run)
            for link in links_done:
                results = by_link.get(link)
                if results is not None:
                    post_slack_results(link, results, fullslackdir, dry_run=dry)


def main(argv: Sequence[str] | None = None) -> None:
    """Console entry point."""
    args = parse_arguments(argv)
    configure_from_parsed_args(args)
    caldb_dir = args.caldb if args.caldb is not None else default_caldb_dir()
    if args.all_to_os is not None:
        jobs = discover_to_os_reduction_jobs(args.all_to_os)
        if not jobs:
            logger.info(
                'No ToOs targets with rawdata or workspace under %s',
                os.path.join(args.all_to_os, 'ToOs'),
            )
            return
        for date_label, od in jobs:
            logger.info('ToOs batch: date=%s outdir=%s', date_label, od)
            run_reduction(
                date_label,
                od,
                caldb_dir=caldb_dir,
                snid=args.snid,
                slack=args.slack,
                slack_dry_run=args.slack_dry_run,
                parallel_jobs=args.jobs,
                frame_extension=args.frame_extension,
                slit_pos=args.slit_pos,
                try_old=args.try_old,
                overwrite=args.overwrite,
            )
    else:
        run_reduction(
            args.date,
            args.outdir,
            caldb_dir=caldb_dir,
            snid=args.snid,
            slack=args.slack,
            slack_dry_run=args.slack_dry_run,
            parallel_jobs=args.jobs,
            frame_extension=args.frame_extension,
            slit_pos=args.slit_pos,
            try_old=args.try_old,
            overwrite=args.overwrite,
        )


if __name__ == '__main__':
    main()
