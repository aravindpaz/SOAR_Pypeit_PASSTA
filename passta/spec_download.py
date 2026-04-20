"""Download SOAR Goodman PASSTA spectra and calibrations from the LCO archive."""

from __future__ import annotations

import argparse
import glob
import logging
import os
import shutil
from typing import Sequence

from astropy import units as u
from astropy.time import Time, TimeDelta

from passta.logging_config import (
    add_logging_arguments,
    configure_from_parsed_args,
)
from passta.lcogt import LCOGT
from passta.versioning import add_version_argument

logger = logging.getLogger(__name__)

_BAD_BASE = ('cfzst', 'wecfzst')


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
    b = basename.lower()
    return b.startswith(_BAD_BASE)


def _filter_bad_basenames(frames: list[dict]) -> list[dict]:
    return [f for f in frames if not _bad_archive_name(f['basename'])]


def _inst_keys(callist: list[dict], obs_by_pid: dict[str, list[dict]]) -> set[str]:
    keys = {o['INSTRUME'].lower() for o in callist}
    for obslist in obs_by_pid.values():
        keys.update(o['INSTRUME'].lower() for o in obslist)
    return keys


def _rmtree_date_if_unused(fulloutdir: str, fullworkdir: str, outdir: str, datedir: str) -> None:
    if not glob.glob(os.path.join(fulloutdir, '*')) and not glob.glob(os.path.join(fullworkdir, '*')):
        shutil.rmtree(os.path.join(outdir, datedir), ignore_errors=True)


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
        return

    lco = LCOGT()
    sdate = Time(date) + TimeDelta(14 * 3600 * u.s)
    edate = sdate + TimeDelta(24 * 3600 * u.s)
    datedir = Time(date).datetime.strftime('%Y%m%d')
    fulloutdir = os.path.join(outdir, datedir, 'rawdata')
    fullworkdir = os.path.join(outdir, datedir, 'workspace')
    os.makedirs(fulloutdir, exist_ok=True)
    os.makedirs(fullworkdir, exist_ok=True)

    telid = '4m0a'
    q = {'sdate': sdate, 'edate': edate, 'telid': telid}

    callist = _filter_bad_basenames(
        lco.get_obslist(**q, propid=['calibrate'], rlevel=[25, 26], obstype='SPECTRUM')
    )
    flatlist = lco.get_obslist(**q, propid=['calibrate'], rlevel=[25, 26], obstype='LAMPFLAT')
    arcs_cal = _filter_bad_basenames(
        lco.get_obslist(**q, propid=['calibrate'], rlevel=[0], obstype='ARC')
    )

    obs_by_pid: dict[str, list[dict]] = {}
    for pid in pids:
        obslist = lco.get_obslist(**q, propid=[pid], rlevel=[25, 26], obstype='SPECTRUM')
        good = [o for o in obslist if not o['filename'].startswith('cfzst')]
        if good:
            obs_by_pid[pid] = list(good)

    if not obs_by_pid:
        _rmtree_date_if_unused(fulloutdir, fullworkdir, outdir, datedir)
        return

    inst = _inst_keys(callist, obs_by_pid)
    flatlist = [f for f in flatlist if f['INSTRUME'].lower() in inst]
    arcs_cal = [f for f in arcs_cal if f['INSTRUME'].lower() in inst]

    dw = max(1, download_workers)
    dl_kw = dict(outrootdir=fulloutdir, skip_header=True, funpack=False, max_download_workers=dw)

    logger.info('Downloading shared calibrators: %s standards, %s flats', len(callist), len(flatlist))
    for lst in (callist, flatlist):
        lco.download_obslist(lst, **dl_kw)

    for pid, obslist in obs_by_pid.items():
        arcs_prop = [
            f for f in lco.get_obslist(**q, propid=[pid], rlevel=[0], obstype='ARC')
            if f['INSTRUME'].lower() in inst and not _bad_archive_name(f['basename'])
        ]
        arclist = _dedupe_frames_by_id(arcs_prop + arcs_cal)
        logger.info('[%s] Need to download %s science, %s arcs', pid, len(obslist), len(arclist))
        for lst in (obslist, arclist):
            lco.download_obslist(lst, **dl_kw)

    _rmtree_date_if_unused(fulloutdir, fullworkdir, outdir, datedir)


def main(argv: Sequence[str] | None = None) -> None:
    """Console entry point."""
    args = parse_arguments(argv)
    configure_from_parsed_args(args)
    download_night(
        args.date,
        args.propid,
        outdir=args.outdir,
        download_workers=args.download_workers,
    )


if __name__ == '__main__':
    main()
