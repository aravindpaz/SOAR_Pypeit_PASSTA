"""Run PypeIt on archived PASSTA Goodman data: setup, reduce, flux, coadd, telluric."""

from __future__ import annotations

import argparse
import astropy
import glob
import logging
import os
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

from pypeit import inputfiles

from passta import slack_utils
from passta.logging_config import add_logging_arguments, configure_from_parsed_args
from passta.lcogt import LCOGT
from passta.versioning import add_version_argument

logger = logging.getLogger(__name__)


def _run_shell(cmd: str) -> None:
    logger.info('Running: %s', cmd)
    os.system(cmd)


def _fits_stem(filename: str) -> str:
    return filename.replace('.fits', '').replace('.fz', '')


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
    with ThreadPoolExecutor(max_workers=w) as pool:
        futs = [pool.submit(fn, x) for x in items]
        for fut in as_completed(futs):
            fut.result()


def _parallel_map_unordered(jobs: int, items: Sequence, fn):
    n = max(1, jobs)
    if n <= 1 or len(items) <= 1:
        return [fn(x) for x in items]
    w = min(n, len(items))
    with ThreadPoolExecutor(max_workers=w) as pool:
        futs = [pool.submit(fn, x) for x in items]
        return [fut.result() for fut in as_completed(futs)]


def default_caldb_dir() -> str:
    """Return the filesystem path to the packaged sensitivity-function database."""
    return str(resources.files('passta').joinpath('caldb'))


def refresh_astropy_cache() -> None:
    """Force-download Astropy observatory/site data with a few retries."""
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
        'date',
        type=str,
        help='Night label matching the dated reduction tree under ``outdir``.',
    )
    parser.add_argument(
        'outdir',
        type=str,
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
    add_logging_arguments(parser)
    args = parser.parse_args(argv)
    args.outdir = os.path.abspath(args.outdir)
    if args.slack_dry_run:
        args.snid = True
    return args


def _pypeit_setup_objdir(spec: str, objdir: str) -> None:
    """Run ``pypeit_setup`` for one reduction directory (subprocess)."""
    cmd = (
        f'pypeit_setup -s {spec} -c all -r {objdir} -d {objdir} '
        f'> {objdir}/setup.log 2> {objdir}/setup.log'
    )
    _run_shell(cmd)


def _pypeit_reduce_objdir(spec: str, objdir: str, caldb_dir: str) -> list[tuple[Any, ...]]:
    """Run ``run_pypeit`` if needed and return science rows from :func:`parse_pypeit_output`."""
    pypeit_file = os.path.join(objdir, f'{spec}_A/{spec}_A.pypeit')
    science_dir = os.path.join(objdir, 'Science')
    if not os.path.exists(pypeit_file):
        logger.error('pypeit_setup failed or missing .pypeit for %s', objdir)
        return []
    pypeit_table = inputfiles.PypeItFile.from_file(pypeit_file)
    if not os.path.exists(science_dir):
        cmd = (
            f'run_pypeit {pypeit_file} -r {objdir} '
            f'> {objdir}/reduction.log 2> {objdir}/reduction.log'
        )
        _run_shell(cmd)
    else:
        logger.info('Reduction already done for %s', objdir)
    return parse_pypeit_output(science_dir, pypeit_table, caldb_dir)


def parse_pypeit_output(
    scidir: str,
    pypeit_table: Any,
    caldb: str,
    slit_pos: float = 505.0,
    pix_tolerance: Sequence[float] = (15.0, 30.0, 60.0),
) -> list[tuple[Any, ...]]:
    """Inspect PypeIt ``Science`` products and build sensitivity functions.

    For each science row, selects the highest-S/N trace near ``slit_pos``. For
    standards, runs ``pypeit_sensfunc`` and appends metadata to ``caldb.txt``.

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
            senslog = os.path.join(scidir, sensfilename).replace('.fits', '.log')
            par_outfile = os.path.join(scidir, 'sensfunc.par')

            os.makedirs(caldb, exist_ok=True)
            os.makedirs(os.path.join(caldb, 'plots'), exist_ok=True)

            cmd = (
                f'pypeit_sensfunc {outfile} -o {fullsensfilename} '
                f'--par_outfile {par_outfile} > {senslog} 2> {senslog}'
            )
            _run_shell(cmd)

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
        f.write('[fluxcalib] \n')
        f.write('  extinct_correct = False \n\n')
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
) -> str | None:
    """Flux-calibrate, coadd, telluric-correct, and symlink the final spectrum.

    Parameters
    ----------
    objdata : list of tuple
        Science rows for a single object from :func:`parse_pypeit_output`.
    specdir : str
        Directory where the telluric-corrected FITS symlink is created.
    caldb_dir : str
        Location of ``caldb.txt`` and sensitivity functions.

    Returns
    -------
    str or None
        Path to the telluric-corrected spectrum if produced, else None.

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
    flux_logname = os.path.join(dirname, 'flux.log')
    if flux_filename is None:
        return None

    _run_shell(
        f'pypeit_flux_calib --par_outfile {flux_filename} > {flux_logname} 2> {flux_logname}'
    )

    out_1dfilename = f'{target}.{datestr}.{dispname}.coadd1d.fits'
    fullout_1dfilename = os.path.join(dirname, out_1dfilename)
    full_tellcorrfilename = fullout_1dfilename.replace('.fits', '_tellcorr.fits')

    coadd_filename = make_coadd_file(objdata, fullout_1dfilename)
    coadd_logname = os.path.join(dirname, 'coadd.log')
    par_outfile = os.path.join(dirname, 'coadd1d.par')
    if coadd_filename is None:
        return None

    if not os.path.exists(fullout_1dfilename):
        _run_shell(
            f'pypeit_coadd_1dspec {coadd_filename} --par_outfile {par_outfile} '
            f'> {coadd_logname} 2> {coadd_logname}'
        )

    if os.path.exists(fullout_1dfilename):
        tellcorr_logname = os.path.join(dirname, 'tellcorr.log')
        par_outfile = os.path.join(dirname, 'tellfit.par')
        os.chdir(dirname)
        if not os.path.exists(full_tellcorrfilename):
            _run_shell(
                f'pypeit_tellfit {fullout_1dfilename} --objmodel poly '
                f'--par_outfile {par_outfile} > {tellcorr_logname} 2> {tellcorr_logname}'
            )

    if os.path.exists(full_tellcorrfilename):
        outlink = os.path.join(specdir, os.path.basename(full_tellcorrfilename))
        if not os.path.exists(outlink):
            os.symlink(full_tellcorrfilename, outlink)
        if os.path.exists(outlink):
            logger.info('Linked %s -> %s', full_tellcorrfilename, outlink)
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
    if res_table is not None and len(res_table) > 0:
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


def run_snid_sage(filepath: str) -> dict[str, Any]:
    """Run external ``sage`` (SNID SAGE) on a 1D spectrum and parse outputs.

    Parameters
    ----------
    filepath : str
        Telluric-corrected FITS spectrum.

    Returns
    -------
    dict
        Keys ``results`` (`~astropy.table.Table` or None) and ``plot`` (path
        string or None).
    """
    dirpath, basename = os.path.split(filepath)
    results_dir = os.path.join(dirpath, 'results')
    result_file = os.path.join(results_dir, basename.replace('.fits', '.output'))
    spec_plot = os.path.join(
        results_dir, basename.replace('.fits', '_flattened_spectrum.png')
    )

    os.makedirs(results_dir, exist_ok=True)
    logpath = os.path.join(results_dir, basename.replace('.fits', '.log'))

    if not os.path.exists(result_file):
        _run_shell(f'sage {filepath} -o {results_dir} > {logpath} 2> {logpath}')

    out: dict[str, Any] = {'results': None, 'plot': None}

    if os.path.exists(result_file):
        data_start = False
        all_data: list[list[str]] = []
        with open(result_file) as outf:
            for line in outf:
                if line.startswith('--------------------'):
                    data_start = True
                elif not data_start:
                    continue
                else:
                    all_data.append(line.split())

        if all_data and all(len(r) > 0 for r in all_data):
            columns = [list(row) for row in zip(*all_data)]
            out['results'] = Table(
                columns,
                names=('idx', 'name', 'type', 'subtype', 'rlap-ccc', 'redshift', 'error', 'age'),
            )

    if os.path.exists(spec_plot):
        out['plot'] = spec_plot

    return out


def run_reduction(
    date: str,
    outdir: str,
    caldb_dir: str | None = None,
    snid: bool = False,
    slack: bool = False,
    slack_dry_run: bool = False,
    parallel_jobs: int = 1,
) -> None:
    """Main driver: symlink raw data, ``pypeit_setup`` / ``run_pypeit``, post-process.

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

    Notes
    -----
    Exits silently if neither ``rawdata`` nor ``workspace`` exists for the date.
    """
    if caldb_dir is None:
        caldb_dir = default_caldb_dir()

    if slack and slack_dry_run:
        logger.warning('Both slack and slack_dry_run set; only Slack dry-run will run.')
        slack = False

    refresh_astropy_cache()

    lco = LCOGT()

    datedir = Time(date).datetime.strftime('%Y%m%d')
    fulloutdir = os.path.join(outdir, datedir, 'rawdata')
    fullworkdir = os.path.join(outdir, datedir, 'workspace')
    fullspecdir = os.path.join(outdir, 'spectra')
    fullslackdir = os.path.join(outdir, 'spectra', 'slack')

    if not os.path.exists(fulloutdir) and not os.path.exists(fullworkdir):
        logger.info('No data to process for %s', date)
        sys.exit(0)

    os.makedirs(fullspecdir, exist_ok=True)
    os.makedirs(fullslackdir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(fulloutdir, '*.fz')))
    all_data: dict[str, list[dict[str, Any]]] = {}
    spec = 'soar_goodman_red'

    for path in files:
        with fits.open(path, mode='update') as hdu:
            wavemode = hdu[1].header['WAVMODE']
            obstype = hdu[1].header['OBSTYPE']
            objname = lco.sanitize_objname(hdu[1].header['OBJECT'])
            hdu[1].header['OBJECT'] = objname
            mjd = Time(hdu[1].header['DATE-OBS']).mjd

            if hdu[1].header['INSTRUME'].upper() == 'GHTS_BLUE':
                spec = 'soar_goodman_blue'

            data = {
                'wavemod': wavemode,
                'obstype': obstype,
                'filename': path,
                'object': objname,
                'mjd': mjd,
            }

            if wavemode not in {'400_M1', '400_M2'}:
                continue

            all_data.setdefault(wavemode, []).append(data)

    reduction_dirs: list[str] = []

    for wavemode, wavedata in all_data.items():
        wavedir = os.path.join(fullworkdir, wavemode)
        os.makedirs(wavedir, exist_ok=True)

        wavedata = sorted(wavedata, key=lambda x: x['mjd'])
        arcdata = [d for d in wavedata if d['obstype'] in {'ARC', 'COMP'}]
        all_objs = np.unique(
            [
                d['object']
                for d in wavedata
                if d['obstype'] in {'SPECTRUM', 'OBJECT'} and d['object'].strip()
            ]
        )

        for objname in all_objs:
            objdir = os.path.join(wavedir, objname)
            os.makedirs(objdir, exist_ok=True)

            mjds = []
            for data in wavedata:
                filename = data['filename']
                basefile = os.path.basename(filename)
                outfile = os.path.join(objdir, basefile)

                if data['obstype'] in {'SPECTRUM', 'OBJECT'}:
                    mjds.append(data['mjd'])
                    if data['object'] != objname:
                        continue
                if data['obstype'] in {'ARC', 'COMP'}:
                    continue

                if not os.path.exists(outfile):
                    logger.info('Linking %s -> %s', filename, outfile)
                    os.symlink(filename, outfile)

            avg_mjd = float(np.average(mjds))
            arcdata_sorted = sorted(arcdata, key=lambda x: abs(x['mjd'] - avg_mjd))
            arcfile = arcdata_sorted[0]['filename']
            basefile = os.path.basename(arcfile)
            outfile = os.path.join(objdir, basefile)

            if not os.path.exists(outfile):
                logger.info('Linking %s -> %s', arcfile, outfile)
                os.symlink(arcfile, outfile)

            reduction_dirs.append(objdir)

    jobs = max(1, parallel_jobs)

    _parallel_each(jobs, reduction_dirs, lambda od: _pypeit_setup_objdir(spec, od))

    raw_obj = _parallel_map_unordered(
        jobs,
        reduction_dirs,
        lambda od: _pypeit_reduce_objdir(spec, od, caldb_dir),
    )
    all_objdata: list[list[tuple[Any, ...]]] = [od for od in raw_obj if len(od) > 0]

    outlinks: list[str | None] = []
    for objdata in all_objdata:
        outlinks.append(pypeit_post_process(objdata, fullspecdir, caldb_dir))

    links_done = [L for L in outlinks if L and os.path.exists(L)]
    if snid and links_done:
        pairs = _parallel_map_unordered(
            jobs,
            links_done,
            lambda link: (link, run_snid_sage(link)),
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
    run_reduction(
        args.date,
        args.outdir,
        caldb_dir=caldb_dir,
        snid=args.snid,
        slack=args.slack,
        slack_dry_run=args.slack_dry_run,
        parallel_jobs=args.jobs,
    )


if __name__ == '__main__':
    main()
