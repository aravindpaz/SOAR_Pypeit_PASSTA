"""Las Cumbres Observatory (LCO) portal client.

Provides :class:`LCOGT` for querying the LCO archive and observation-request API,
downloading frames, and building API v3 observation requests.

Version history
---------------
**1.00** (2019-05-01): Initial code.

**1.01** (2019-09-20): Updated to LCO API v3.

**1.02** (2025-08-25): Added NEWFIRM.

Environment variables
-----------------------
``LCOGT_USERNAME`` and ``LCOGT_PASSWORD`` are required to authenticate.

Author
------
C. D. Kilpatrick
"""

from __future__ import annotations

import copy
import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
import shutil
import time
import warnings
from datetime import datetime, timedelta

import numpy as np
import requests
from astropy import units as u
from astropy.io import fits
from astropy.time import Time, TimeDelta

warnings.filterwarnings("ignore")

logger = logging.getLogger(__name__)


def _archive_download_to_path(url: str, dest_path: str, max_tries: int, verbose: bool) -> None:
    """GET ``url`` to ``dest_path``; retry on connection errors; log failures."""
    if verbose:
        logger.debug('Downloading from: %s', url)
    logger.info('Downloading LCOGT file: %s', dest_path)
    response = None
    for _ in range(max_tries):
        try:
            response = requests.get(url)
            break
        except requests.exceptions.ConnectionError:
            time.sleep(30)
    if response is not None and response.status_code == 200:
        try:
            with open(dest_path, 'wb') as f:
                f.write(response.content)
        except OSError:
            logger.error('Failed to download: %s', dest_path)


class LCOGT:
    """Client for LCO archive queries, downloads, and observation requests."""

    def __init__(self):
        """Read ``LCOGT_*`` credentials and initialize archive/request endpoints."""
        uri_base = 'https://observe.lco.global/'
        archive_base = 'https://archive-api.lco.global/'
        self.params = {
            'uri': {
                'archive': {
                    'authorization': archive_base + 'api-token-auth/',
                    'frames': archive_base + 'frames/',
                    'token': None
                },
                'request': {
                    'authorization': uri_base + 'api/api-token-auth/',
                    'request': uri_base + 'api/requestgroups/',
                    'token': None
                }
            },
            'date_format': '%Y-%m-%d %H:%M:%S',
            'constraints': {
                'max_airmass': 2.5,
                'min_lunar_distance': 15
            },
            'strategy': {
                'default': {
                    'type': 'default',
                    'proposal': [
                    {'name':'default','obstype':'NORMAL'},
                    ],
                    'filters': ['up', 'gp', 'rp', 'ip'],
                    'min_exposure': {'default':45, 'up':150},
                    'max_exposure': 540,
                    # SNR strategy are pairwise mag, snr values.
                    # If current source magnitude is < mag, then use given snr.
                    'snr': [[16,40],[18,20],[99,10]],
                    'cadence': 6,
                    'telescope_class': '1m0',
                    'instrument_type': '1M0-SCICAM-SINISTRO',
                    'acquisition_config': 'OFF',
                    'guiding_config': 'ON',
                    'window': 1.0,
                    'ipp': 1.0
                },
                'fast': {
                    'type': 'photometry',
                    'proposal': [
                    {'name':'default','obstype':'NORMAL'},
                    ],
                    'filters': ['up', 'gp', 'rp', 'ip'],
                    'min_exposure': {'default':45, 'up':150},
                    'max_exposure': 540,
                    # SNR strategy are pairwise mag, snr values.
                    # If current source magnitude is < mag, then use given snr.
                    'snr': [[16,40],[18,20],[99,10]],
                    'cadence': 4,
                    'telescope_class': '1m0',
                    'instrument_type': '1M0-SCICAM-SINISTRO',
                    'acquisition_config': 'OFF',
                    'guiding_config': 'ON',
                    'window': 1.0,
                    'ipp': 1.0
                },
                'spectroscopy': {
                    'type': 'spectroscopy',
                    'proposal': [
                    {'name':'default','obstype':'NORMAL'},
                    ],
                    'min_exposure': 300,
                    'max_exposure': 2400,
                    'telescope_class': '2m0',
                    'instrument_type': '2M0-FLOYDS-SCICAM',
                    'slit': 'slit_1.6as',
                    'acquisition_config': 'OFF',
                    'guiding_config': 'ON',
                    'window': 1.0,
                    'ipp': 1.0
                },
                'photometry': {
                    'type': 'photometry',
                    'proposal': [
                    {'name':'default','obstype':'NORMAL'},
                    ],
                    'filters': ['up', 'gp', 'rp', 'ip'],
                    'min_exposure': {'default':45, 'up':150},
                    'max_exposure': 540,
                    # SNR strategy are pairwise mag, snr values.  first is mag
                    # and second is snr.  If mag_source < mag, then use snr.
                    'snr': [[16,40],[18,20],[99,10]],
                    'cadence': 6,
                    'telescope_class': '1m0',
                    'instrument_type': '1M0-SCICAM-SINISTRO',
                    'acquisition_config': 'OFF',
                    'guiding_config': 'ON',
                    'window': 1.0,
                    'ipp': 1.01
                },
                'template': {
                    'type': 'photometry',
                    'proposal': [
                    {'name':'default','obstype':'NORMAL'},
                    ],
                    'filters': ['up', 'gp', 'rp', 'ip','U','B','V','R','I','zs'],
                    'min_exposure': {'default':600, 'up':600},
                    'max_exposure': 600,
                    # SNR strategy are pairwise mag, snr values.  first is mag
                    # and second is snr.  If mag_source < mag, then use snr.
                    'snr': [[16,40],[18,20],[99,10]],
                    'cadence': 0,
                    'telescope_class': '1m0',
                    'instrument_type': '1M0-SCICAM-SINISTRO',
                    'acquisition_config': 'OFF',
                    'guiding_config': 'ON',
                    'window': 1.0,
                    'ipp': 1.2
                },
                'newfirm-survey': {
                    'type': 'newfirm',
                    'proposal': [
                    {'name':'default','obstype':'NORMAL'},
                    ],
                    'filters': ['JX', 'HX', 'KXs'],
                    'min_exposure': {'default':20},
                    'max_exposure': 600,
                    # SNR strategy are pairwise mag, snr values.  first is mag
                    # and second is snr.  If mag_source < mag, then use snr.
                    'snr': [[16,40],[18,20],[99,10]],
                    'cadence': 0,
                    'telescope_class': '4m0',
                    'instrument_type': 'BLANCO_NEWFIRM',
                    'acquisition_config': 'MANUAL',
                    'guiding_config': 'ON',
                    'window': 1.0,
                    'ipp': 1.05,
                },
                'photometry-spectral': {
                    'type': 'photometry-spectral',
                    'proposal': [
                    {'name':'default','obstype':'NORMAL'}
                    ],
                    'filters': ['up', 'gp', 'rp', 'ip'],
                    'min_exposure': {'default':45, 'up':150},
                    'max_exposure': 600,
                    # SNR strategy are pairwise mag, snr values.  first is mag
                    # and second is snr.  If mag_source < mag, then use snr.
                    'snr': [[16,40],[18,20],[99,10]],
                    'cadence': 4,
                    'telescope_class': '2m0',
                    'instrument_type': '2M0-SCICAM-SPECTRAL',
                    'acquisition_config': 'OFF',
                    'guiding_config': 'ON',
                    'window': 1.0,
                    'ipp': 1.01
                },
                'photometry-muscat': {
                    'type': 'photometry-muscat',
                    'proposal': [
                    {'name':'default','obstype':'NORMAL'}
                    ],
                    'filters': ['all'],
                    'min_exposure': {'default':45, 'up':150},
                    'max_exposure': 600,
                    # SNR strategy are pairwise mag, snr values.  first is mag
                    # and second is snr.  If mag_source < mag, then use snr.
                    'snr': [[16,40],[18,20],[99,10]],
                    'cadence': 4,
                    'telescope_class': '2m0',
                    'instrument_type': '2M0-SCICAM-MUSCAT',
                    'acquisition_config': 'OFF',
                    'guiding_config': 'ON',
                    'window': 1.0,
                    'ipp': 1.01
                },
                'photometry-frb-time-critical': {
                    'type': 'photometry-frb-time-critical',
                    'proposal': [
                    {'name':'default','obstype':'TIME_CRITICAL'}
                    ],
                    'filters': ['rp'],
                    'min_exposure': {'default':60},
                    'max_exposure': 60,
                    # SNR strategy are pairwise mag, snr values.  first is mag
                    # and second is snr.  If mag_source < mag, then use snr.
                    'cadence': 10,
                    'telescope_class': '1m0',
                    'instrument_type': '1M0-SCICAM-SINISTRO',
                    'acquisition_config': 'OFF',
                    'guiding_config': 'ON',
                    'window': 15.0/(60*24),
                    'ipp': 1.2
                },
            }
        }

        self.constants = {
            'zpt': {
                'up': 20.5665,
                'gp': 23.2249,
                'rp': 23.1314,
                'ip': 22.8465,
                'all': 23.2249
            }
        }

        if 'LCOGT_USERNAME' not in os.environ or 'LCOGT_PASSWORD' not in os.environ:
            raise RuntimeError(
                'LCOGT_USERNAME and LCOGT_PASSWORD must be set in the environment.'
            )

        self.username = os.environ['LCOGT_USERNAME']
        self.password = os.environ['LCOGT_PASSWORD']

        self.format = {'token': 'Token {token}'}

    def get_username_password(self):
        """Return stored LCO credentials if both are set.

        Returns
        -------
        tuple[str | None, str | None]
            ``(username, password)`` or ``(None, None)`` if either is missing.
        """
        if self.username is not None and self.password is not None:
            return (self.username, self.password)
        return (None, None)

    def get_token_header(self, username, password, auth_type='archive',
        max_tries=4):
        """Build an ``Authorization`` header for the archive or request API.

        Parameters
        ----------
        username : str or None
            LCO username; if None, no header is produced.
        password : str or None
            LCO password.
        auth_type : {'archive', 'request'}, optional
            Which API token to use (default ``'archive'``).
        max_tries : int, optional
            Retries on connection errors when fetching a new token (default 4).

        Returns
        -------
        dict or None
            Header mapping ``{'Authorization': 'Token ...'}``, or None if
            authentication fails or ``auth_type`` is invalid.
        """
        # Check that we're getting the right token
        if auth_type not in self.params['uri'].keys():
            return(None)

        params = self.params['uri'][auth_type]

        # Check if header token has already been defined
        if params['token']:
            # Return the heade that we need
            fmt = self.format['token']
            header = {'Authorization': fmt.format(token=params['token'])}
            return(header)
        else:
            # Need to generate a new token
            data = {'username': username, 'password': password}
            uri = params['authorization']

            # Now run a request
            tries = 0
            response = None
            while tries < max_tries:
                try:
                    response = requests.post(uri, data=data).json()
                    break
                except requests.exceptions.ConnectionError:
                    time.sleep(30)
                    tries += 1

            # Check that the request worked
            if response and 'token' in response.keys():
                self.params['uri'][auth_type]['token'] = response['token']
                fmt = self.format['token']
                header = {'Authorization': fmt.format(token=response['token'])}
                return(header)
            else:
                # There was some problem with authentication/connection
                return None

    def _archive_frames_request(self, params, headers, max_tries=5):
        """GET the archive ``frames`` endpoint with retries on connection errors.

        Parameters
        ----------
        params : dict
            Query parameters for ``/frames/``.
        headers : dict or None
            Request headers including authorization.
        max_tries : int, optional
            Maximum attempts (default 5).

        Returns
        -------
        requests.Response or None
            HTTP response, or None if every attempt raised ``ConnectionError``.
        """
        url = self.params['uri']['archive']['frames']
        for i in range(max_tries):
            try:
                return requests.get(url, params=params, headers=headers)
            except requests.exceptions.ConnectionError:
                logger.warning('Timeout on archive frames request, try=%s.', i + 1)
                if i + 1 < max_tries:
                    logger.info('Sleeping 30s before retry...')
                    time.sleep(30)
        return None

    def get_spectral_calibrations(self, date, telid, site, outrootdir='',
        funpack=True, max_tries=5):
        """Download spectral calibration frames near a date for a telescope/site.

        Queries LAMPFLAT and ARC frames in a ±1 day window, then downloads them.

        Parameters
        ----------
        date : `astropy.time.Time`
            Reference date.
        telid : str
            Telescope identifier (e.g. ``'4m0a'``).
        site : str
            Site identifier passed to the archive API.
        outrootdir : str, optional
            Output directory for downloads (default current directory).
        funpack : bool, optional
            Passed to :meth:`download_obslist` (default True).
        max_tries : int, optional
            Retries per archive GET (default 5).

        Notes
        -----
        Requires valid LCO credentials and archive authorization.
        """
        username, password = self.get_username_password()
        headers = self.get_token_header(username, password)

        params = {'limit': 100, 'TELID': telid, 'SITEID': site}
        delta = TimeDelta(1, format='jd')
        fmt = self.params['date_format']
        params['start'] = (date - delta).datetime.strftime(fmt)
        params['end'] = (date + delta).datetime.strftime(fmt)

        results = []
        for obstype in ('LAMPFLAT', 'ARC'):
            params['OBSTYPE'] = obstype
            response = self._archive_frames_request(params, headers, max_tries=max_tries)
            if response is None:
                continue
            if response.status_code != 200:
                logger.warning('Archive HTTP error body: %s', response.text[:4000])
            else:
                results += response.json()['results']

        self.download_obslist(
            results,
            outrootdir=outrootdir,
            use_basename=True,
            skip_header=True,
            funpack=funpack,
        )

    @staticmethod
    def _frame_row_matches(row, telid, obstype, rlevel):
        """Return True if an archive frame row passes telid/obstype/rlevel filters."""
        if telid is not None and row['TELID'] not in telid:
            return False
        if obstype is not None and row['configuration_type'] != obstype:
            return False
        if rlevel is None:
            return True
        if isinstance(rlevel, list):
            return row['RLEVEL'] in rlevel
        return row['RLEVEL'] == rlevel

    def _extend_filtered_results(self, results, payload, telid, obstype, rlevel):
        """Append filtered frame dicts from a decoded JSON ``payload`` to ``results``."""
        for row in payload.get('results', []):
            if self._frame_row_matches(row, telid, obstype, rlevel):
                results.append(row)

    def get_obslist(self, propid=None, sdate=None, edate=None, telid=None,
        obstype=None, rlevel=91, obj=None, reqnum=None, ra=None, dec=None,
        public=False, max_tries=5):
        """Query the LCO archive for frames matching optional constraints.

        Parameters
        ----------
        propid : sequence of str or str or None
            Proposal ID(s). A single string is treated as one ID. If omitted,
            all proposals matching other filters are returned (subject to API
            limits).
        sdate, edate : `astropy.time.Time` or None
            Optional start/end window (UTC, formatted per LCO API).
        telid : str, sequence, or None
            If set, keep only rows whose ``TELID`` is contained in this value
            (same membership semantics as the original implementation).
        obstype : str or None
            ``OBSTYPE`` query parameter and/or filter on ``configuration_type``.
        rlevel : int, sequence of int, or None
            Reduction level filter; if a list, any listed level matches.
        obj : str or None
            Object name filter (``OBJECT`` query parameter).
        reqnum : str or None
            Request number filter.
        public : bool, optional
            If True, request public data only (default False).
        max_tries : int, optional
            Retries per HTTP GET on connection errors (default 5).

        Returns
        -------
        list of dict
            Archive frame metadata records passing the filters.
        """
        username, password = self.get_username_password()
        headers = self.get_token_header(username, password)

        params = {'limit': 5000}
        fmt = self.params['date_format']
        results = []

        if propid is None:
            propids = []
        elif isinstance(propid, (str, bytes)):
            propids = [propid]
        else:
            propids = list(propid)

        if sdate is not None:
            params['start'] = sdate.datetime.strftime(fmt)
        if edate is not None:
            params['end'] = edate.datetime.strftime(fmt)
        if obstype is not None:
            params['OBSTYPE'] = obstype
        if obj is not None:
            params['OBJECT'] = obj
        if reqnum is not None:
            params['REQNUM'] = reqnum
        if public:
            params['public'] = True
        if ra is not None and dec is not None:
            params['covers'] = 'POINT({0} {1})'.format(ra, dec)

        if propids:
            for pid in propids:
                q = dict(params)
                q['PROPID'] = pid
                response = self._archive_frames_request(q, headers, max_tries=max_tries)
                if response is None:
                    continue
                if response.status_code != 200:
                    logger.warning('Archive HTTP error body: %s', response.text[:4000])
                    continue
                self._extend_filtered_results(results, response.json(), telid, obstype, rlevel)
        else:
            response = self._archive_frames_request(params, headers, max_tries=max_tries)
            if response is None:
                return results
            if response.status_code != 200:
                logger.warning('Archive HTTP error body: %s', response.text[:4000])
            else:
                self._extend_filtered_results(results, response.json(), telid, obstype, rlevel)

        return results

    def get_standardobs(self, sdate=None, telid=None, rlevel=None,
        max_tries=5):
        """Fetch recent public spectrophotometric standard spectra from the archive.

        Parameters
        ----------
        sdate : `astropy.time.Time` or None
            If None, defaults to 14 days before "now" (with a UTC offset applied
            as in the original implementation).
        telid : optional
            Reserved for compatibility; not applied to results in this method.
        rlevel : int or None
            ``RLEVEL`` filter; if falsy, uses ``0``.
        max_tries : int, optional
            Retries per archive GET (default 5).

        Returns
        -------
        list of dict
            Concatenated archive results for a fixed list of standard-star names.

        Notes
        -----
        ``telid`` is accepted for API compatibility but frames are not re-filtered
        by telescope here.
        """
        username, password = self.get_username_password()
        headers = self.get_token_header(username, password)

        params = {'OBSTYPE': 'SPECTRUM', 'public': True}

        delta = TimeDelta(14, format='jd')
        date = Time(datetime.utcnow()) + TimeDelta(7 * 3600 * u.s)
        start = (date - delta).datetime
        fmt = self.params['date_format']

        if not sdate:
            params['start'] = start.strftime(fmt)
        else:
            params['start'] = sdate.datetime.strftime(fmt)

        if not rlevel:
            params['RLEVEL'] = 0
        else:
            params['RLEVEL'] = rlevel

        objs = ['EGGR274', 'L745-46A', 'FEIGE110', 'HZ 44', 'BD+284211', 'GD71']

        results = []
        for obj in objs:
            params['OBJECT'] = obj
            response = self._archive_frames_request(params, headers, max_tries=max_tries)
            if response is None:
                continue
            if response.status_code != 200:
                logger.warning('Archive HTTP error body: %s', response.text[:4000])
            else:
                results += response.json()['results']

        return results


    def get_requestgroups(self, propid=None, sdate=None, edate=None,
        obstype=None, itype=None):
        """Query submitted observation request groups from the LCO request API.

        Parameters
        ----------
        propid : str, sequence of str, or None
            Proposal id(s); a bare string is treated as a single id.
        sdate, edate : `astropy.time.Time` or None
            Optional ``created_after`` / ``created_before`` filters.
        obstype : optional
            Accepted for API compatibility; unused in the current implementation.
        itype : str or sequence of str or None
            If set, only results whose first configuration matches one of these
            ``instrument_type`` values are kept.

        Returns
        -------
        list of dict
            Aggregated ``results`` entries from the API.
        """
        _ = obstype
        username, password = self.get_username_password()

        # Get authorization token
        headers = self.get_token_header(username, password, auth_type='request')

        if itype is not None:
            if ',' in itype:
                itype = itype.split(',')
            else:
                itype = [itype]

        params = {}
        results = []
        if propid is None:
            propids = []
        elif isinstance(propid, (str, bytes)):
            propids = [propid]
        else:
            propids = list(propid)
        if sdate is not None:
            start = sdate.datetime.strftime(self.params['date_format'])
            params['created_after'] = start
        if edate is not None:
            end = edate.datetime.strftime(self.params['date_format'])
            params['created_before'] = end

        for pid in propids:
            params['proposal'] = pid
            params['limit'] = 500

            # Now do request
            response = requests.get(self.params['uri']['request']['request'],
                params=params, headers=headers).json()

            # If instrument_type, remove values that do not conform
            if itype is not None:
                for r in copy.copy(response['results']):
                    test = r['requests'][0]['configurations'][0]['instrument_type']
                    if test not in itype:
                        response['results'].remove(r)

            results += response['results']

        return(results)

    def sanitize_objname(self, target):
        """Normalize a source name for filesystem-safe strings (PASSTA convention).

        Parameters
        ----------
        target : str
            Raw ``OBJECT`` string.

        Returns
        -------
        str
            Lowercased name with common transient prefixes removed, whitespace
            stripped, and characters ``' '``, ``'('``, ``')'``, ``'/'`` replaced
            or removed.
        """
        target = target.lower()

        if target.startswith('at'):
            target = target[2:]
        elif target.startswith('sn'):
            target = target[2:]

        target = target.strip().replace(' ', '_').translate(str.maketrans('', '', '()/'))

        return target

    def _process_obslist_frame(
        self,
        frame,
        outrootdir,
        use_basename,
        skip_header,
        funpack,
        force_objname,
        verbose,
        max_tries,
    ):
        """Download, funpack, and optionally trim headers for one archive frame."""
        filename = ''
        if use_basename:
            filename = frame['basename'] + '.fits.fz'
            target = self.sanitize_objname(frame['OBJECT'])
        else:
            if force_objname:
                target = force_objname
            else:
                target = self.sanitize_objname(frame['OBJECT'])

            if not target:
                target = 'field' + str(frame['id'])

            filt = frame['FILTER']
            idnum = str(frame['id'])
            date = Time(frame['DATE_OBS']).datetime.strftime('ut%y%m%d')

            if filt.lower() in ['opaque']:
                return
            if filt.lower() == 'no_filter':
                filename = f'{target}.{date}.{idnum}.fits.fz'
            else:
                filename = f'{target}.{date}.{filt}.{idnum}.fits.fz'

        fullfilename = os.path.join(outrootdir, filename)
        unpacked = fullfilename.removesuffix('.fz')
        os.makedirs(outrootdir, exist_ok=True)

        file_exists = False
        url = frame['url']
        if funpack:
            if os.path.exists(unpacked):
                logger.info('LCOGT file: %s already exists!', unpacked)
                file_exists = True
            else:
                _archive_download_to_path(url, fullfilename, max_tries, verbose)
                if os.path.exists(fullfilename):
                    file_exists = True

        elif not os.path.exists(fullfilename):
            _archive_download_to_path(url, fullfilename, max_tries, verbose)
            if os.path.exists(fullfilename):
                file_exists = True
        else:
            logger.info('LCOGT file: %s already exists!', fullfilename)
            file_exists = True

        if not os.path.exists(unpacked) and funpack:
            cmd = 'funpack {file}'
            os.system(cmd.format(file=fullfilename))
        if os.path.exists(fullfilename) and funpack:
            os.remove(fullfilename)

        if not skip_header and not file_exists:
            hdulist = fits.open(unpacked)
            newhdu = fits.HDUList()
            hdu = hdulist['SCI']
            hdu.header['OBSTYPE'] = 'OBJECT'
            hdu.header['OBJECT'] = target
            newhdu.append(hdu)
            newhdu.writeto(unpacked, overwrite=True)

    def download_obslist(self, obslist, outrootdir='', use_basename=False,
        skip_header=False, funpack=True, force_objname='', verbose=False,
        max_tries=4, max_download_workers=1):
        """Download archive frames to disk, optionally funpack and trim headers.

        Parameters
        ----------
        obslist : list of dict
            Archive ``results`` rows (must include ``url``, ``OBJECT``, etc.).
        outrootdir : str, optional
            Output directory (created if missing).
        use_basename : bool, optional
            If True, name outputs ``<basename>.fits.fz`` from the archive record.
        skip_header : bool, optional
            If False, rewrite a single-extension ``SCI`` HDU after funpack.
        funpack : bool, optional
            If True, expand ``.fz`` with ``funpack`` and delete the compressed file.
        force_objname : str, optional
            Override object-based filename stems when ``use_basename`` is False.
        verbose : bool, optional
            Print download URLs.
        max_tries : int, optional
            Retries per file on connection errors.
        max_download_workers : int, optional
            If greater than 1, process frames concurrently (I/O bound HTTP downloads).
            ``funpack`` / FITS post-processing still run inside each worker thread;
            use modest values (for example 4--8) to avoid throttling by the archive.
        """
        if not obslist:
            return
        os.makedirs(outrootdir, exist_ok=True)

        def _one(frame):
            return self._process_obslist_frame(
                frame,
                outrootdir,
                use_basename,
                skip_header,
                funpack,
                force_objname,
                verbose,
                max_tries,
            )

        if max_download_workers <= 1:
            for frame in obslist:
                _one(frame)
            return

        workers = min(max_download_workers, len(obslist))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(_one, frame) for frame in obslist]
            for fut in as_completed(futures):
                fut.result()

    def make_location(self, telescope):
        """Build the ``location`` object (telescope class) for a request group."""
        return {'telescope_class': telescope}

    def make_constraints(self):
        """Return default airmass and lunar separation constraints from ``self.params``."""
        max_airmass = self.params['constraints']['max_airmass']
        min_lunar_distance = self.params['constraints']['min_lunar_distance']
        constraints = {
            'max_airmass': max_airmass,
            'min_lunar_distance': min_lunar_distance
        }
        return constraints

    def make_target(self, name, ra, dec):
        """Return an ICRS ``target`` dict for API v3 (fixed proper motion/parallax)."""
        target = {
            'type': 'ICRS',
            'name': name,
            'ra': ra,
            'dec': dec,
            'proper_motion_ra': 0.0,
            'proper_motion_dec': 0.0,
            'parallax': 0.0,
            'epoch': 2000.0
        }
        return target

    def make_instrument_configs(self, filt, exptime, strat, readout='SLOW'):
        """Return a one-element list of ``instrument_configs`` for the strategy type.

        Parameters
        ----------
        filt : str
            Filter name, ``'spec'`` for spectroscopy templates, or NEWFIRM band.
        exptime : float or str
            Exposure time in seconds (string for NEWFIRM ``exposure_time`` field).
        strat : dict
            Strategy entry from ``self.params['strategy']``.
        readout : {'SLOW', 'FAST'}, optional
            MuSCAT readout mode selector.

        Returns
        -------
        list of dict
            Instrument configuration block consumed by :meth:`make_configurations`.
        """
        if strat['type']=='default':
            configuration = {
                'instrument_name': '1M0-SCICAM-SINISTRO',
                'optical_elements': {'filter': filt},
                'mode': 'full_frame',
                'exposure_time': exptime,
                'exposure_count': 1,
                'bin_x': 1,
                'bin_y': 1,
                'extra_params': {'defocus': 0.0}
            }
        elif strat['type']=='photometry':
            configuration = {
                'instrument_name': '1M0-SCICAM-SINISTRO',
                'optical_elements': {'filter': filt},
                'mode': 'full_frame',
                'exposure_time': exptime,
                'exposure_count': 1,
                'bin_x': 1,
                'bin_y': 1,
                'extra_params': {'defocus': 0.0}
            }
        elif strat['type']=='photometry-frb-time-critical':
            configuration = {
                'instrument_name': '1M0-SCICAM-SINISTRO',
                'optical_elements': {'filter': filt},
                'mode': 'full_frame',
                'exposure_time': exptime,
                'exposure_count': 10,
                'bin_x': 1,
                'bin_y': 1,
                'extra_params': {'defocus': 0.0}
            }
        elif strat['type']=='spectroscopy':
            configuration = {
                "bin_x": 1,
                "bin_y": 1,
                "exposure_count": 1,
                "exposure_time": exptime,
                "mode": "default",
                "rotator_mode": "VFLOAT",
                "extra_params": {},
                "optical_elements": {
                    "slit": strat['slit']
                }
            }
        elif strat['type']=='photometry-muscat':
            if readout=='SLOW':
                mode='MUSCAT_SLOW'
            else:
                mode='MUSCAT_FAST'
            configuration = {
                "exposure_count": 1,
                "exposure_time": exptime,
                "mode": mode,
                "rotator_mode": "",
                "extra_params": {
                    "defocus": 0,
                    "exposure_time_g": exptime,
                    "exposure_time_r": exptime,
                    "exposure_time_i": exptime,
                    "exposure_time_z": exptime,
                    "exposure_mode": "SYNCHRONOUS"
                },
                "optical_elements": {
                    "diffuser_r_position": "out",
                    "diffuser_i_position": "out",
                    "diffuser_g_position": "out",
                    "diffuser_z_position": "out"
                }
            }
        elif strat['type']=='photometry-spectral':
            configuration = {
                'optical_elements': {'filter': filt},
                'mode': 'default',
                'exposure_time': exptime,
                'exposure_count': 1,
                'rotator_mode': '',
                'extra_params': {'defocus': 0.0}
            }
        elif strat['type']=='newfirm':
            # Increase exposure by repeating dither sequence so that all 
            # exposures remain 5s x 4 coadds for calibration purposes
            sequence_repeats = int(np.ceil(exptime/80.0))
            configuration = {
                "exposure_count": 1,
                "exposure_time": "5.0",
                "mode": "fowler1",
                "rotator_mode": "",
                    "extra_params": {
                        "coadds": "4",
                        "sequence_repeats": sequence_repeats,
                        "offset_ra": 0,
                        "offset_dec": 0
                },
                "optical_elements": {
                    "filter": str(filt).lower()
                }
            }
        return [configuration]

    def make_acquisition_config(self, strat, mode=None):
        """Return acquisition configuration; defaults to the strategy template value."""
        if not mode:
            mode = strat['acquisition_config']
        return {'mode': mode}

    def make_guiding_config(self, strat):
        """Return guiding configuration; spectroscopy strategies set ``optional`` False."""
        mode = strat['guiding_config']
        optional = True
        if 'spec' in strat['type']:
            optional=False
        return {'mode': mode, 'optional': optional}

    def make_window(self, start, duration):
        """Define a scheduling window for LCO API v3 (start plus duration in days).

        Parameters
        ----------
        start : datetime.datetime
            Earliest time the observation may execute.
        duration : float
            Window length in days.

        Returns
        -------
        list of dict
            One window dict with LCO-formatted ``start`` and ``end`` strings.
        """
        fmt = self.params['date_format']
        end = start + timedelta(days=duration)
        window = [{
            'start': start.strftime(fmt),
            'end': end.strftime(fmt)
        }]
        return window

    def post_user_request(self, user_request):
        """POST a completed user observation request JSON to the LCO request API.

        Parameters
        ----------
        user_request : dict
            Top-level observation request payload.

        Returns
        -------
        requests.Response
            Raw HTTP response from ``requests.post``.
        """
        username, password = self.get_username_password()

        # Get authorization token
        header = self.get_token_header(username, password, auth_type='request')
        uri = self.params['uri']['request']['request']
        response = requests.post(uri, json=user_request, headers=header)
        return response

    def get_exposure_time(self, filt, mag, strat):
        """Estimate exposure time in seconds from magnitude, filter, and strategy.

        Parameters
        ----------
        filt : str
            Filter name or ``'spec'`` for FLOYDS-style spectroscopy.
        mag : float
            Source magnitude used in the SNR heuristic.
        strat : dict
            Strategy dict containing ``snr`` ladder, exposure limits, etc.

        Returns
        -------
        float or None
            Recommended exposure in seconds, or None if above max exposure or
            spectroscopy rules exclude the source.
        """
        if filt=='spec':
            if mag < 14:
                return(300)
            elif mag < 15:
                return(600)
            elif mag < 16:
                return(900)
            elif mag < 17:
                return(1500)
            elif mag < 17.5:
                return(1800)
            elif mag < 18.0:
                return(2100)
            elif mag < 18.5:
                return(2400)
            else:
                return(None)

        # Handle difference in exposure time for 2m versus 1m
        if 'muscat' in strat['type'] or 'spectral' in strat['type']:
            const = -0.75257498916
        else:
            const = 0.0

        snr = 10.
        if 'snr' in strat.keys():
            for pair in strat['snr']:
                if mag < pair[0]:
                    snr = pair[1]
        term1 = 20. * snr**2
        term2 = 0.4 * (mag - self.constants['zpt'][filt] + const)
        exptime = term1 * 10**term2

        min_exposure = 45.0
        if filt in strat['min_exposure'].keys():
            min_exposure = strat['min_exposure'][filt]
        else:
            min_exposure = strat['min_exposure']['default']
        if exptime < min_exposure:
            if mag < 13.0:
                exptime = min_exposure/2.0
            else:
                exptime = min_exposure

        if exptime > strat['max_exposure']:
            return(None)

        return(exptime)

    def get_max_allowable_ipp(self, obs_requests):
        """Query the API for the maximum allowable IPP for a draft request.

        Parameters
        ----------
        obs_requests : dict
            Observation request JSON as accepted by the ``max_allowable_ipp`` endpoint.

        Returns
        -------
        dict
            Parsed JSON response from the LCO request API.
        """
        username, password = self.get_username_password()

        uri = self.params['uri']['request']['request']+'max_allowable_ipp/'
        header = self.get_token_header(username, password, auth_type='request')

        response = requests.post(uri, json=obs_requests, headers=header)

        return(response.json())

    def make_configurations(self, obj, ra, dec, mag, strat, inst_config=[]):
        """Build LCO API v3 ``configurations`` for one target and strategy.

        Assembles constraints, instrument/acquisition/guiding settings, target
        metadata, and exposure type per the selected ``strat`` template.
        """
        # Iterate over each filter and append to configurations list
        configurations = []

        # If default, construct configurations for photometry strategy
        if (strat['type']=='default' or 'photometry' in strat['type']):
            for i,filt in enumerate(strat['filters']):
                constraints = self.make_constraints()

                iconfig = {}
                for element in inst_config:
                    if 'filter' in element.keys() and element['filter']==filt:
                        iconfig = element
                        break

                # Get exposure time and make instrument_config
                if (iconfig and 'exptime' in iconfig.keys() and
                    'filter' in iconfig.keys()):
                    if filt!=iconfig['filter']: continue
                    exptime = iconfig['exptime']
                    instrument_configs = self.make_instrument_configs(filt,
                        exptime, strat)
                else:
                    exptime = self.get_exposure_time(filt, mag, strat)
                    if not exptime:
                        continue
                    instrument_configs = self.make_instrument_configs(filt,
                        exptime, strat)

                # Make acquisition and guiding config with strat
                acquisition_config = self.make_acquisition_config(strat)
                guiding_config = self.make_guiding_config(strat)

                # Make a tar object
                target = self.make_target(obj, ra, dec)

                # Compile everything into configuration
                configuration = {
                    'constraints': constraints,
                    'instrument_configs': instrument_configs,
                    'acquisition_config': acquisition_config,
                    'guiding_config': guiding_config,
                    'target': target,
                    'instrument_type': strat['instrument_type'],
                    'type': 'EXPOSE',
                    'priority': i+1
                }
                configurations.append(configuration)

        elif strat['type']=='spectroscopy':
            # Need to construct LAMP FLAT, ARC, SPECTRUM, ARC, LAMP FLAT
            for obstype in ['LAMP_FLAT', 'ARC', 'SPECTRUM', 'ARC', 'LAMP_FLAT']:

                # Use default constraint values
                constraints = self.make_constraints()

                exptime = self.get_exposure_time('spec', mag, strat)
                if not exptime:
                    continue

                if obstype=='LAMP_FLAT':
                    exptime = 50
                if obstype=='ARC':
                    exptime = 60

                # Make acquisition and guiding config with strat
                instrument_configs = self.make_instrument_configs('spec', exptime, strat)
                if obstype=='SPECTRUM':
                    acquisition_config = self.make_acquisition_config(strat, mode='WCS')
                else:
                    acquisition_config = self.make_acquisition_config(strat)
                guiding_config = self.make_guiding_config(strat)

                # Make a target object
                target = self.make_target(obj, ra, dec)

                # Compile everything into configuration
                configuration = {
                    'type': obstype,
                    'constraints': constraints,
                    'instrument_configs': instrument_configs,
                    'acquisition_config': acquisition_config,
                    'guiding_config': guiding_config,
                    'target': target,
                    'instrument_type': strat['instrument_type']
                }
                configurations.append(configuration)

        elif strat['type']=='newfirm':

            instrument_configs = []
            for i,filt in enumerate(strat['filters']):

                iconfig = {}
                for element in inst_config:
                    if 'filter' in element.keys() and element['filter']==filt:
                        iconfig = element
                        break

                # Get exposure time and make instrument_config
                if (iconfig and 'exptime' in iconfig.keys() and
                    'filter' in iconfig.keys()):
                    if filt!=iconfig['filter']: continue
                    exptime = iconfig['exptime']
                    instrument_configs.extend(self.make_instrument_configs(filt,
                        exptime, strat))
                else:
                    exptime = self.get_exposure_time(filt, mag, strat)
                    if not exptime:
                        continue
                    instrument_configs.extend(self.make_instrument_configs(filt,
                        exptime, strat))

            # Make acquisition and guiding config with strat
            acquisition_config = self.make_acquisition_config(strat)
            guiding_config = self.make_guiding_config(strat)
            constraints = self.make_constraints()

            # Make a tar object
            target = self.make_target(obj, ra, dec)

            # Compile everything into configuration
            configuration = {
                    'constraints': constraints,
                    'instrument_configs': instrument_configs,
                    'acquisition_config': acquisition_config,
                    'guiding_config': guiding_config,
                    'target': target,
                    'instrument_type': strat['instrument_type'],
                    'type': 'EXPOSE',
                    'priority': i+1,
                    "extra_params": {
                        "dither_value": 80,
                        "dither_sequence": "2x2",
                        "detector_centering": "det_3",
                        "dither_sequence_random_offset": True,
                },
            }
            
            configurations.append(configuration)

        return(configurations)


    def make_requests(self, obj, ra, dec, mag, strat, start=None,
        inst_config={}):
        """Build the ``requests`` list for an LCO API v3 observation group."""
        location = self.make_location(strat['telescope_class'])
        if start:
            dt = start
            dt = dt.datetime
        else:
            dt = Time(datetime.utcnow())+TimeDelta(0.5*3600*u.s)
            dt = dt.datetime

        window = self.make_window(dt, strat['window'])
        configurations = self.make_configurations(obj, ra, dec, mag, strat,
            inst_config=inst_config)

        requests = [{
            'location': location,
            'windows': window,
            'configurations': configurations,
            'optimization_type': "AIRMASS",
            "acceptability_threshold": 90,
            "configuration_repeats": 1,
        }]

        return(requests)

    def make_obs_request(self, obj, ra, dec, mag, strategy = 'default',
        propidx=0, start=None, recalculate_ipp=False, inst_config=[]):
        """Assemble and optionally POST the outermost LCO observation request.

        Wraps proposal metadata, IPP handling, and nested ``requests`` produced by
        :meth:`make_requests`.
        """

        # Get params - strategy data
        strat = self.params['strategy'][strategy]

        if propidx > len(strat['proposal'])-1:
            return(None)

        proposal = strat['proposal'][propidx]

        # Make the obs_request dictionary
        obs_request = {
            'name': obj,
            'proposal': proposal['name'],
            'ipp_value': strat['ipp'],
            'operator': 'SINGLE',
            'observation_type': proposal['obstype']
        }

        # Iterate through the next level of request
        requests = self.make_requests(obj, ra, dec, mag, strat, start=start,
            inst_config=inst_config)
        obs_request['requests']=requests

        if recalculate_ipp and len(requests[0]['configurations'])>0:
            logger.info('Recalculating IPP...')
            response = self.get_max_allowable_ipp(obs_request)
            logger.debug('max_allowable_ipp response: %s', response)
            try:
                if 'errors' in response.keys() and 'requests' in response['errors'].keys():
                    data = response['errors']['requests']
                    if 'non_field_errors' in data[0].keys():
                        logger.error('Invalid obs request: %s', obs_request)
                        return(None)
                elif 'errors' in response.keys():
                    data = response['errors']
                    if 'non_field_errors' in data.keys():
                        logger.error('Invalid obs request: %s', obs_request)
                        return(None)
            except Exception:
                logger.debug('IPP error parsing skipped', exc_info=True)
            max_ipp = None
            for key in response.keys():
                r = response[key]
                for key in r.keys():
                    ipp_data = r[key]
                    logger.debug('IPP data: %s', ipp_data)
                    if 'max_allowable_ipp_value' in ipp_data.keys():
                        max_ipp = ipp_data['max_allowable_ipp_value']

            if max_ipp:
                curr_ipp = obs_request['ipp_value']
                if max_ipp < curr_ipp:
                    logger.info('Adjusting IPP value to: %s', max_ipp)
                    obs_request['ipp_value']=max_ipp
                else:
                    m = 'No need to adjust IPP, {0} < {1}'.format(curr_ipp,
                        max_ipp)
                    logger.info('%s', m)


        if not requests[0]['configurations']:
            return(None)
        else:
            response = self.post_user_request(obs_request)
            response = response.json()
            return(response)
