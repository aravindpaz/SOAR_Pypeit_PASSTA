# v1.00 2019-05-01: basic code
# v1.01 2019-09-20: updated to LCO API v3
# v1.02 2025-08-25: Added NEWFIRM

"""

Code for viewing and posting observation requests and viewing and downloading
data from the LCO portal.

Author: C. D. Kilpatrick

"""

# PACKAGES
import requests
import warnings
import os
import stat
import json
import shutil
import sys
import copy
import time
import shibboleth
from astropy.time import Time, TimeDelta
from astropy.io import ascii, fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from dateutil.parser import parse
from datetime import datetime, timedelta
import numpy as np

import warnings
warnings.filterwarnings('ignore')

class lcogt(object):
    def __init__(self):

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

        if ('LCOGT_USERNAME' not in os.environ.keys()) or ('LCOGT_PASSWORD' not in os.environ.keys()):
            raise Exception('LCOGT_USERNAME and LCOGT_PASSWORD need to be environment variables.')

        self.username = os.environ['LCOGT_USERNAME']
        self.password = os.environ['LCOGT_PASSWORD']

        self.format = {'token': 'Token {token}'}

    def get_username_password(self):
        if (self.username is not None and self.password is not None):
            return(self.username, self.password)
        else:
            # Get username password out of shibboleth
            return(None, None)

    def get_token_header(self, username, password, auth_type='archive',
        max_tries=4):

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
                return(None)

    # Download spectral calibration files for a specific telid and date
    def get_spectral_calibrations(self, date, telid, site, outrootdir='',
        funpack=True, max_tries=5):

        # Get username and password
        username,password = self.get_username_password()

        # Get authorization token
        headers = self.get_token_header(username, password)

        params = {'limit': 100}
        results = []
        params['TELID'] = telid
        params['SITEID'] = site
        delta = TimeDelta(1, format='jd')
        params['start'] = (date-delta).datetime.strftime(self.params['date_format'])
        params['end'] = (date+delta).datetime.strftime(self.params['date_format'])

        # First check for LAMPFLATs
        params['OBSTYPE'] = 'LAMPFLAT'
        i=0
        while True:
            i=i+1
            try:
                response = requests.get(self.params['uri']['archive']['frames'],
                    params=params, headers=headers)
                break
            except requests.exceptions.ConnectionError:
                # Could be the API throttling us - let's retry up to 5 times
                print(f'Timeout on requesting spectral calibrations, try={i}.')
                if i<max_tries:
                    print('Sleep for 30s and try again...')
                    time.sleep(30)
                else:
                    print(f'Reached max tries={max_tries}')
                    break

        if response.status_code != 200:
            print(response.text)
        else:
            data = response.json()
            results += data['results']

        # Now check for ARCs
        params['OBSTYPE'] = 'ARC'
        i=0
        while True:
            i=i+1
            try:
                response = requests.get(self.params['uri']['archive']['frames'],
                    params=params, headers=headers)
                break
            except requests.exceptions.ConnectionError:
                # Could be the API throttling us - let's retry up to 5 times
                print(f'Timeout on requesting spectral calibrations, try={i}.')
                if i<max_tries:
                    print('Sleep for 30s and try again...')
                    time.sleep(30)
                else:
                    print(f'Reached max tries={max_tries}')
                    break

        if response.status_code != 200:
            print(response.text)
        else:
            data = response.json()
            results += data['results']

        # Now try to download the obslist
        self.download_obslist(results, outrootdir=outrootdir,
            use_basename=True, skip_header=True, funpack=funpack)

    # Get a json object with complete list of observations.  Optional to define
    # a program or date range to narrow search
    def get_obslist(self, propid=None, sdate=None, edate=None, telid=None,
        obstype=None, rlevel=91, obj=None, reqnum=None, ra=None, dec=None,
        public=False, max_tries=5):

        # Get username and password
        username,password = self.get_username_password()

        # Get authorization token
        headers = self.get_token_header(username, password)

        # Get parameters for request
        params = {'limit': 5000}
        fmt=self.params['date_format']
        results = []
        propids = []
        #if rlevel is not None:
        #    params['RLEVEL'] = rlevel
        #else:
        #    params['RLEVEL'] = 91
        if propid is not None:
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
            params['public']=True
        if ra and dec:
            params['covers']='POINT({0} {1})'.format(ra, dec)

        if len(propids)>0:
            for pid in propids:
                params['PROPID'] = pid

                # Now do request
                i=0
                response = None
                while True:
                    i=i+1
                    try:
                        response = requests.get(self.params['uri']['archive']['frames'],
                            params=params, headers=headers)
                        break
                    except requests.exceptions.ConnectionError:
                        # Could be the API throttling us - let's retry up to 5 times
                        print(f'Timeout on requesting spectral calibrations, try={i}.')
                        if i<max_tries:
                            print('Sleep for 30s and try again...')
                            time.sleep(30)
                        else:
                            print(f'Reached max tries={max_tries}')
                            break

                if response and response.status_code != 200:
                    print(response.text)
                else:
                    data = response.json()
                    results = []

                    for r in data['results']:
                        if telid is not None and r['TELID'] not in telid:
                            continue
                        if (obstype is not None and r['configuration_type']!=obstype):
                            continue
                        if rlevel is not None:
                            if isinstance(rlevel, list):
                                if r['RLEVEL'] not in rlevel:
                                    continue
                            elif r['RLEVEL']!=rlevel:
                                continue
                            
                        results.append(r)

        else:
            i=0
            response = None
            while True:
                i=i+1
                try:
                    response = requests.get(self.params['uri']['archive']['frames'],
                        params=params, headers=headers)
                    break
                except requests.exceptions.ConnectionError:
                    # Could be the API throttling us - let's retry up to 5 times
                    print(f'Timeout on requesting spectral calibrations, try={i}.')
                    if i<max_tries:
                        print('Sleep for 30s and try again...')
                        time.sleep(30)
                    else:
                        print(f'Reached max tries={max_tries}')
                        break
            if response and response.status_code!=200:
                print(response.text)
            else:
                data = response.json()
                for r in data['results']:
                    if telid is not None and r['TELID'] not in telid:
                        continue
                    if (obstype is not None and r['configuration_type']!=obstype):
                        continue
                    if rlevel is not None:
                        if isinstance(rlevel, list):
                            if r['RLEVEL'] not in rlevel:
                                continue
                        elif r['RLEVEL']!=rlevel:
                            continue
                        
                    results.append(r)

        return(results)

    def get_standardobs(self, sdate=None, telid=None, rlevel=None,
        max_tries=5):

        # Get username and password
        username,password = self.get_username_password()

        # Get authorization token
        headers = self.get_token_header(username, password)

        # Get parameters for request
        params = {'OBSTYPE': 'SPECTRUM', 'public': True}

        delta = TimeDelta(14, format='jd')
        date = Time(datetime.utcnow())+TimeDelta(7*3600*u.s)
        start = (date-delta).datetime
        fmt=self.params['date_format']

        if not sdate: params['start']=start.strftime(fmt)
        else: params['start']=sdate.datetime.strftime(fmt)

        if not rlevel: params['RLEVEL']=0
        else: params['RLEVEL']=rlevel

        objs=['EGGR274','L745-46A','FEIGE110','HZ 44','BD+284211','GD71']

        results=[]
        for obj in objs:
            params['OBJECT']=obj
            i=0
            response = None
            while True:
                i=i+1
                try:
                    response = requests.get(self.params['uri']['archive']['frames'],
                        params=params, headers=headers)
                    break
                except requests.exceptions.ConnectionError:
                    # Could be the API throttling us - let's retry up to 5 times
                    print(f'Timeout on requesting spectral calibrations, try={i}.')
                    if i<max_tries:
                        print('Sleep for 30s and try again...')
                        time.sleep(30)
                    else:
                        print(f'Reached max tries={max_tries}')
                        break

            if response and response.status_code != 200:
                print(response.text)
            else:
                data = response.json()
                results += data['results']

        return(results)


    def get_requestgroups(self, propid=None, sdate=None, edate=None,
        obstype=None, itype=None):

        # Get username and password
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
        propids = []
        if propid is not None:
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

        target = target.lower()

        if target.startswith('at'):
            target = target[2:]
        elif target.startswith('sn'):
            target = target[2:]

        target = target.strip()

        # Need to sanitize target, e.g., for spaces
        target = target.replace(' ','_')
        target = target.replace('(','')
        target = target.replace(')','')
        target = target.replace('/','')

        return(target)


    # Given an input obslist from the LCOGT archive, download the associated
    # files for each element
    def download_obslist(self, obslist, outrootdir='', use_basename=False,
        skip_header=False, funpack=True, force_objname='', verbose=False,
        max_tries=4):

        for frame in obslist:
            filename = ''
            if use_basename:
                filename = frame['basename'] + '.fits.fz'
            else:
                if force_objname:
                    target = force_objname
                else:
                    target = self.sanitize_objname(frame['OBJECT'])

                if not target:
                    target = 'field'+str(frame['id'])

                filt = frame['FILTER']
                idnum = str(frame['id'])
                date = Time(frame['DATE_OBS']).datetime.strftime('ut%y%m%d')

                if filt.lower() in ['opaque']: 
                    continue
                elif filt.lower()=='no_filter':
                    filename = f'{target}.{date}.{idnum}.fits.fz'
                else:
                    filename = f'{target}.{date}.{filt}.{idnum}.fits.fz'

            fullfilename = os.path.join(outrootdir, filename)
            if not os.path.exists(outrootdir):
                shutil.os.makedirs(outrootdir)

            tries = 0
            file_exists = False
            if funpack:
                if not os.path.exists(fullfilename.strip('.fz')):
                    message = 'Downloading LCOGT file: {file}'
                    print(message.format(file=fullfilename))
                    url = frame['url']
                    if verbose:
                        print(f'Downloading from: {url}')

                    r = None
                    while tries < max_tries:
                        try:
                            r = requests.get(url)
                            break
                        except requests.exceptions.ConnectionError:
                            time.sleep(30)
                            tries += 1

                    if r is not None and r.status_code==200:
                        with open(fullfilename, 'wb') as f:
                            try:
                                f.write(r.content)
                            except:
                                print(f'Failed to download: {fullfilename}')
                    
                    if os.path.exists(fullfilename):
                        file_exists = True
                else:
                    message = 'LCOGT file: {file} already exists!'
                    print(message.format(file=fullfilename.strip('.fz')))
                    file_exists = True

            elif not os.path.exists(fullfilename):
                message = 'Downloading LCOGT file: {file}'
                print(message.format(file=fullfilename))
                url=frame['url']
                if verbose:
                    print(f'Downloading from: {url}')

                r = None
                while tries < max_tries:
                    try:
                        r = requests.get(url)
                        break
                    except requests.exceptions.ConnectionError:
                        time.sleep(30)
                        tries += 1

                if r is not None and r.status_code==200:
                    with open(fullfilename, 'wb') as f:
                        try:
                            f.write(r.content)
                        except:
                            print(f'Failed to download: {fullfilename}')
                    
                if os.path.exists(fullfilename):
                    file_exists = True
            else:
                message = 'LCOGT file: {file} already exists!'
                print(message.format(file=fullfilename.strip('.fz')))
                file_exists = True


            # funpack
            if not os.path.exists(fullfilename.strip('.fz')) and funpack:
                cmd = 'funpack {file}'
                os.system(cmd.format(file=fullfilename))
            if os.path.exists(fullfilename) and funpack:
                os.remove(fullfilename)

            # Remove extraneous extension
            if not skip_header and not file_exists:
                hdulist = fits.open(fullfilename.strip('.fz'))
                newhdu = fits.HDUList()
                hdu = hdulist['SCI']
                hdu.header['OBSTYPE'] = 'OBJECT'
                hdu.header['OBJECT'] = target
                newhdu.append(hdu)
                newhdu.writeto(fullfilename.strip('.fz'), overwrite=True)

    def make_location(self, telescope):
        return({'telescope_class': telescope})

    def make_constraints(self):
        max_airmass = self.params['constraints']['max_airmass']
        min_lunar_distance = self.params['constraints']['min_lunar_distance']
        constraints = {
            'max_airmass': max_airmass,
            'min_lunar_distance': min_lunar_distance
        }
        return(constraints)

    # This sets up a target based on name, ra, dec
    def make_target(self, name, ra, dec):
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
        return(target)

    # This sets up the instrument parameters for an obs
    # Imaging specific
    def make_instrument_configs(self, filt, exptime, strat, readout='SLOW'):
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
        return([configuration])

    def make_acquisition_config(self, strat, mode=None):
        if not mode:
            mode = strat['acquisition_config']
        return({'mode': mode})

    def make_guiding_config(self, strat):
        mode = strat['guiding_config']
        optional = True
        if 'spec' in strat['type']:
            optional=False
        return({'mode': mode, 'optional': optional})

    """
    Window object for LCOGT API v3.
    start time = earliest the obs can be executed
    duration = amount of time in which request can be executed
    """
    def make_window(self, start, duration):
        fmt = self.params['date_format']
        end = start + timedelta(days=duration)
        window = [{
            'start': start.strftime(fmt),
            'end': end.strftime(fmt)
        }]
        return(window)

    def post_user_request(self, user_request):
        # Get username and password
        username, password = self.get_username_password()

        # Get authorization token
        header = self.get_token_header(username, password, auth_type='request')
        uri = self.params['uri']['request']['request']
        response = requests.post(uri, json=user_request, headers=header)
        return(response)

    def get_exposure_time(self, filt, mag, strat):
        # Get preferred snr given input magnitude
        # Handle 'spec' case
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

    # Get the maximum allowable IPP for a requestgroup
    def get_max_allowable_ipp(self, obs_requests):

        # Get username and password
        username, password = self.get_username_password()

        uri = self.params['uri']['request']['request']+'max_allowable_ipp/'
        header = self.get_token_header(username, password, auth_type='request')

        response = requests.post(uri, json=obs_requests, headers=header)

        return(response.json())

    """
    LCO API v3 configurations.
    configurations consist of constraints, instrument configuration, acquisition
    configuration, guiding configuration, target, instrument type, exposure
    type, and observation priority.  This method will construct a series of
    configurations for a single target.
    """
    def make_configurations(self, obj, ra, dec, mag, strat, inst_config=[]):
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


    """
    LCO API v3 request.
    requests is a list of requests each defined by a location, configuration,
    and a window.
    """
    def make_requests(self, obj, ra, dec, mag, strat, start=None,
        inst_config={}):

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

    """
    LCO API v3 observation request.
    Each method constructs one layer of the request.  Here we need to make the
    outermost obs request with name, proposal, ipp value, operator, observation
    type, and the individual requests.  Send the
    """
    def make_obs_request(self, obj, ra, dec, mag, strategy = 'default',
        propidx=0, start=None, recalculate_ipp=False, inst_config=[]):

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
            print('Recalculating IPP...')
            response = self.get_max_allowable_ipp(obs_request)
            print(response)
            try:
                if 'errors' in response.keys() and 'requests' in response['errors'].keys():
                    data = response['errors']['requests']
                    if 'non_field_errors' in data[0].keys():
                        print('ERROR in obs request',obs_request)
                        return(None)
                elif 'errors' in response.keys():
                    data = response['errors']
                    if 'non_field_errors' in data.keys():
                        print('ERROR in obs request',obs_request)
                        return(None)
            except:
                pass
            max_ipp = None
            for key in response.keys():
                r = response[key]
                for key in r.keys():
                    ipp_data = r[key]
                    print(ipp_data)
                    if 'max_allowable_ipp_value' in ipp_data.keys():
                        max_ipp = ipp_data['max_allowable_ipp_value']

            if max_ipp:
                curr_ipp = obs_request['ipp_value']
                if max_ipp < curr_ipp:
                    print('Adjusting IPP value to:',max_ipp)
                    obs_request['ipp_value']=max_ipp
                else:
                    m = 'No need to adjust IPP, {0} < {1}'.format(curr_ipp,
                        max_ipp)
                    print(m)


        if not requests[0]['configurations']:
            return(None)
        else:
            response = self.post_user_request(obs_request)
            response = response.json()
            return(response)
