from lcogt import lcogt
import astropy
import urllib
from astropy import units as u
from astropy.io import fits
from astropy.time import Time
import numpy as np
import glob
import sys
import os
import argparse
import copy
import time

lco = lcogt('/home/ckilpatrick/scripts/tokens/shibboleth')

def parse_arguments(usage=''):

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('date', type=str, default='.',
                        help='Input date to download data.')
    parser.add_argument('outdir', type=str, default='.',
                        help='Base directory for pypeit reductions.')

    args = parser.parse_args()

    return(args)

def refresh_astropy_cache():

    tries = 0
    while tries<5:
        try:
            astropy.coordinates.EarthLocation.get_site_names(refresh_cache=True)
            break
        except urllib.error.URLError:
            tries += 1
            time.sleep(5)



def main(date, outdir):

    basedir = outdir
    datedir = Time(date).datetime.strftime('%Y%m%d')
    fulloutdir = os.path.join(outdir, datedir, 'rawdata')
    fullworkdir = os.path.join(outdir, datedir, 'workspace')

    if not os.path.exists(fullworkdir):
        os.makedirs(fullworkdir)

    # Create one reduction directory per object so they have the correct arc
    files = sorted(glob.glob(os.path.join(fulloutdir, '*.fz')))
    all_data = {}

    spec = 'soar_goodman_red'

    for file in files:
        hdu = fits.open(file, mode='update')

        wavemode = hdu[1].header['WAVMODE']
        obstype = hdu[1].header['OBSTYPE']
        objname = lco.sanitize_objname(hdu[1].header['OBJECT'])
        hdu[1].header['OBJECT']=objname
        mjd = Time(hdu[1].header['DATE-OBS']).mjd

        if hdu[1].header['INSTRUME']=='GHTS_BLUE':
            spec = 'soar_goodman_blue'

        data = {'wavemod': wavemode, 'obstype': obstype, 'filename': file,
            'object': objname, 'mjd': mjd}

        # pypeit only handles these two setups for now
        if wavemode not in ['400_M1','400_M2']:
            continue

        if wavemode not in all_data.keys():
            all_data[wavemode]=[data]
        else:
            all_data[wavemode].append(data)

        hdu.close()


    reduction_dirs = []

    for wavemode in all_data.keys():
        wavedir = os.path.join(fullworkdir, wavemode)
        if not os.path.exists(wavedir):
            os.makedirs(wavedir)

        wavedata = all_data[wavemode]
        wavedata = sorted(wavedata, key=lambda x: x['mjd'])
        arcdata = [d for d in wavedata if d['obstype']=='ARC']
        all_objs = np.unique([d['object'] for d in wavedata if d['obstype']=='SPECTRUM'])


        for objname in all_objs:
            objdir = os.path.join(wavedir, objname)
            if not os.path.exists(objdir):
                os.makedirs(objdir)

            mjds = []
            for data in wavedata:
                filename = data['filename']
                basefile = os.path.basename(filename)
                outfile = os.path.join(objdir, basefile)

                if data['obstype']=='SPECTRUM':
                    mjds.append(data['mjd'])
                    if data['object']!=objname: continue
                # Will identify the best arc for this reduction after
                if data['obstype']=='ARC': continue

                if not os.path.exists(outfile):
                    print(f'Linking {filename}->{outfile}')
                    os.symlink(filename, outfile)

            # Identify a single arc closest in time to associate with these files
            # It is best for pypeit to use one arc then perform flexure correction
            # to the sky lines
            avg_mjd = np.average(mjds)
            arcdata = sorted(arcdata, key=lambda x: np.abs(x['mjd']-avg_mjd))
            filename = arcdata[0]['filename']
            basefile = os.path.basename(filename)
            outfile = os.path.join(objdir, basefile)

            if not os.path.exists(outfile):
                print(f'Linking {filename}->{outfile}')
                os.symlink(filename, outfile)

            reduction_dirs.append(objdir)

    
    for objdir in reduction_dirs:

        refresh_astropy_cache()
        cmd = f'pypeit_setup -s {spec} -c all -r {objdir} -d {objdir} > {objdir}/setup.log 2> {objdir}/setup.log'
        print(cmd)
        os.system(cmd)

    for objdir in reduction_dirs:
        pypeit_file = os.path.join(objdir, f'{spec}_A/{spec}_A.pypeit')
        science_dir = os.path.join(objdir, 'Science')

        if os.path.exists(pypeit_file):
            if not os.path.exists(science_dir):
                refresh_astropy_cache()
                cmd = f'run_pypeit {pypeit_file} -r {objdir} > {objdir}/reduction.log 2> {objdir}/reduction.log'
                print(cmd)
                os.system(cmd)
            else:
                print(f'Reduction already done for {objdir}')
        else:
            print(f'ERROR: pypeit_setup failed for {objdir}')


if __name__=="__main__":
    args = parse_arguments()
    main(args.date, args.outdir)
