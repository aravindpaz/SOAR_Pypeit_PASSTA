from lcogt import lcogt
import astropy
import urllib
from astropy import units as u
from astropy.io import fits
from astropy.io import ascii
from astropy.time import Time
from astropy.table import Table
from astropy.table import vstack
from astropy.table import unique
import numpy as np
import glob
import sys
import os
import argparse
import copy
import time
import shutil

def refresh_astropy_cache():

    tries = 0
    while tries<5:
        try:
            astropy.coordinates.EarthLocation.get_site_names(refresh_cache=True)
            time.sleep(0.5)
            astropy.coordinates.EarthLocation.get_site_names()
            break
        except urllib.error.URLError:
            tries += 1
            time.sleep(5)

refresh_astropy_cache()
from pypeit import inputfiles

def parse_arguments(usage=''):

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('date', type=str, default='.',
                        help='Input date to download data.')
    parser.add_argument('outdir', type=str, default='.',
                        help='Base directory for pypeit reductions.')

    args = parser.parse_args()

    return(args)


def parse_pypeit_output(scidir, pypeit_table, caldb, slit_pos=505, pix_tolerance=[15, 30, 60]):
    """
    Takes the output from a science directory and parses the output to produce
    list of 1D files or standards that need to be converted to sensitivity functions.
    """

    objdata = []
    for row in pypeit_table.data:
        if row['frametype']=='science':

            filename = row['filename']
            target = row['target'].lower()
            ra = row['ra']
            dec = row['dec']
            mjd = row['mjd']
            airmass = row['airmass']
            dispname = row['mode']
            slitname = row['decker']

            framebase = filename.replace('.fits','').replace('.fz','')
            outfiles = glob.glob(os.path.join(scidir, f'spec1d_{framebase}*.txt'))

            if len(outfiles)==0:
                print(f'WARNING: no traces detected for {framebase}')
                continue
            else:
                outfile = outfiles[0]

            spectable = ascii.read(outfile, delimiter='|', header_start=0)

            # Get the objects that are within tolerance of expected slit position
            # Do this iteratively for different values of pix_tolerance
            for pix_tol in pix_tolerance:
                mask = np.abs(spectable['spat_pixpos']-slit_pos) < pix_tol
                subtable = spectable[mask]

                print(subtable,f'for slit_pos={slit_pos} and pix_tol={pix_tol}')

                if len(subtable)==0:
                    continue
                else:
                    break

            if len(subtable)==0:
                print(f'WARNING: could not find a good trace in: {filename}')
                continue

            # Take the highest signal source at this position
            idx = np.argmax(subtable['s2n'].data)

            objname = subtable[idx]['name']
            objdata.append((outfile, objname, target, mjd, dispname, slitname))

        elif row['frametype']=='standard':

            framebase = row['filename'].replace('.fits','').replace('.fz','')
            outfile = glob.glob(os.path.join(scidir, f'spec1d_{framebase}*.fits'))[0]

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
            senslog = os.path.join(scidir, sensfilename).replace('.fits','.log')
            par_outfile = os.path.join(scidir, 'sensfunc.par')

            if not os.path.exists(caldb):
                os.makedirs(caldb)
            if not os.path.exists(os.path.join(caldb, 'plots')):
                os.makedirs(os.path.join(caldb, 'plots'))

            cmd = f'pypeit_sensfunc {outfile} -o {fullsensfilename} --par_outfile {par_outfile} > {senslog} 2> {senslog}'
            print(cmd)
            os.system(cmd)

            for file in glob.glob(os.path.join(caldb, '*.pdf')):
                newfile = os.path.join(caldb, 'plots', os.path.basename(file))
                if not os.path.exists(newfile):
                    shutil.move(file, os.path.join(caldb, 'plots'))
                else:
                    os.remove(file)

            if os.path.exists(fullsensfilename):

                caldbfile = os.path.join(caldb, 'caldb.txt')
                if not os.path.exists(caldbfile):
                    names = ('filename','target','ra','dec','mjd','dispname','slitname','airmass')
                    dtype = (str, str, np.float64, np.float64, np.float64, str, str, np.float64)
                    caldb_table = Table(names=names, dtype=dtype)

                    caldb_table.add_row((os.path.basename(sensfilename), target, ra, dec, mjd, dispname, slitname, airmass))
                    caldb_table = unique(caldb_table)

                    caldb_table.write(caldbfile, format='ascii.ecsv', overwrite=True)
                else:
                    caldb_table = ascii.read(caldbfile)
                    caldb_table.add_row((os.path.basename(sensfilename), target, ra, dec, mjd, dispname, slitname, airmass))
                    caldb_table = unique(caldb_table)

                    caldb_table.write(caldbfile, format='ascii.ecsv', overwrite=True)

    return(objdata)


def make_flux_file(objdata, fullsensfile):

    outfile = objdata[0][0]
    target = objdata[0][2]
    dirname, _ = os.path.split(outfile)

    flux_filename = os.path.join(dirname, f'{target}.flux')

    with open(flux_filename, 'w') as f:

        f.write('[fluxcalib] \n')
        f.write('  extinct_correct = False \n\n')
        f.write('flux read \n')
        f.write(f'  path {dirname} \n')
        f.write('filename | sensfile \n')

        for obj in objdata:
            filename = obj[0]
            if filename.endswith('.txt'):
                filename = filename.replace('.txt','.fits')
            f.write(f'{filename} | {fullsensfile} \n')

        f.write('flux end \n')

    if os.path.exists(flux_filename):
        return(flux_filename)

def make_coadd_file(objdata, fullout_1dfilename):

    outfile = objdata[0][0]
    target = objdata[0][2]
    dirname, _ = os.path.split(outfile)

    coadd_filename = os.path.join(dirname, f'{target}.coadd1d')

    with open(coadd_filename, 'w') as f:

        f.write('[coadd1d] \n')
        f.write(f'  coaddfile = {fullout_1dfilename} \n')
        f.write(f'  wave_method = linear \n\n')
        f.write('coadd1d read \n')
        f.write(f'  path {dirname} \n')
        f.write('filename | obj_id \n')

        for obj in objdata:
            filename = obj[0]
            if filename.endswith('.txt'):
                filename = filename.replace('.txt','.fits')
            objname = obj[1]
            f.write(f'{filename} | {objname} \n')

        f.write('coadd1d end \n')

    if os.path.exists(coadd_filename):
        return(coadd_filename)



def pypeit_post_process(objdata, specdir, caldb_dir):

    caldb_file = os.path.join(caldb_dir, 'caldb.txt')
    if not os.path.exists(caldb_file):
        raise Exception(f'Could not locate calibration database: {caldb_file}')

    caldb = ascii.read(caldb_file)

    outlinks = []

    firstobj = objdata[0]
    
    outfile, objname, target, mjd, dispname, slitname = firstobj
    dirname, _ = os.path.split(outfile)

    datestr = Time(mjd, format='mjd').datetime.strftime('%Y%m%d')

    mask = (caldb['dispname']==dispname) & (caldb['slitname']==slitname)
    goodcals = caldb[mask]

    # Get closest flux standard in time from good calibrations 
    idx = np.argmin(np.abs(goodcals['mjd']-mjd))

    fullsensfile = os.path.join(caldb_dir, goodcals[idx]['filename'])

    flux_filename = make_flux_file(objdata, fullsensfile)
    flux_logname = os.path.join(dirname, 'flux.log')
    par_outfile = os.path.join(dirname, 'fluxing.par')

    # Run pypeit flux calibration on flux file
    cmd = f'pypeit_flux_calib --par_outfile {flux_filename} > {flux_logname} 2> {flux_logname}'
    print(cmd)
    os.system(cmd)

    # Now that files are fluxed, need to coadd into one file
    out_1dfilename = f'{target}.{datestr}.{dispname}.coadd1d.fits'
    fullout_1dfilename = os.path.join(dirname, out_1dfilename)
    full_tellcorrfilename = fullout_1dfilename.replace('.fits','_tellcorr.fits')

    coadd_filename = make_coadd_file(objdata, fullout_1dfilename)
    coadd_logname = os.path.join(dirname, 'coadd.log')
    par_outfile = os.path.join(dirname, 'coadd1d.par')

    cmd = f'pypeit_coadd_1dspec {coadd_filename} --par_outfile {par_outfile} > {coadd_logname} 2> {coadd_logname}'
    print(cmd)
    os.system(cmd)

    if os.path.exists(fullout_1dfilename):
        tellcorr_logname = os.path.join(dirname, 'tellcorr.log')
        par_outfile = os.path.join(dirname, 'tellfit.par')

        # Run telluric correction
        os.chdir(dirname)
        cmd = f'pypeit_tellfit {fullout_1dfilename} --objmodel poly --par_outfile {par_outfile} > {tellcorr_logname} 2> {tellcorr_logname}'
        print(cmd)
        os.system(cmd)

    if os.path.exists(full_tellcorrfilename):
        outlink = os.path.join(specdir, os.path.basename(full_tellcorrfilename))

        if not os.path.exists(outlink):
            os.symlink(full_tellcorrfilename, outlink)

        if os.path.exists(outlink):
            print(f'Linked {full_tellcorrfilename}->{outlink}')
            return(outlink)

def main(date, outdir, caldb_dir):
    """
    Main driver script for the pypeit reduction pipeline.  Takes downloaded files
    from goodman_spec_download.py and processes them with pypeit to produce 
    calibrated 1D spectra.
    """

    lco = lcogt()

    basedir = outdir
    datedir = Time(date).datetime.strftime('%Y%m%d')
    fulloutdir = os.path.join(outdir, datedir, 'rawdata')
    fullworkdir = os.path.join(outdir, datedir, 'workspace')
    fullspecdir = os.path.join(outdir, 'spectra')

    # These should already exist from goodman_spec_download.py, if not then exit
    if not os.path.exists(fulloutdir) and not os.path.exists(fullworkdir):
        print(f'No data to process for {date}...')
        sys.exit()
    
    # Make output spectrum directory if it does not exist
    if not os.path.exists(fullspecdir):
        os.makedirs(fullspecdir)

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

        tries = 0
        while tries<5:
            try:
                refresh_astropy_cache()
                cmd = f'pypeit_setup -s {spec} -c all -r {objdir} -d {objdir} > {objdir}/setup.log 2> {objdir}/setup.log'
                print(cmd)
                os.system(cmd)
                break
            except astropy.coordinates.errors.UnknownSiteException:
                tries+=1

        if tries==5:
            raise Exception('Could not get astropy cache to refresh...')

    all_objdata = []
    for objdir in reduction_dirs:
        pypeit_file = os.path.join(objdir, f'{spec}_A/{spec}_A.pypeit')
        pypeit_table = inputfiles.PypeItFile.from_file(pypeit_file)
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

        objdata = parse_pypeit_output(science_dir, pypeit_table, caldb_dir)
        if len(objdata)>0:
            all_objdata.append(objdata)

    outlinks = []
    for objdata in all_objdata:
        scidir, _ = os.path.split(objdata[0][0])
        outlinks.append(pypeit_post_process(objdata, fullspecdir, caldb_dir))


if __name__=="__main__":
    args = parse_arguments()
    file_dir, _ = os.path.split(__file__)
    caldb_dir = os.path.join(file_dir, 'caldb')
    main(args.date, args.outdir, caldb_dir)
