from astropy import units as u
from astropy.io import fits
from astropy.time import Time
from astropy.time import TimeDelta
import glob
import sys
import os
import argparse
import copy
import shutil

def parse_arguments(usage=''):

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('date', type=str, default='.',
                        help='Input date to download data.')
    parser.add_argument('input_dir', type=str, default='.',
                        help='Path to the directory of input data.')
    parser.add_argument('--outdir', type=str, default='.',
                        help='Base output directory for downloading data.')

    args = parser.parse_args()

    return(args)

def main(date, input_dir, outdir='.'):

    args = parse_arguments()

    basedir = outdir    
    datedir = Time(date).datetime.strftime('%Y%m%d')
    fulloutdir = os.path.join(outdir, datedir, 'rawdata')
    fullworkdir = os.path.join(outdir, datedir, 'workspace')

    if not os.path.exists(fulloutdir):
        os.makedirs(fulloutdir)
    if not os.path.exists(fullworkdir):
        os.makedirs(fullworkdir)

    input_files = glob.glob(os.path.join(input_dir, '*.fits'))
    input_files += glob.glob(os.path.join(input_dir, '*.fz'))

    for file in input_files:

        hdu = fits.open(file)

        # Edit headers
        obstype = hdu[0].header['OBSTYPE']
        if obstype.lower()=='object':
            hdu[0].header['OBSTYPE']='SPECTRUM'

        outfile = os.path.join(fulloutdir, os.path.basename(file))

        hdu.writeto(outfile, overwrite=True, output_verify='silentfix')

        if outfile.endswith('.fits'):
            packfile = outfile+'.fz'

            if not os.path.exists(packfile):
                cmd = f'fpack {outfile}'
                print(cmd)
                os.system(cmd)

            os.remove(outfile)


if __name__=='__main__':
    args = parse_arguments()
    main(args.date, args.input_dir, outdir=args.outdir)
