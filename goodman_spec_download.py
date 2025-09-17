from lcogt import lcogt
from astropy import units as u
from astropy.io import fits
from astropy.time import Time
from astropy.time import TimeDelta
import glob
import sys
import os
import argparse
import copy

def parse_arguments(usage=''):

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('date', type=str, default='.',
                        help='Input date to download data.')
    parser.add_argument('propid', type=str, default='',
                        help='Proposal ID from which to download data '+\
                        '(can be comma-separated list).')
    parser.add_argument('--outdir', type=str, default='.',
                        help='Base output directory for downloading data.')

    args = parser.parse_args()

    return(args)

def main(date, propid, outdir='.'):

    lco = lcogt()

    args = parse_arguments()

    # Start time is noon in Chile on given date
    sdate = Time(date) + TimeDelta(14 * 3600 * u.s)
    # End time is noon in Chile on the next day
    edate = sdate + TimeDelta(24 * 3600 * u.s)

    basedir = outdir
    datedir = Time(date).datetime.strftime('%Y%m%d')
    fulloutdir = os.path.join(outdir, datedir, 'rawdata')
    fullworkdir = os.path.join(outdir, datedir, 'workspace')

    if not os.path.exists(fulloutdir):
        os.makedirs(fulloutdir)
    if not os.path.exists(fullworkdir):
        os.makedirs(fullworkdir)

    for propid in propid.split(','):

        telid='4m0a'

        obslist = lco.get_obslist(sdate=sdate, edate=edate,
                propid=[propid], telid=telid, rlevel=[25],
                obstype='SPECTRUM')

        # Dont know of any other way to distinguish raw from reduced files, so
        # use this for now until they update RLEVEL parameter
        outobslist = []
        for obs in obslist:
            filename = obs['filename']
            if filename.startswith('cfzst'): continue

            outobslist.append(obs)

        obslist = copy.copy(outobslist)

        # Get any standard star calibration files
        callist = lco.get_obslist(sdate=sdate, edate=edate,
            propid=['calibrate'], telid=telid, rlevel=[25],
            obstype='SPECTRUM')

        instlist = []
        for obs in obslist:
            if obs['INSTRUME'].lower() not in instlist:
                instlist.append(obs['INSTRUME'].lower())
        for obs in callist:
            if obs['INSTRUME'].lower() not in instlist:
                instlist.append(obs['INSTRUME'].lower())

        if len(instlist)==0:
            print('No observations for today.  Exiting...')
            sys.exit()

        # Get associated calibration images
        flatlist = lco.get_obslist(sdate=sdate, edate=edate, telid=telid, 
            propid=['calibrate'], rlevel=[25], obstype='LAMPFLAT')

        # Get associated arc lamp images
        arclist = []
        for pid in [propid,'calibrate']:
            arcs = lco.get_obslist(sdate=sdate, edate=edate, telid=telid, 
                propid=[pid],
                #propid=[propid,'calibrate'], 
                rlevel=[0], obstype='ARC')
            arclist.extend(arcs)

        flatlist = [f for f in flatlist if f['INSTRUME'].lower() in instlist]
        arclist = [f for f in arclist if f['INSTRUME'].lower() in instlist]
        arclist = [f for f in arclist if not f['basename'].lower().startswith('cfzst')]
        arclist = [f for f in arclist if not f['basename'].lower().startswith('wecfzst')]
        callist = [f for f in callist if not f['basename'].lower().startswith('cfzst')]
        callist = [f for f in callist if not f['basename'].lower().startswith('wecfzst')]

        nsci = len(obslist)
        ncal = len(callist)
        nfla = len(flatlist)
        narc = len(arclist)

        print(f'Need to download {nsci} science, {ncal} standards, {nfla} flats, {narc} arcs')

        lco.download_obslist(obslist, outrootdir=fulloutdir, skip_header=True, funpack=False)
        lco.download_obslist(callist, outrootdir=fulloutdir, skip_header=True, funpack=False)
        lco.download_obslist(flatlist, outrootdir=fulloutdir, skip_header=True, funpack=False)
        lco.download_obslist(arclist, outrootdir=fulloutdir, skip_header=True, funpack=False)

if __name__=='__main__':
    args = parse_arguments()
    main(args.date, args.propid, outdir=args.outdir)
