
import os
import argparse

import pandas as pd

from astropy.coordinates import Angle

from speedyfit import photometry_query
from speedyfit.default_setup import default_double, default_single


def process_objects(filename):

    def convert_ra(ra):
        if ' ' in str(ra).strip() or ':' in str(ra).strip():
            ra = str(ra).strip()
            ra.replace(' ', ':')
        else:
            a = Angle(ra, unit='degree').hms
            ra = '{:.0f}:{:.0f}:{:.2f}'.format(*a)
        return ra

    def convert_dec(dec):
        if ' ' in str(dec).strip() or ':' in str(dec).strip():
            dec = str(dec).strip()
            dec.replace(' ', ':')
        else:
            a = Angle(ra, unit='degree').dms
            dec = '{:.0f}:{:.0f}:{:.2f}'.format(*a)
            pass

        return dec

    data = pd.read(filename)

    if 'name' in data.columns.values:
        object_list = data['name'].apply(lambda x: x.replace(' ', '_')).values
    else:
        # deal with coordinates
        ra = data['ra'].apply(convert_ra)
        dec = data['dec'].apply(convert_dec)



        pass

    return object_list


def read_setup(setup):

    if os.path.isfile(setup):
        infile = os.open(setup, 'r')
        setup = infile.read_lines()
        infile.close()
        setup = '/n'.join(setup)

    elif setup == 'single':
        setup = default_single
    elif setup == 'double':
        setup = default_double
    else:
        print('Setup not recognized! Using single as default')
        setup = default_single

    return setup


def prepare_setup(object_list, basedir, default_setup):

    for i, objectname in object_list.iter_lines():

        if not os.path.isdir(basedir+'/'+objectname):
            os.mkdir(basedir+'/'+objectname)

        plx, e_plx = photometry_query.get_parallax(objectname)

        out = default_setup.copy()
        if '<photfilename>' in out:
            out = out.replace('<photfilename>', objectname + '.phot')
        if '<objectname>' in out:
            out = out.replace('<objectname>', objectname)
        if '<plx>' in out:
            out = out.replace('<plx>', str(plx))
            out = out.replace('<e_plx>', str(e_plx))

        filename = basedir + '/' + objectname + '/' + objectname + '_setup.yaml'

        ofile = open(filename, 'w')
        ofile.write(out)
        ofile.close()


def prepare_photometry(object_list, basedir, skip_existing=True):

    for i, objectname in object_list.iter_lines():

        if os.path.isfile(basedir+'/'+objectname+'/'+objectname+'.phot') and skip_existing:
            continue

        if not os.path.isdir(basedir+'/'+objectname):
            os.mkdir(basedir+'/'+objectname)

        try:
            photometry = photometry_query.get_photometry(objectname, filename=objectname + '.phot')
        except Exception as e:
            print("Failed to obtain photometry for {}".format(objectname))
            print(e)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', action="store", type=str,
                        help='file containing a list with the names or coordinates of the systems to fit')
    parser.add_argument("-setup", dest='setup', type=str, default='single',
                        help="The path to the setup file to use, or a standard: 'single' or 'double'. "
                             "(default = 'single'")
    parser.add_argument("-basedir", dest='basedir', type=str, default='./',
                        help="The directory where all photometry, setup and results will be stored. (default = ./)")
    parser.add_argument("--fit", dest='fit', action='store_true',
                        help="Fit the systems")
    parser.add_argument('--phot', dest='phot', action='store_true',
                        help='Obtain photometry of the systems')
    parser.add_argument('--noplot', dest='noplot', action='store_true',
                        help="Don't show any plots, only store to disk.")
    args, variables = parser.parse_known_args()

    basedir = args.basedir
    default_setup = read_setup(args.setup)

    if not os.path.isdir(basedir):
        os.mkdir(basedir)

    object_list = process_objects(args.filename)

    # create the setup file for all systems
    prepare_setup(object_list, basedir, default_setup)

    # try to obtain photometry for all systems if requested
    if args.phot:
        prepare_photometry(object_list, basedir)

    # fit all systems if requested
    if args.fit:
        pass