import copy
import os
import argparse

import yaml

import pandas as pd

from astropy.coordinates import Angle

from speedyfit import photometry_query
from speedyfit.speedyfit import fit_sed, get_observations, write_results, plot_results
from speedyfit.default_setup import default_double, default_single


def process_objects(data):

    def convert_ra(ra):
        if ' ' in str(ra).strip() or ':' in str(ra).strip():
            ra = str(ra).strip()
            ra = ra.replace(' ', '').replace(':', '')
        else:
            a = Angle(ra, unit='degree').hms
            ra = '{:02.0f}{:02.0f}{:05.2f}'.format(*a)
        return ra

    def convert_dec(dec):
        if ' ' in str(dec).strip() or ':' in str(dec).strip():
            dec = str(dec).strip()
            dec = dec.replace(' ', '').replace(':', '')
            if not '-' in dec and not '+' in dec:
                dec = '+' + dec
        else:
            a = Angle(dec, unit='degree').dms
            dec = '{:+03.0f}{:02.0f}{:05.2f}'.format(a[0], abs(a[1]), abs(a[2]))
        return dec

    if 'name' in data.columns.values:
        object_list = data['name'].apply(lambda x: x.replace(' ', '_')).values
    else:
        # deal with coordinates
        ra_ = data['ra'].apply(convert_ra)
        dec_ = data['dec'].apply(convert_dec)

        name = ['J{}{}'.format(r, d) for r, d in zip(ra_, dec_)]

        object_list = pd.Series(data=name)

    return object_list


def read_setup(setup):

    if os.path.isfile(setup):
        infile = open(setup, 'r')
        setup = infile.readlines()
        infile.close()
        setup = ''.join(setup)

    elif setup == 'single':
        setup = default_single
    elif setup == 'double':
        setup = default_double
    else:
        print('Setup not recognized! Using single as default')
        setup = default_single

    return setup


def prepare_setup(object_list, basedir, default_setup):

    for objectname in object_list:

        if not os.path.isdir(basedir+'/'+objectname):
            os.mkdir(basedir+'/'+objectname)

        plx, e_plx = photometry_query.get_parallax(objectname)

        out = copy.copy(default_setup)
        if '<photfilename>' in out:
            photfilename_ = basedir + '/' + objectname + '/' + objectname + '.phot'
            out = out.replace('<photfilename>', photfilename_)
        if '<objectname>' in out:
            objectname_ = basedir + '/' + objectname + '/' + objectname
            out = out.replace('<objectname>', objectname_)
        if '<plx>' in out:
            out = out.replace('<plx>', str(plx))
            out = out.replace('<e_plx>', str(e_plx))

        filename = basedir + '/' + objectname + '/' + objectname + '_setup.yaml'

        ofile = open(filename, 'w')
        ofile.write(out)
        ofile.close()


def prepare_photometry(object_list, basedir, skip_existing=True):

    for objectname in object_list:

        if os.path.isfile(basedir+'/'+objectname+'/'+objectname+'.phot') and skip_existing:
            continue

        if not os.path.isdir(basedir+'/'+objectname):
            os.mkdir(basedir+'/'+objectname)

        try:
            filename = basedir + '/' + objectname + '/' + objectname + '.phot'
            photometry = photometry_query.get_photometry(objectname, filename=filename)
        except Exception as e:
            print("Failed to obtain photometry for {}".format(objectname))
            print(e)


def fit_seds(object_list, basedir):

    all_results = {}

    for objectname in object_list:

        # read the setup
        filename = basedir + '/' + objectname + '/' + objectname + '_setup.yaml'
        try:
            infile = open(filename)
            setup = yaml.safe_load(infile)
            infile.close()
        except Exception as e:
            print("Could not read setupfile for {}".format(objectname))
            print(e)
            continue

        try:
            # get the photometry
            photbands, obs, obs_err = get_observations(setup)

            # do the SED fit
            results, samples, constraints, gridnames = fit_sed(setup, photbands, obs, obs_err)

            # add results to results dic
            all_results.setdefault('objectname', []).append(objectname)
            for k, v in results.items():
                all_results.setdefault(k, []).append(v[0])
                all_results.setdefault(k+'_err', []).append(( v[2] + v[3] ) / 2.0)

            # write the results
            write_results(setup, results, samples, obs, obs_err, photbands)

            # make any requested plots
            plot_results(setup, results, samples, constraints, gridnames, obs, obs_err, photbands)

        except Exception as e:
            print("Could not fit {}".format(objectname))
            print(e)

    all_results = pd.DataFrame(data=all_results)
    all_results.to_csv('results.csv')
    print(all_results)



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

    object_data = pd.read_csv(args.filename)
    object_list = process_objects(object_data)

    # create the setup file for all systems
    prepare_setup(object_list, basedir, default_setup)

    # try to obtain photometry for all systems if requested
    if args.phot:
        prepare_photometry(object_list, basedir)

    # fit all systems if requested
    if args.fit:
        fit_seds(object_list, basedir)