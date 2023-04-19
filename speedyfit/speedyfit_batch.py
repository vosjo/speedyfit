import copy
import os
import argparse

import yaml

import pandas as pd
import numpy as np

from speedyfit import photometry_query
from speedyfit.main import fit_sed, get_observations, write_results, plot_results
from speedyfit.default_setup import default_tmap, default_binary
from speedyfit.utils import get_name_from_coordinates

def extract_constraint_names(data):
    columns = data.columns.values
    constraint_names = []
    for name in columns:
        if 'constraint_' in name:
            cname = name.split('_')[1]
            constraint_names.append(cname)

    return constraint_names


def process_objects(data):

    # check if we are working with names or if we need to make names from the coordinates
    if 'name' in data.columns.values:
        object_list = data['name'].apply(lambda x: x.replace(' ', '_').strip()).values
        object_list = pd.DataFrame({'name': object_list})
    else:
        # convert coordinates into a name if no name was given.
        name = get_name_from_coordinates(data['ra'].values, data['dec'].values)
        object_list = pd.DataFrame({'name': name})

    if 'ra' in data and 'dec' in data:
        object_list['ra'] = data['ra']
        object_list['dec'] = data['dec']

    # now check for constraints
    constraint_cols = extract_constraint_names(data)
    for cname in constraint_cols:
        object_list[cname] = data['constraint_'+cname]

    return object_list


def read_setup(setup):

    if os.path.isfile(setup):
        infile = open(setup, 'r')
        setup = infile.readlines()
        infile.close()
        setup = ''.join(setup)

    elif setup == 'single':
        setup = default_tmap
    elif setup == 'binary':
        setup = default_binary
    else:
        print('Setup not recognized! Using single as default')
        setup = default_tmap

    return setup


def prepare_setup(object_list, basedir, default_setup):

    constraint_names = list(object_list.columns.values)
    constraint_names.remove('name')
    constraint_names.remove('ra')
    constraint_names.remove('dec')

    for i, row in object_list.iterrows():
        objectname = row['name']
        ra, dec = row['ra'], row['dec']
        print(f'Preparing setup files for object: {objectname} ({ra}  {dec})')

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
        if '<constraints>' in out:
            constraints = f"\n  parallax: [{plx:0.4f}, {e_plx:0.4f}]"
            for cname in constraint_names:
                val = row[cname]
                try:
                    val = val.replace('(', '').replace(')', '')
                except:
                    continue
                if val.strip() == '':
                    continue
                constraints += f"\n  {cname}: [{val}]"
            out = out.replace('<constraints>', constraints)

        filename = basedir + '/' + objectname + '/' + objectname + '_setup.yaml'

        ofile = open(filename, 'w')
        ofile.write(out)
        ofile.close()


def prepare_photometry(object_list, basedir, skip_existing=True):

    for i, row  in object_list.iterrows():
        objectname = row['name']

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


def fit_seds(object_list, basedir, start_index=0, stop_index=None):

    all_results = {}

    stop_index = len(object_list) if stop_index is None else stop_index
    object_list_sel = object_list[start_index:stop_index+1]

    for i, row in object_list_sel.iterrows():
        objectname = row['name']
        print(f"Starting fit {i}/{len(object_list)} for object: {objectname}")

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
                        help="The path to the setup file to use, or a standard: 'single' or 'binary'. "
                             "(default = 'single'")
    parser.add_argument("-basedir", dest='basedir', type=str, default='./',
                        help="The directory where all photometry, setup and results will be stored. (default = ./)")
    parser.add_argument("--fit", dest='fit', action='store_true',
                        help="Fit the systems")
    parser.add_argument('--phot', dest='phot', action='store_true',
                        help='Obtain photometry of the systems')
    parser.add_argument('-start', dest='start_index', type=int, default=0,
                        help='Index of the object with which to start the fit, defaults to 0')
    parser.add_argument('-stop', dest='stop_index', type=int, default=None,
                        help='Index of the object with which to end the fitting, defaults to None (all objects are fitted).')
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
        fit_seds(object_list, basedir, start_index=args.start_index, stop_index=args.stop_index)