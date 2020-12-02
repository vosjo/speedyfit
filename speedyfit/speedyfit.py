import sys
import yaml
import argparse
import numpy as np
import pylab as pl
import corner

from astropy.io import ascii

from numpy.lib.recfunctions import append_fields, repack_fields

from speedyfit import mcmc, model, plotting, fileio, filters, photometry_query
from speedyfit.default_setup import default_double, default_single



def select_photometry(photbands, obs, obs_err, remove_nan=True, remove_color=True, include=None, exclude=None,
                      verbose=True):
    """
    Function to select the wanted photometry.

    - removes photometry with nan values in value or error
    - removed colors
    - selects photometry based on include and exclude arrays
    """

    # -- remove photometry with nan values in measurement or error.
    if remove_nan:
        nani = np.isnan(obs) | np.isnan(obs_err)
        if any(nani) and verbose:
            print("Warning: there are NaN values in the following photometric bands:")
            for p in photbands[nani]:
                print("\t {}".format(p))
            print("They were removed")
        obs, obs_err, photbands = obs[~nani], obs_err[~nani], photbands[~nani]

    # -- remove colors
    if remove_color:
        color = np.array([filters.is_color(p) for p in photbands])
        s = np.where(~color)
        if any(color) and verbose:
            print("Warning: The following colors were removed:")
            for p in photbands[s]:
                print("\t {}".format(p))
        photbands, obs, obs_err = photbands[s], obs[s], obs_err[s]

    # -- only include bands that are requested based on the include keyword, or exclude all non wanted photometry
    #    based on the exclude keyword.
    if include is not None:
        photsys_include = [x for x in include if not '.' in x]
        photband_include = [x for x in include if '.' in x]
        incband = []
        for i, photband in enumerate(photbands):
            if photband in photband_include or photband.split('.')[0] in photsys_include:
                incband.append(i)
        incband = (np.array(incband),)
        photbands, obs, obs_err = photbands[incband], obs[incband], obs_err[incband]

    elif exclude is not None:
        photsys_exclude = [x for x in exclude if not '.' in x]
        photband_exclude = [x for x in exclude if '.' in x]
        incband = []
        for i, photband in enumerate(photbands):
            if photband not in photband_exclude and photband.split('.')[0] not in photsys_exclude:
                incband.append(i)
        incband = (np.array(incband),)
        photbands, obs, obs_err = photbands[incband], obs[incband], obs_err[incband]

    return photbands, obs, obs_err


def get_observations(setup):
    # -- parse photometry
    if isinstance(setup['photband_index'], str):
        data = ascii.read(setup['photometryfile'], format='fixed_width')
        photbands = np.array(data['band'])
        obs = np.array(data['flux'])
        obs_err = np.array(data['eflux'])
    else:
        data = ascii.read(setup['photometryfile'], data_start=0, header_start=None)

        setup['photband_index'] = data.colnames[setup['photband_index']]
        setup['obs_index'] = data.colnames[setup['obs_index']]
        setup['err_index'] = data.colnames[setup['err_index']]

        photbands = np.array(data[setup['photband_index']])
        obs = np.array(data[setup['obs_index']])
        obs_err = np.array(data[setup['err_index']])

    # -- select the requested photometry
    photbands, obs, obs_err = select_photometry(photbands, obs, obs_err, remove_nan=True, remove_color=True,
                                                include=setup.get('photband_include', None),
                                                exclude=setup.get('photband_exclude', None))

    return photbands, obs, obs_err


def fit_sed(setup, photbands, obs, obs_err):

    # -- pars limits
    pnames = setup['pnames']
    limits = np.array(setup['limits'])

    # -- pars constraints
    constraints = setup['constraints']
    for con, val in list(constraints.items()):
        if len(val) == 2:
            constraints[con] = [val[0], val[1], val[1]]

    if 'parallax' in constraints:
        p, pm, pp = constraints.pop('parallax')
        constraints['distance'] = [1000. / p, 1000. * pm / p ** 2, 1000. * pp / p ** 2]

    print("Applied constraints: ")
    for con, val in list(constraints.items()):
        print("\t {} = {} - {} + {}".format(con, val[0], val[1], val[2]))

    if 'distance' in constraints:
        # convert pc to Rsol
        constraints['distance'] = [44365810.04823812 * constraints['distance'][0],
                                   44365810.04823812 * constraints['distance'][1],
                                   44365810.04823812 * constraints['distance'][2], ]

    # -- pars limits on derived properties
    derived_limits = setup['derived_limits']

    # -- pars grid
    gridnames = setup['grids']
    grids = model.load_grids(gridnames, pnames, limits, photbands)

    # -- switch logg to g for a binary system
    if 'q' in constraints:
        if 'logg' in pnames:
            limits[pnames.index('logg')] = 10 ** limits[pnames.index('logg')]
            pnames[pnames.index('logg')] = 'g'

        if 'logg2' in pnames:
            limits[pnames.index('logg2')] = 10 ** limits[pnames.index('logg2')]
            pnames[pnames.index('logg2')] = 'g2'

        if 'logg' in constraints:
            g = 10 ** constraints['logg'][0]
            g_el = g * constraints['logg'][1] * np.log(10)
            g_eu = g * constraints['logg'][2] * np.log(10)
            constraints['g'] = [g, g_el, g_eu]

        if 'logg2' in constraints:
            g = 10 ** constraints['logg2'][0]
            g_el = g * constraints['logg2'][1] * np.log(10)
            g_eu = g * constraints['logg2'][2] * np.log(10)
            constraints['g2'] = [g, g_el, g_eu]

    # -- check for variables that are kept fixed
    fixed = np.where(limits[:, 0] == limits[:, 1])
    varied = np.where(limits[:, 0] != limits[:, 1])

    pnames = np.array(pnames)
    fixed_variables = {}
    for par, val in zip(pnames[fixed], limits[:, 0][fixed]):
        fixed_variables[par] = val

    pnames = list(pnames[varied])
    limits = limits[varied]

    # -- pars mcmc setup
    nwalkers = setup.get('nwalkers', 100)
    nsteps = setup.get('nsteps', 2000)
    nrelax = setup.get('nrelax', 500)
    a = setup.get('a', 10)

    # -- MCMC
    results, samples = mcmc.MCMC(obs, obs_err, photbands,
                                 pnames, limits, grids,
                                 fixed_variables=fixed_variables,
                                 constraints=constraints, derived_limits=derived_limits,
                                 nwalkers=nwalkers, nsteps=nsteps, nrelax=nrelax,
                                 a=a)

    # -- add fixed variables to results dictionary
    for par, val in list(fixed_variables.items()):
        results[par] = [val, val, 0, 0]

    # -- deal with the switch back to logg
    if 'g' in samples.dtype.names:
        samples = append_fields(samples, 'logg', np.log10(samples['g']), usemask=False)

    if 'g' in results:
        results['logg'] = np.log10(results['g'])
        results.pop('g')

    if 'g2' in samples.dtype.names:
        samples = append_fields(samples, 'logg2', np.log10(samples['g2']), usemask=False)

    if 'g2' in results:
        results['logg2'] = np.log10(results['g2'])
        results.pop('g2')

    _ = constraints.pop('g', None)
    _ = constraints.pop('g2', None)

    names = list(samples.dtype.names)
    if 'g' in names: names.remove('g')
    if 'g2' in names: names.remove('g2')
    samples = samples[names]

    percentiles = setup.get('percentiles', [16, 50, 84])
    pc = np.percentile(samples.view(np.float64).reshape(samples.shape + (-1,)), percentiles, axis=0)
    pars = {}
    for p, v, e1, e2 in zip(samples.dtype.names, pc[1], pc[1] - pc[0], pc[2] - pc[1]):
        results[p] = [results[p], v, e1, e2]
        pars[p] = v

    return results, samples, constraints, gridnames


def write_results(setup, results, samples, obs, obs_err, photbands):

    outpars, outvals = [], []
    for par in ['teff', 'logg', 'L', 'rad']:
        outpars.append(par)
        outpars.append(par + '_err')
        outvals.append(results[par][1])
        outvals.append(np.average([results[par][2], results[par][3]]))

    if 'teff2' in results:
        for par in ['teff2', 'logg2', 'L2', 'rad2']:
            outpars.append(par)
            outpars.append(par + '_err')
            outvals.append(results[par][1])
            outvals.append(np.average([results[par][2], results[par][3]]))

    resultfile = setup.get('resultfile', None)
    if resultfile is not None:
        import pandas as pd
        data = pd.DataFrame(data=[outvals], columns=outpars)
        data.to_csv(resultfile, index=False)

    datafile = setup.get('datafile', None)
    if not datafile is None:
        # -- get plain text settings to retain comments and ordering
        setupfile = open(args.filename)
        setup_str = "".join(setupfile.readlines())
        setupfile.close()

        fileio.write2fits(samples, datafile, setup=setup_str)

    h5file = setup.get('h5file', None)
    if h5file is not None:
        fileio.write_summary2hdf5(setup['objectname'], samples, obs, obs_err, photbands, pars=results,
                                  grids=setup['grids'], filename=h5file)


def plot_results(setup, results, samples, constraints, gridnames, obs, obs_err, photbands):

    # check for 10 possible plots. Should be enough for now.
    for i in range(10):

        pindex = 'plot' + str(i)
        if not pindex in setup: continue

        if setup[pindex]['type'] == 'sed_fit':

            res = setup[pindex].get('result', 'best')

            pl.figure(i)
            pl.subplots_adjust(wspace=0.25)
            plotting.plot_fit(obs, obs_err, photbands, pars=results, constraints=constraints, grids=setup['grids'],
                              gridnames=gridnames, result=res)

            if not setup[pindex].get('path', None) is None:
                pl.savefig(setup[pindex].get('path', 'sed_fit.png'))

        if setup[pindex]['type'] == 'constraints':

            pl.figure(i, figsize=(2 * len(constraints), 6))
            pl.subplots_adjust(wspace=0.40, left=0.07, right=0.98)

            plotting.plot_constraints(constraints, samples, results)

            if not setup[pindex].get('path', None) is None:
                pl.savefig(setup[pindex].get('path', 'constraints.png'))

        if setup[pindex]['type'] == 'distribution':

            pars1 = []
            for p in setup[pindex].get('parameters', ['teff', 'rad', 'L', 'd']):
                if p in samples.dtype.names: pars1.append(p)

            data = repack_fields(samples[pars1])

            if setup[pindex].get('show_best', False):
                truths = [results[p][0] for p in data.dtype.names]
            else:
                truths = None

            fig = corner.corner(data.view(np.float64).reshape(data.shape + (-1,)),
                                labels=data.dtype.names,
                                quantiles=setup[pindex].get('quantiles', [0.025, 0.16, 0.5, 0.84, 0.975]),
                                levels=setup[pindex].get('levels', [0.393, 0.865, 0.95]),
                                truths=truths,
                                show_titles=True, title_kwargs={"fontsize": 12}, )

            if not setup[pindex].get('path', None) is None:
                pl.savefig(setup[pindex].get('path', 'distribution.png'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', action="store", type=str, help='use setup given in this file')
    parser.add_argument("-empty", dest='empty', type=str, default=None,
                        help="Create empty setup file ('single' or 'double')")
    parser.add_argument('--phot', dest='photometry', action='store_true',
                        help='When creating a new setupfile, use this option to also download photometry from Vizier and Tap archives.')
    parser.add_argument('--noplot', dest='noplot', action='store_true',
                        help="Don't show any plots, only store to disk.")
    args, variables = parser.parse_known_args()

    if args.empty is not None:

        objectname = args.filename
        filename = objectname + '_single.yaml' if args.empty == 'single' else objectname + '_binary.yaml'

        plx, e_plx = photometry_query.get_parallax(objectname)

        out = default_single if args.empty == 'single' else default_double
        out = out.replace('<photfilename>', objectname + '.phot')
        out = out.replace('<objectname>', objectname)
        out = out.replace('<plx>', str(plx))
        out = out.replace('<e_plx>', str(e_plx))

        ofile = open(filename, 'w')
        ofile.write(out)
        ofile.close()

        if args.photometry:
            photometry = photometry_query.get_photometry(objectname, filename=objectname + '.phot')

        sys.exit()

    if args.filename is None:
        print("Nothing to do")
        sys.exit()

    # -- load the setup file
    setupfile = open(args.filename)
    setup = yaml.safe_load(setupfile)
    setupfile.close()

    # -- obtain the observations
    photbands, obs, obs_err = get_observations(setup)

    # -- perform the SED fit
    results, samples, constraints, gridnames = fit_sed(setup, photbands, obs, obs_err)

    # -- write the results
    write_results(setup, results, samples, obs, obs_err, photbands)

    # -- create plots
    plot_results(setup, results, samples, constraints, gridnames, obs, obs_err, photbands)

    print("================================================================================")
    print("")
    print("Resulting parameter values and errors:")
    print("   Par             Best        Pc       emin       emax")
    for p in samples.dtype.names:
        print("   {:10s} = {}   {}   -{}   +{}".format(p, *plotting.format_parameter(p, results[p])))

    if not args.noplot:
        pl.show()


if __name__ == "__main__":
    main()
