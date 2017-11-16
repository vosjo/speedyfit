
import sys
import yaml
import argparse
import numpy as np
import pylab as pl

import mcmc, model, plotting, fileio

from ivs.io import ascii

default = """
# photometry file with index to the columns containing the photbands, observations and errors
photometryfile: path/to/file.dat
photband_index: 0
obs_index: 1
err_index: 2
# parameters to fit and the limits on them in same order as parameters
pnames: [teff, logg, rad, teff2, logg2, rad2, ebv]
limits:
- [3500, 6000]
- [3.5, 5.0]
- [0.7, 1.5]
- [20000, 40000]
- [4.5, 6.5]
- [0.05, 0.3]
- [0, 0.02]
# constraints on distance and mass ratio if known
constraints:
  q: [3.03, 0.2]
  distance: [600, 50] # in parsec
# added constraints on derived properties as mass, luminosity, luminosity ratio
derived_limits:
  mass: [0.5, 1.0]
  mass2: [0.1, 1.0]
# path to the model grids with integrated photometry
grids: 
- path/to/grid/1.fits
- path/to/grid/2.fits
# setup for the MCMC algorithm
nwalkers: 100    # total number of walkers
nsteps: 2000     # steps taken by each walker (not including burn-in)
nrelax: 500      # burn-in steps taken by each walker
a: 10            # relative size of the steps taken
# output options
datafile: none   # filepath to write results of all walkers
"""

if __name__=="__main__":

   parser = argparse.ArgumentParser()
   parser.add_argument("-f", type=str, dest='filename', default=None,
                       help="use setup given in filename")
   parser.add_argument("-empty", type=str, dest='empty', default=None,
                       help="When used, create an empty setup file with given filename")
   args, variables = parser.parse_known_args()
   
   if not args.empty is None:
      
      ofile = open(args.empty, 'w')
      ofile.write(default)
      ofile.close()
      
      sys.exit()
   
   setupfile = open(args.filename)
   setup = yaml.safe_load(setupfile)
   setupfile.close()
   
   #-- parse photometry
   data = ascii.read2array(setup['photometryfile'], dtype=str).T
   photbands = data[setup['photband_index']]
   obs = np.array(data[setup['obs_index']], dtype=float)
   obs_err = np.array(data[setup['err_index']], dtype=float)
   
   #-- pars limits
   pnames = setup['pnames']
   limits = np.array(setup['limits'])
   
   #-- pars constraints
   constraints = setup['constraints']
   
   #-- pars limits on derived properties
   derived_limits = setup['derived_limits']
   
   #-- pars grid
   gridnames = setup['grids']
   
   grids = []
   for i, name in enumerate(gridnames):
      
      ind = '' if i == 0 else str(i + 1)
      
      axis_values, grid_pars, pixelgrid, grid_names = model.prepare_grid(photbands, name,
            teffrange=limits[pnames.index('teff'+ind)],
            loggrange=limits[pnames.index('logg'+ind)],
            ebvrange =limits[pnames.index('ebv')],
            variables=['teff','logg','ebv'])

      grids.append([axis_values, pixelgrid])
   
   #-- pars mcmc setup
   nwalkers = setup.get('nwalkers', 100)
   nsteps = setup.get('nsteps', 2000)
   nrelax = setup.get('nrelax', 500)
   a = setup.get('a', 10)
   
   #-- MCMC
   results, samples = mcmc.MCMC(obs, obs_err, photbands, 
                                 pnames, limits, grids, 
                                 constraints=constraints, derived_limits=derived_limits,
                                 nwalkers=nwalkers, nsteps=nsteps, nrelax=nrelax,
                                 a=a, percentiles=[16, 50, 84])
   
   datafile = setup.get('datafile', None)
   if not datafile is None:
      fileio.write2fits(samples, datafile, setup=setup)
   
   
   pars = {}
   for par, v in zip(pnames, results):
      pars[par] = v[0]
   plotting.plot_fit(obs, obs_err, photbands, pars=pars, constraints=constraints)
   pl.show()