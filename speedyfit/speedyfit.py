
import sys
import yaml
import argparse
import numpy as np
import pylab as pl
import corner

from astropy.io import ascii

from numpy.lib.recfunctions import append_fields, repack_fields

from . import mcmc, model, plotting, fileio, filters

default_single = """
# photometry file with index to the columns containing the photbands, observations and errors
objectname: <objectname>
photometryfile: <photfilename>
photband_index: band
obs_index: flux
err_index: eflux
# parameters to fit and the limits on them in same order as parameters
pnames: [teff, logg, rad, ebv]
limits:
- [20000, 60000]
- [5.80, 5.80]
- [0.01, 0.5]
- [0, 0.10]
# constraints on distance and mass ratio is known
constraints: 
  parallax: [<plx>, <e_plx>]
# added constraints on derived properties as mass, luminosity, luminosity ratio
derived_limits: {}
# path to the model grids with integrated photometry
grids: 
- tmap
# setup for the MCMC algorithm
nwalkers: 100    # total number of walkers
nsteps: 1000     # steps taken by each walker (not including burn-in)
nrelax: 250      # burn-in steps taken by each walker
a: 10            # relative size of the steps taken
# set the percentiles for the error determination 
percentiles: [0.2, 50, 99.8] # 16 - 84 corresponds to 1 sigma
# output options
resultfile: <objectname>_results_single.csv   # filepath to write results
plot1:
 type: sed_fit
 result: pc
 path: <objectname>_sed_single.png
plot2:
 type: distribution
 show_best: true
 path: <objectname>_distribution_single.png
 parameters: ['teff', 'rad', 'L', 'ebv', 'd', 'mass']
"""

default_double = """
# photometry file with index to the columns containing the photbands, observations and errors
objectname: <objectname>
photometryfile: <photfilename>
photband_index: band
obs_index: flux
err_index: eflux
# parameters to fit and the limits on them in same order as parameters
pnames: [teff, logg, rad, teff2, logg2, rad2, ebv]
limits:
- [3500, 10000]
- [4.31, 4.31]
- [0.01, 2.5]
- [20000, 50000]
- [5.8, 5.8]
- [0.01, 0.5]
- [0, 0.10]
# constraints on distance and mass ratio if known
constraints:
  parallax: [<plx>, <e_plx>]
# added constraints on derived properties as mass, luminosity, luminosity ratio
derived_limits: {}
# path to the model grids with integrated photometry
grids: 
- kurucz2
- tmap
# setup for the MCMC algorithm
nwalkers: 100    # total number of walkers
nsteps: 500     # steps taken by each walker (not including burn-in)
nrelax: 250      # burn-in steps taken by each walker
a: 10            # relative size of the steps taken
# set the percentiles for the error determination 
percentiles: [16, 50, 84] # 16 - 84 corresponds to 1 sigma
# output options
resultfile: <objectname>_results_binary.csv   # filepath to write results
plot1:
 type: sed_fit
 path: <objectname>_sed_binary.png
plot2:
 type: distribution
 path: <objectname>_distribution_primary.png
 parameters: ['teff', 'rad', 'teff2', 'rad2', 'ebv', 'd']
#plot3:
# type: distribution
# #path: <objectname>_distribution_derived.png
# parameters: ['mass', 'L', 'mass2', 'L2', 'q']
"""


def main():

   parser = argparse.ArgumentParser()
   parser.add_argument('filename', action="store", type=str, help='use setup given in this file')
   parser.add_argument("-empty", dest='empty', type=str, default=None,
                       help="Create empty setup file ('single' or 'double')")
   parser.add_argument('--phot', dest='photometry', action='store_true', help='When creating a new setupfile, use this option to also download photometry from Vizier and Tap archives.')
   parser.add_argument('--noplot', dest='noplot', action='store_true',
                       help="Don't show any plots, only store to disk.")
   args, variables = parser.parse_known_args()
   
   if args.empty is not None:
      
      from . import photometry_query
      
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
         photometry = photometry_query.get_photometry(objectname, filename =objectname + '.phot')
      
      sys.exit()
   
   if args.filename is None:
      
      print("Nothing to do")
      sys.exit()
   
   #-- load the setup file
   setupfile = open(args.filename)
   setup = yaml.safe_load(setupfile)
   setupfile.close()
   
   #-- parse photometry
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

   nani = np.isnan(obs) | np.isnan(obs_err)
   if any(nani):
      print("Warning: there are NaN values in the following photometric bands:")
      for p in photbands[nani]:
         print("\t {}".format(p))
   obs, obs_err, photbands = obs[~nani], obs_err[~nani], photbands[~nani]
   
   #-- remove colors
   color = np.array([filters.is_color(p) for p in photbands])
   s = np.where(~color)
   photbands, obs, obs_err = photbands[s], obs[s], obs_err[s]
   
   #-- pars limits
   pnames = setup['pnames']
   limits = np.array(setup['limits'])
   
   #-- pars constraints
   constraints = setup['constraints']
   for con, val in list(constraints.items()):
      if len(val) == 2:
         constraints[con] = [val[0], val[1], val[1]]
   
   if 'parallax' in constraints:
      p, pm, pp = constraints.pop('parallax')
      constraints['distance'] = [1000./p, 1000.*pm/p**2, 1000.*pp/p**2]

   print ("Applied constraints: ")
   for con, val in list(constraints.items()):
      print("\t {} = {} - {} + {}".format(con, val[0], val[1], val[2]))

   if 'distance' in constraints:
      # convert pc to Rsol
      constraints['distance'] = [44365810.04823812 * constraints['distance'][0], 
                                 44365810.04823812 * constraints['distance'][1],
                                 44365810.04823812 * constraints['distance'][2],]

   #-- pars limits on derived properties
   derived_limits = setup['derived_limits']
   
   #-- pars grid
   gridnames = setup['grids']
   grids = model.load_grids(gridnames, pnames, limits, photbands)
   
   #-- switch logg to g for a binary system
   if 'q' in constraints:
      if 'logg' in pnames:
         limits[pnames.index('logg')] = 10**limits[pnames.index('logg')]
         pnames[pnames.index('logg')] = 'g'
      
      if 'logg2' in pnames:
         limits[pnames.index('logg2')] = 10**limits[pnames.index('logg2')]
         pnames[pnames.index('logg2')] = 'g2'
      
      if 'logg' in constraints:
         g = 10**constraints['logg'][0]
         g_el = g * constraints['logg'][1] * np.log(10)
         g_eu = g * constraints['logg'][2] * np.log(10)
         constraints['g'] = [g, g_el, g_eu]
         
      if 'logg2' in constraints:
         g = 10**constraints['logg2'][0]
         g_el = g * constraints['logg2'][1] * np.log(10)
         g_eu = g * constraints['logg2'][2] * np.log(10)
         constraints['g2'] = [g, g_el, g_eu]
   
   #-- check for variables that are kept fixed
   fixed = np.where(limits[:,0] == limits[:,1])
   varied = np.where(limits[:,0] != limits[:,1])
   
   pnames = np.array(pnames)
   fixed_variables = {}
   for par, val in zip(pnames[fixed], limits[:,0][fixed]):
      fixed_variables[par] = val
      
   pnames = list(pnames[varied])
   limits = limits[varied]
   
   #-- pars mcmc setup
   nwalkers = setup.get('nwalkers', 100)
   nsteps = setup.get('nsteps', 2000)
   nrelax = setup.get('nrelax', 500)
   a = setup.get('a', 10)
   percentiles = setup.get('percentiles', [16, 50, 84])
   
   
   #-- MCMC
   results, samples = mcmc.MCMC(obs, obs_err, photbands,
                                pnames, limits, grids,
                                fixed_variables=fixed_variables,
                                constraints=constraints, derived_limits=derived_limits,
                                nwalkers=nwalkers, nsteps=nsteps, nrelax=nrelax,
                                a=a)
   
   #-- add fixed variables to results dictionary
   for par, val in list(fixed_variables.items()):
      results[par] = [val, val, 0, 0]
   
   #-- deal with the switch back to logg
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
   
   print("================================================================================")
   print("")
   print("Resulting parameter values and errors:")
      
   pc  = np.percentile(samples.view(np.float64).reshape(samples.shape + (-1,)), percentiles, axis=0)
   pars = {}
   for p, v, e1, e2 in zip(samples.dtype.names, pc[1], pc[1]-pc[0], pc[2]-pc[1]):
      results[p] = [results[p], v, e1, e2]
      pars[p] = v
   
   ##-- calculate the Chi2 of the 50th percentile results
   #y_syn, extra_drv = model.get_itable(**pars)
   #chi2, _, _ = stat_func(obs,
                         #obs_err,
                         #colors, y_syn, pars,
                         #constraints=constraints)
   #results['chi2'] = [results['chi2'][0], chi2, results['chi2'][2], results['chi2'][3]]
   
   
   print("   Par             Best        Pc       emin       emax")
   for p in samples.dtype.names:
      print("   {:10s} = {}   {}   -{}   +{}".format(p, *plotting.format_parameter(p, results[p])))
   
   
   # out = ""
   # out += "{:0.0f}\t{:0.0f}\t".format(results['teff'][1], np.average([results['teff'][2],results['teff'][3]]))
   # for par in ['logg', 'L', 'rad']:
   #    out += "{:0.3f}\t{:0.3f}\t".format(results[par][1],
   #                                       np.average([results[par][2],results[par][3]]))
   # if 'teff2' in results:
   #    out += "{:0.0f}\t{:0.0f}\t".format(results['teff2'][1], np.average([results['teff2'][2],results['teff2'][3]]))
   #    for par in ['logg2', 'L2', 'rad2']:
   #       out += "{:0.3f}\t{:0.3f}\t".format(results[par][1],
   #                                        np.average([results[par][2],results[par][3]]))
   # out += "{:0.0f}\t{:0.0f}\t".format(results['d'][1], np.average([results['d'][2],results['d'][3]]))
   # print out

   outpars, outvals = [], []
   for par in ['teff', 'logg', 'L', 'rad']:
      outpars.append(par)
      outpars.append(par+'_err')
      outvals.append(results[par][1])
      outvals.append(np.average([results[par][2],results[par][3]]))

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
      
      #-- get plain text settings to retain comments and ordering
      setupfile = open(args.filename)
      setup_str = "".join(setupfile.readlines())
      setupfile.close()
      
      fileio.write2fits(samples, datafile, setup=setup_str)


   h5file = setup.get('h5file', None)
   if h5file is not None:
      fileio.write_summary2hdf5(setup['objectname'], samples, obs, obs_err, photbands, pars=results,
                                grids=setup['grids'], filename=h5file)
   
   #-- Plotting 
   
   for i in range(10):
      
      pindex = 'plot'+str(i)
      if not pindex in setup: continue
      
      if setup[pindex]['type'] == 'sed_fit':
         
         res = setup[pindex].get('result', 'best')
         
         pl.figure(i)
         pl.subplots_adjust(wspace=0.25)
         plotting.plot_fit(obs, obs_err, photbands, pars=results, constraints=constraints, grids=setup['grids'], gridnames=gridnames, result=res)
         
         if not setup[pindex].get('path', None) is None:
            pl.savefig(setup[pindex].get('path', 'sed_fit.png'))
   
      
      if setup[pindex]['type'] == 'constraints':
         
         pl.figure(i, figsize=(2*len(constraints), 6) )
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
                       labels = data.dtype.names,
                       quantiles=setup[pindex].get('quantiles', [0.025, 0.16, 0.5, 0.84, 0.975]),
                       levels=setup[pindex].get('levels', [0.393, 0.865, 0.95]),
                       truths=truths,
                       show_titles=True, title_kwargs={"fontsize": 12},)
         
         if not setup[pindex].get('path', None) is None:
            pl.savefig(setup[pindex].get('path', 'distribution.png'))

   if not args.noplot:
      pl.show()


if __name__ == "__main__":
    main()
