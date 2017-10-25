import sys
import numpy as np

import emcee

import filters, statfunc, model

from ivs.io import ascii
 
def lnlike(theta, derived_properties, y, yerr, **kwargs):
   """
   log likelihood function
   
   Calculates the chi2 of the model defined by theta compared to the observed magnitudes
   and colors. Will also take possible constraints on q, lr and d into account.
   """
   model_func = kwargs.pop('model_func', model.get_itable)
   stat_func = kwargs.pop('stat_func', statfunc.stat_chi2)
   colors = kwargs.get('colors', [False for i in y])
   constraints = kwargs.pop('constraints', {})
   
   #-- create keyword parameters from theta
   pars = {}
   for name, value in zip(kwargs['pnames'], theta):
      pars[name]=value
   
   
   #-- calculate synthetic magnitudes **kwargs contains infor about which grid to use
   kwargs.update(pars)
   y_syn, Labs = model_func(**kwargs)
   
   
   chi2, scales, e_scales = stat_func(y,
                                      yerr,
                                      colors, y_syn, 
                                      constraints_syn=derived_properties,
                                      **constraints)
   
   return np.log( np.exp(-chi2/2.) )
   
def lnprior(theta, derived_properties, limits, **kwargs):
   """
   Simple uniform (flat) prior on all parameters if the parameters 
   are within their range, and the derived properties (q, m, ..) are also 
   within their limits.
   
   if all parameters are within the provided limits, the the returned 
   log probability is 0, otherwise it is -inf.
   
   :param theta: list of model parameters
   :type theta: list
   :param limits: limits on the model parameters
   :type limits: list of tuples
   
   :return: logarithm of the probability of the parameters (theta) given the 
            model limits
   :rtype: float
   """
   
   derived_limits = kwargs.pop('derived_limits', {})
   
   #-- check if all parameters are within their limits
   if any(theta < limits[:,0]) or any(theta > limits[:,1]):
      return -np.inf
      
   #-- check that all derived properties are within limits
   for lim in derived_limits.keys():
      if derived_properties[lim] < derived_limits[lim][0] or\
         derived_properties[lim] > derived_limits[lim][1]:
         return -np.inf
   
   return 0
   
def lnprob(theta, y, yerr, limits, **kwargs):
   """
   full log probability function combining the prior and the likelihood
   
   will return -inf if any of :py:func:`lnprior` or :py:func:`lnlikelyhood` is 
   infite, otherwise it will return the sum of both functions.
   
   :param theta: list of model parameters (normaly mass, fe/h and age)
   :type theta: list
   :param y: 1D array of observables
   :type y: array
   :param yerr: 1D array containing errors on every observable
   :type yerr: array
   :param limits: limits on the model parameters
   :type limits: list of tuples
   
   :return: the sum of the log prior and log likelihood
   :rtype: float
   """
   
   syn_drv = statfunc.get_derived_properties(theta, kwargs['pnames'])
   
   lp = lnprior(theta, syn_drv, limits, **kwargs)
   if not np.isfinite(lp):
      return -np.inf
   
   ll = lnlike(theta, syn_drv, y, yerr, **kwargs)
   if not np.isfinite(ll):
      return -np.inf
   
   return lp + ll

def MCMC(photfilename):
   
   #-- setup parameter and limits
   pnames = ['teff', 'logg', 'rad', 'teff2', 'logg2', 'rad2', 'ebv']
   limits = [[4000, 7000],
             [2.5, 4.0],
             [0.5, 5.0],
             [20000, 40000],
             [5.0, 6.5],
             [0.05, 0.5],
             [0, 0.02]]
   limits = np.array(limits)
   
   d = 1000 / 1.475 * 44365810 # in Rsun
   constraints = {'q':(2.35, 0.3), 'distance':(d, d/4)}
   
   #-- load observed photometry
   master = ascii.read2recarray(photfilename)
   photbands = master['photband']
   colors = np.array([filters.is_color(photband) for photband in photbands],bool)
   
   obs = master['cmeas']
   obs_err = master['e_cmeas']
   
   #-- Prepare the grids
   gridfilename = '/home/joris/Python/ivsdata/sedtables/modelgrids/ikurucz93_z0.0_k2odfnew_sed_lawfitzpatrick2004_Rv3.10.fits'
   
   axis_values, grid_pars, pixelgrid, grid_names = model.prepare_grid(photbands, gridfilename,
               teffrange=(4000, 7000),loggrange=(2.5, 4.0),
               ebvrange=(0.0, 0.02),
               variables=['teff','logg','ebv'])
   
   grid1 = [axis_values, pixelgrid]
   
   gridfilename = '/home/joris/Python/ivsdata/sedtables/modelgrids/iTMAP2012_sdB_extended_lawfitzpatrick2004_Rv3.10.fits'
      
   axis_values, grid_pars, pixelgrid, grid_names = model.prepare_grid(photbands, gridfilename,
               teffrange=(20000, 40000),loggrange=(5.0, 6.5),
               ebvrange=(0.0, 0.02),
               variables=['teff','logg','ebv'])
   
   grid2 = [axis_values, pixelgrid]
   
   grids = [grid1, grid2]
   
   #-- initialize the walkers
   nwalkers = 100
   wlimits = [[4000, 5500],
             [2.7, 3.2],
             [0.75, 2.0],
             [26000, 31000],
             [5.5, 6.0],
             [0.10, 0.20],
             [0, 0.02]]
   pos = [ np.random.uniform(lim[0], lim[1], nwalkers) for lim in limits]
   pos = np.array(pos).T
   
   #-- setup the sampler
   ndim = len(pnames)
   sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, a=10, args=(obs, obs_err, limits), kwargs={'pnames':pnames, 'colors':colors, 'grid':grids, 'constraints':constraints})
   
   
   #-- run the sampler
   nsteps = 2000
   for i, result in enumerate(sampler.sample(pos, iterations=nsteps)):
      if (i+1) % 100 == 0:
         print("{0:5.1%}".format(float(i) / nsteps))
   
   #-- remove first 50 steps and combine the results from the individual walkers 
   samples = sampler.chain[:, 750:, :].reshape((-1, ndim))
   
   percentiles=[16, 50, 84]
   pc  = np.percentile(samples, percentiles, axis=0)
   results = [(v, e1, e2) for v, e1, e2 in zip(pc[1], pc[1]-pc[0], pc[2]-pc[1])]
   
   print results
   return results, samples
   
   
if __name__=="__main__":
   
   results, samples = MCMC('BD-7_5977.phot')
   
   import pylab as pl
   import corner
   
   fig = corner.corner(samples, 
                       labels = ['teff', 'logg', 'rad', 'teff2', 'logg2', 'rad2', 'ebv'],
                       quantiles=[0.16, 0.5, 0.84],
                       show_titles=True, title_kwargs={"fontsize": 12})
                     
   pl.show()