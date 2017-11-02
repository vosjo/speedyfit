import sys
import numpy as np

import emcee

import statfunc, model

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
   y_syn, extra_drv = model_func(**kwargs)
   
   
   chi2, scales, e_scales = stat_func(y,
                                      yerr,
                                      colors, y_syn, 
                                      constraints_syn=derived_properties,
                                      **constraints)
   
   #-- add distance to extra derived parameter (which already contains luminosities)
   extra_drv['d'] = np.sqrt(1/scales)/44365810.04823812 
   
   return -chi2/2, extra_drv
   
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
   syn_drv['d'] = 0
   
   lp = lnprior(theta, syn_drv, limits, **kwargs)
   if not np.isfinite(lp):
      return -np.inf, syn_drv
   
   ll, extra_drv = lnlike(theta, syn_drv, y, yerr, **kwargs)
   syn_drv.update(extra_drv)
   if not np.isfinite(ll):
      return -np.inf, syn_drv
   
   return lp + ll, syn_drv
   

def MCMC(obs, obs_err, photbands, 
         pnames, limits, grids, 
         constraints={}, derived_limits={},
         nwalkers=100, nsteps=1000, nrelax=150, a=15, percentiles=[16, 50, 84]):
   
   #-- check which bands are colors
   colors = np.array([model.is_color(photband) for photband in photbands],bool)
   
   #-- initialize the walkers
   np.random.seed(1)
   pos = [ np.random.uniform(lim[0], lim[1], nwalkers) for lim in limits]
   pos = np.array(pos).T
   
   #-- setup the sampler
   ndim = len(pnames)
   kwargs = {'pnames':pnames, 
             'colors':colors, 
             'grid':grids, 
             'constraints':constraints, 
             'derived_limits':derived_limits}
   
   sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, a=a, 
                                   args=(obs, obs_err, limits), kwargs=kwargs)
   
   #================
   # MCMC part
   
   #-- burn in (let walkers relax before starting to store results)
   print "\nBurn In"
   for i, result in enumerate(sampler.sample(pos, iterations=nrelax, storechain=False)):
      if (i+1) % 100 == 0:
         print("{0:5.1%}".format(float(i) / nrelax))
   pos = result[0]
   
   print "\nRun"
   #-- run the sampler
   for i, result in enumerate(sampler.sample(pos, iterations=nsteps)):
      if (i+1) % 100 == 0:
         print("{0:5.1%}".format(float(i) / nsteps))
   
   
   #-- remove first nrelax steps and combine the results from the individual walkers 
   samples = sampler.flatchain
   blobs = np.array(sampler.blobs).flatten()
   probabilities = sampler.flatlnprobability
   
   #-- remove all steps that are not accepted (lnprob = -inf)
   accept = np.where(np.isfinite(probabilities))
   samples = samples[accept]
   blobs = blobs[accept]
   
   #-- calculate results
   pc  = np.percentile(samples, percentiles, axis=0)
   results = [(v, e1, e2) for v, e1, e2 in zip(pc[1], pc[1]-pc[0], pc[2]-pc[1])]
   results = np.array(results)
   
   #-- convert to recarrays
   dtypes = [(n, 'f8') for n in pnames]
   samples = np.array([tuple(s) for s in samples], dtype=dtypes)
   
   names = blobs[0].keys()
   pars = []
   for b in blobs:
      if len(b.keys()) < 7:
         print 'err', b
      pars.append(tuple([b[n] for n in names]))
   dtypes = [(n, 'f8') for n in names]
   blobs = np.array(pars, dtype=dtypes)
   
   
   return results, samples, blobs
   
   

   
   