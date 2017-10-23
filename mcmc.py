 
import numpy as np
 
def get_derived_properties(theta, pnames):
   """
   Function that will derive several properties based on the chosen models
   
   Currently the following properties are calculated:
   - mass
   - mass1
   - mass2
   - q
   
   returns dictionary of all properties that could be calculated.
   """
   
   derived_properties = {}
   
   GG = 6.67384e-08
   Rsol = 69550800000.0
   Msol = 1.988547e+33
   
   #-- derive masses
   if 'rad' in pnames and 'logg' in pnames:
      mass = 10**theta[pnames.index('logg')] * (theta[pnames.index('rad')] * Rsol)**2 / GG
      derived_properties['mass'] = mass / Msol
   
   if 'rad1' in pnames and 'logg1' in pnames:
      mass = 10**theta[pnames.index('logg1')] * (theta[pnames.index('rad1')] * Rsol)**2 / GG
      derived_properties['mass1'] = mass / Msol
      
   if 'rad2' in pnames and 'logg2' in pnames:
      mass = 10**theta[pnames.index('logg2')] * (theta[pnames.index('rad2')] * Rsol)**2 / GG
      derived_properties['mass2'] = mass / Msol
   
   #-- derive mass ratio
   if 'rad' in pnames and 'logg' in pnames and 'rad2' in pnames and 'logg2' in pnames:
      m1 = ( theta[pnames.index('rad')]**2 * 10**theta[pnames.index('logg')] )
      m2 = ( theta[pnames.index('rad2')]**2 * 10**theta[pnames.index('logg2')] )
      derived_properties['q'] = m1 / m2
      
   if 'rad1' in pnames and 'logg1' in pnames and 'rad2' in pnames and 'logg2' in pnames:
      m1 = ( theta[pnames.index('rad1')]**2 * 10**theta[pnames.index('logg1')] )
      m2 = ( theta[pnames.index('rad2')]**2 * 10**theta[pnames.index('logg2')] )
      derived_properties['q'] = m1 / m2
   
   #-- derive luminocity ratio
   if 'rad' in pnames and 'teff' in pnames and 'rad2' in pnames and 'teff2' in pnames:
      l1 = ( theta[pnames.index('rad')]**2 * theta[pnames.index('teff')]**4 )
      l2 = ( theta[pnames.index('rad2')]**2 * theta[pnames.index('teff2')]**4 )
      derived_properties['lr'] = l1 / l2
      
   if 'rad1' in pnames and 'teff1' in pnames and 'rad2' in pnames and 'teff2' in pnames:
      l1 = ( theta[pnames.index('rad1')]**2 * theta[pnames.index('teff1')]**4 )
      l2 = ( theta[pnames.index('rad2')]**2 * theta[pnames.index('teff2')]**4 )
      derived_properties['lr'] = l1 / l2
   
   return derived_properties
 
def lnlike(theta, derived_properties, y, yerr, **kwargs):
   """
   log likelihood function
   
   Calculates the chi2 of the model defined by theta compared to the observed magnitudes
   and colors. Will also take possible constraints on q, lr and d into account.
   """
   
   model_func = kwargs.pop('model_func',model.get_itable_pix)
   stat_func = kwargs.pop('stat_func',stat_chi2)
   photbands = kwargs.pop(photbands)
   constraints = kwargs.pop(constraints, {})
   grids = kwargs.pop('grids', [])
   
   #-- check which filter is a color
   colors = np.array([filters.is_color(photband) for photband in photbands],bool)
   
   #-- calculate synthetic magnitudes **kwargs contains infor about which grid to use
   y_syn = model_func(theta, photbands=photbands, **kwargs)
   
   
   chisqs,scales,e_scales = stat_func(y,
                                      yerr,
                                      colors, syn_flux, 
                                      constraints_syn=constraints_syn,
                                      **constraints)
   
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
   
   syn_drv = get_derived_properties(theta, kwargs['pnames'])
   
   lp = lnprior(theta, syn_drv, limits, **kwargs)
   if not np.isfinite(lp):
      return -np.inf
   
   ll = lnlike(theta, syn_drv, y, yerr)
   if not np.isfinite(ll):
      return -np.inf
   
   return lp + ll