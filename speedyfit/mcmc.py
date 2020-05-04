import numpy as np
from numpy.lib.recfunctions import merge_arrays

import emcee

from . import statfunc, model, filters

def lnlike(pars, derived_properties, y, yerr, **kwargs):
    """
    log likelihood function

    Calculates the chi2 of the model defined by theta compared to the observed magnitudes
    and colors. Will also take possible constraints on q, lr and d into account.
    """
    model_func = kwargs.pop('model_func', model.get_itable)
    stat_func = kwargs.pop('stat_func', statfunc.stat_chi2)
    colors = kwargs.get('colors', np.array([False for i in y], bool))
    constraints = kwargs.pop('constraints', {})


    #-- calculate synthetic magnitudes **kwargs contains info about which grid to use
    kwargs.update(pars)
    y_syn, extra_drv = model_func(**kwargs)
    if 'L' in extra_drv and 'L2' in extra_drv:
        #-- update light ratio with correct value from models
        extra_drv['lr'] = extra_drv['L'] / extra_drv['L2']
    derived_properties.update(extra_drv)


    chi2, scales, e_scales = stat_func(y,
                                       yerr,
                                       colors, y_syn, pars,
                                       derived_properties=derived_properties,
                                       constraints=constraints)

    #-- add distance to extra derived parameter (which already contains luminosities)
    #   distance is converted from Rsol to pc
    extra_drv['d'] = np.sqrt(1/scales)/44365810.04823812
    extra_drv['scale'] = scales
    extra_drv['chi2'] = chi2

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
    for lim in list(derived_limits.keys()):
        if derived_properties[lim] < derived_limits[lim][0] or \
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

    #-- create keyword parameters from theta
    pars = {}
    for name, value in zip(kwargs['pnames'], theta):
        pars[name]=value

    #-- add extra variables which are not fitted to pars.
    pars.update(kwargs.pop('fixed_variables', {}))

    #-- get derived properties
    prop_func = kwargs.pop('prop_func', statfunc.get_derived_properties)
    syn_drv = prop_func(**pars)

    if 'rad' in syn_drv:
        pars['rad']=syn_drv['rad']
    if 'rad2' in syn_drv:
        pars['rad2']=syn_drv['rad2']

    #-- calculate prior probability
    lp = lnprior(theta, syn_drv, limits, **kwargs)
    if not np.isfinite(lp):
        return -np.inf, syn_drv

    #-- calculate likelihood
    ll, extra_drv = lnlike(pars, syn_drv, y, yerr, **kwargs)
    syn_drv.update(extra_drv)
    if not np.isfinite(ll):
        return -np.inf, syn_drv

    return lp + ll, syn_drv


def MCMC(obs, obs_err, photbands,
         pnames, limits, grids,
         fixed_variables={}, constraints={}, derived_limits={},
         nwalkers=100, nsteps=1000, nrelax=150, a=10, pos=None):

    #-- check which bands are colors
    colors = np.array([filters.is_color(photband) for photband in photbands], bool)

    #-- initialize the walkers randomly is no starting positions are given
    if pos is None:
        pos = [ np.random.uniform(lim[0], lim[1], nwalkers) for lim in limits]
        pos = np.array(pos).T
    else:
        nwalkers = pos.shape[0]

    #-- setup the sampler
    ndim = len(pnames)
    kwargs = {'pnames':pnames,
              'colors':colors,
              'grid':grids,
              'fixed_variables':fixed_variables,
              'constraints':constraints,
              'derived_limits':derived_limits,
              'prop_func': statfunc.get_derived_properties}

    # TODO: storing the blobs as dictionary with dtype object and then later converting to recarray is inefficient.
    # This needs to be addressed: provide the correct dtypes here and let emcee directly store them in recarray
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, a=a,
                                    args=(obs, obs_err, limits), kwargs=kwargs, blobs_dtype=[('blob','O')])

    #================
    # MCMC part

    #-- run the sampler, both burn in and actual run in one
    sampler.run_mcmc(pos, nsteps+nrelax, progress=True)

    #-- combine the results from the individual walkers discarding the burn in steps
    samples = sampler.get_chain(discard=nrelax, thin=1, flat=True)
    blobs = sampler.get_blobs(discard=nrelax, thin=1, flat=True)
    probabilities = sampler.get_log_prob(discard=nrelax, thin=1, flat=True)

    #-- clear the samples to save memory
    sampler.reset()

    #-- remove all steps that are not accepted (lnprob == -inf)
    accept = np.where(np.isfinite(probabilities))
    samples = samples[accept]
    blobs = blobs[accept]
    probabilities = probabilities[accept]

    if len(samples) == 0:
        raise ValueError('No models were accepted, all probabilities were -inf')

    #-- convert to recarrays
    dtypes = [(n, 'f8') for n in pnames]
    samples = np.array([tuple(s) for s in samples], dtype=dtypes)

    names = list(blobs[0][0].keys())
    pars = []
    for b in blobs:
        pars.append(tuple([b[0][n] for n in names]))
    dtypes = [(n, 'f8') for n in names]
    blobs = np.array(pars, dtype=dtypes)

    #-- remove all steps where model creation failed (d == 0)
    accept = np.where(blobs['d'] > 0)
    samples = samples[accept]
    blobs = blobs[accept]
    probabilities = probabilities[accept]

    #-- merge all results in 1 recarray and select best model
    data = merge_arrays((samples, blobs), asrecarray=True, flatten=True)
    best = np.where(probabilities == np.max(probabilities))

    results = {}
    for n, v in zip(data.dtype.names, data[best][0]):
        results[n] = v


    return results, data
   
   

   
   
