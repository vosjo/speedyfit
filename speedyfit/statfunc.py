# -*- coding: utf-8 -*-
import numpy as np

import warnings


def get_derived_properties(**pars):
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

    # -- derive masses
    if 'rad' in pars and 'logg' in pars:
        mass = 10 ** pars['logg'] * (pars['rad'] * Rsol) ** 2 / GG
        derived_properties['mass'] = mass / Msol

    if 'rad' in pars and 'g' in pars:
        mass = pars['g'] * (pars['rad'] * Rsol) ** 2 / GG
        derived_properties['mass'] = mass / Msol

    if 'rad2' in pars and 'logg2' in pars:
        mass = 10 ** pars['logg2'] * (pars['rad2'] * Rsol) ** 2 / GG
        derived_properties['mass2'] = mass / Msol

    if 'rad2' in pars and 'g2' in pars:
        mass = pars['g2'] * (pars['rad2'] * Rsol) ** 2 / GG
        derived_properties['mass2'] = mass / Msol

        # -- derive radii
    if 'mass' in pars and 'logg' in pars:
        rad = np.sqrt(GG * pars['mass'] * Msol / 10 ** pars['logg'])
        derived_properties['rad'] = rad / Rsol

    if 'mass' in pars and 'g' in pars:
        rad = np.sqrt(GG * pars['mass'] * Msol / pars['g'])
        derived_properties['rad'] = rad / Rsol

    if 'mass2' in pars and 'logg2' in pars:
        rad = np.sqrt(GG * pars['mass2'] * Msol / 10 ** pars['logg2'])
        derived_properties['rad2'] = rad / Rsol

    if 'mass2' in pars and 'g2' in pars:
        rad = np.sqrt(GG * pars['mass2'] * Msol / pars['g2'])
        derived_properties['rad2'] = rad / Rsol

    # -- derive mass ratio
    if 'mass' in derived_properties and 'mass2' in derived_properties:
        derived_properties['q'] = derived_properties['mass'] / derived_properties['mass2']

    if 'mass' in pars and 'mass2' in pars:
        derived_properties['q'] = pars['mass'] / pars['mass2']

    # -- derive luminosity ratio
    if 'rad' in pars and 'teff' in pars and 'rad2' in pars and 'teff2' in pars:
        l1 = (pars['rad'] ** 2 * pars['teff'] ** 4)
        l2 = (pars['rad2'] ** 2 * pars['teff2'] ** 4)
        derived_properties['lr'] = l1 / l2

    if 'mass' in pars and 'teff' in pars and 'mass2' in pars and 'teff2' in pars:
        l1 = (derived_properties['rad'] ** 2 * pars['teff'] ** 4)
        l2 = (derived_properties['rad2'] ** 2 * pars['teff2'] ** 4)
        derived_properties['lr'] = l1 / l2

    # -- derive radius ratio
    if 'rad' in pars and 'rad2' in pars:
        derived_properties['rr'] = pars['rad'] / pars['rad2']

    if 'rad' in derived_properties and 'rad2' in derived_properties:
        derived_properties['rr'] = derived_properties['rad'] / derived_properties['rad2']

    # -- add empty values for luminosity and distance to prevent problems with
    #   failed models
    derived_properties.update({'d': 0, 'L': 0})
    if 'rad2' in pars:
        derived_properties['L2'] = 0
    if hasattr(pars['teff'], '__iter__'):
        derived_properties['d'] = np.zeros_like(pars['teff'])
        derived_properties['L'] = np.zeros_like(pars['teff'])
        if 'rad2' in pars:
            derived_properties['L2'] = np.zeros_like(pars['teff'])

    return derived_properties


# def get_derived_properties(theta, pnames):
# """
# Function that will derive several properties based on the chosen models

# Currently the following properties are calculated:
# - mass
# - mass1
# - mass2
# - q

# returns dictionary of all properties that could be calculated.
# """

# derived_properties = {}

# GG = 6.67384e-08
# Rsol = 69550800000.0
# Msol = 1.988547e+33

##-- derive masses
# if 'rad' in pnames and 'logg' in pnames:
# mass = 10**theta[pnames.index('logg')] * (theta[pnames.index('rad')] * Rsol)**2 / GG
# derived_properties['mass'] = mass / Msol

# if 'rad' in pnames and 'g' in pnames:
# mass = theta[pnames.index('g')] * (theta[pnames.index('rad')] * Rsol)**2 / GG
# derived_properties['mass'] = mass / Msol

# if 'rad2' in pnames and 'logg2' in pnames:
# mass = 10**theta[pnames.index('logg2')] * (theta[pnames.index('rad2')] * Rsol)**2 / GG
# derived_properties['mass2'] = mass / Msol

# if 'rad2' in pnames and 'g2' in pnames:
# mass = theta[pnames.index('g2')] * (theta[pnames.index('rad2')] * Rsol)**2 / GG
# derived_properties['mass2'] = mass / Msol

##-- derive radii
# if 'mass' in pnames and 'logg' in pnames:
# rad = np.sqrt(GG * theta[pnames.index('mass')] * Msol / 10**theta[pnames.index('logg')])
# derived_properties['rad'] = rad / Rsol

# if 'mass' in pnames and 'g' in pnames:
# rad = np.sqrt(GG * theta[pnames.index('mass')] * Msol / theta[pnames.index('g')])
# derived_properties['rad'] = rad / Rsol

# if 'mass2' in pnames and 'logg2' in pnames:
# rad = np.sqrt(GG * theta[pnames.index('mass2')] * Msol / 10**theta[pnames.index('logg2')])
# derived_properties['rad2'] = rad / Rsol

# if 'mass2' in pnames and 'g2' in pnames:
# rad = np.sqrt(GG * theta[pnames.index('mass2')] * Msol / theta[pnames.index('g2')])
# derived_properties['rad2'] = rad / Rsol

##-- derive mass ratio
# if 'mass' in derived_properties and 'mass2' in derived_properties:
# derived_properties['q'] = derived_properties['mass'] / derived_properties['mass2']

# if 'mass' in pnames and 'mass2' in pnames:
# derived_properties['q'] = theta[pnames.index('mass')] / theta[pnames.index('mass2')]

##-- derive luminosity ratio
# if 'rad' in pnames and 'teff' in pnames and 'rad2' in pnames and 'teff2' in pnames:
# l1 = ( theta[pnames.index('rad')]**2 * theta[pnames.index('teff')]**4 )
# l2 = ( theta[pnames.index('rad2')]**2 * theta[pnames.index('teff2')]**4 )
# derived_properties['lr'] = l1 / l2

# if 'mass' in pnames and 'teff' in pnames and 'mass2' in pnames and 'teff2' in pnames:
# l1 = ( derived_properties['rad']**2 * theta[pnames.index('teff')]**4 )
# l2 = ( derived_properties['rad2']**2 * theta[pnames.index('teff2')]**4 )
# derived_properties['lr'] = l1 / l2

##-- add empty values for luminosity and distance to prevent problems with
##   failed models
# derived_properties.update({'d':0, 'L':0})
# if 'rad2' in pnames:
# derived_properties['L2'] = 0

# return derived_properties

def get_derived_properties_binary(theta, pnames):
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

    # -- derive masses
    if 'g' in pnames:
        m1 = theta[pnames.index('g')] * (theta[pnames.index('rad')] * Rsol) ** 2 / GG
    else:
        m1 = 10 ** theta[pnames.index('logg')] * (theta[pnames.index('rad')] * Rsol) ** 2 / GG
    derived_properties['mass'] = m1 / Msol

    if 'g2' in pnames:
        m2 = theta[pnames.index('g2')] * (theta[pnames.index('rad2')] * Rsol) ** 2 / GG
    else:
        m2 = 10 ** theta[pnames.index('logg2')] * (theta[pnames.index('rad2')] * Rsol) ** 2 / GG
    derived_properties['mass2'] = m2 / Msol

    # -- derive mass ratio
    derived_properties['q'] = m1 / m2

    # -- derive luminosity ratio
    l1 = (theta[pnames.index('rad')] ** 2 * theta[pnames.index('teff')] ** 4)
    l2 = (theta[pnames.index('rad2')] ** 2 * theta[pnames.index('teff2')] ** 4)
    derived_properties['lr'] = l1 / l2

    # -- add empty values for luminosity and distance to prevent problems with
    #   failed models
    derived_properties.update({'d': 0, 'L': 0, 'L2': 0})

    return derived_properties


def get_derived_properties_single(theta, pnames):
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

    # -- derive masses
    if 'g' in pnames:
        mass = theta[pnames.index('g')] * (theta[pnames.index('rad')] * Rsol) ** 2 / GG
    else:
        mass = 10 ** theta[pnames.index('logg')] * (theta[pnames.index('rad')] * Rsol) ** 2 / GG
    derived_properties['mass'] = mass / Msol

    # -- add empty values for luminosity and distance to prevent problems with
    #   failed models
    derived_properties.update({'d': 0, 'L': 0})

    return derived_properties


def stat_chi2(meas, e_meas, colors, syn, pars, **kwargs):
    """
    Calculate Chi2 and compute angular diameter.
   
    Colors and absolute fluxes are used to compute the Chi2, only absolute
    fluxes are used to compute angular diameter. If no absolute fluxes are
    given, the angular diameter is set to 0.

    If constraints are given in the kwargs they are included in the Chi2
    calculation. Accepted constrains are distance and q (mass ratio). Both
    should be given as a (value, error) tuple.

    @param meas: array of measurements
    @type meas: 1D array
    @param e_meas: array containing measurements errors
    @type e_meas: 1D array
    @param colors: boolean array separating colors (True) from absolute fluxes (False)
    @type colors: 1D boolean array
    @param syn: synthetic fluxes and colors
    @type syn: 1D array
    @param full_output: set to True if you want individual chisq
    @type full_output: boolean
    @return: chi-square, scale, e_scale
    @rtype: float,float,float
    """

    # check for NaNs in the observations
    nani = np.isnan(meas) | np.isnan(e_meas)
    if any(nani):
        warnings.warn('{} of the observed magnitudes are NaN values! (NaN values are ignored)'.format(np.sum(nani)))
    meas, e_meas, colors, syn = meas[~nani], e_meas[~nani], colors[~nani], syn[~nani]

    # First deal with Chi2 of the observations
    # =========================================
    # -- if syn represents only one measurement
    if sum(~colors) > 0:
        ratio = (meas / syn)[~colors]
        weights = (meas / e_meas)[~colors]
        # -- weighted average and standard deviation
        scale = np.average(ratio, weights=weights)
        e_scale = np.sqrt(np.dot(weights, (ratio - scale) ** 2) / weights.sum())
    else:
        scale, e_scale = 0, 0
    # -- we don't need to scale the colors, only the absolute fluxes
    chisq = np.where(colors, (syn - meas) ** 2 / e_meas ** 2, (syn * scale - meas) ** 2 / e_meas ** 2)

    # Then add Chi2 of derived properties as distance, mass ratio, ...
    # =================================================================
    derived_properties = kwargs.get('derived_properties', {})
    constraints = kwargs.get('constraints', {})

    for con, (c, c_m, c_p) in list(constraints.items()):
        syn_c = None
        if con == 'distance' and sum(~colors) > 0:
            # distance can only be constrained if there is at least 1 absolute measurement
            c, c_m, c_p = 1. / c ** 2, 2. * c_m / c ** 3, 2. * c_p / c ** 3
            syn_c = scale

        elif con in pars:
            # constrained on a fitted parameter
            syn_c = pars[con]

        elif con in derived_properties:
            # constrained on a derived parameter
            syn_c = derived_properties[con]

        if not syn_c is None:
            chi2_c = np.where(syn_c < c, (syn_c - c) ** 2 / c_m ** 2, (syn_c - c) ** 2 / c_p ** 2)

            # append to chisq array
            chisq = np.append(chisq, chi2_c)

    return chisq.sum(axis=0), scale, e_scale


def stat_chi2_multiple(meas, e_meas, colors, syn, pars, **kwargs):
    """
    Calculate Chi2 and compute angular diameter.

    Colors and absolute fluxes are used to compute the Chi2, only absolute
    fluxes are used to compute angular diameter. If no absolute fluxes are
    given, the angular diameter is set to 0.

    If constraints are given in the kwargs they are included in the Chi2
    calculation. Accepted constrains are distance and q (mass ratio). Both
    should be given as a (value, error) tuple.

    @param meas: array of measurements
    @type meas: 1D array
    @param e_meas: array containing measurements errors
    @type e_meas: 1D array
    @param colors: boolean array separating colors (True) from absolute fluxes (False)
    @type colors: 1D boolean array
    @param syn: synthetic fluxes and colors
    @type syn: 1D array
    @param full_output: set to True if you want individual chisq
    @type full_output: boolean
    @return: chi-square, scale, e_scale
    @rtype: float,float,float
    """

    # First deal with Chi2 of the observations
    # =========================================
    # -- if syn represents a grid of measurements
    if sum(~colors) > 0:
        ratio = (meas / syn)[~colors]
        weights = (meas / e_meas)[~colors]
        # -- weighted average and standard deviation
        scale = np.average(ratio, weights=weights.reshape(-1), axis=0)
        e_scale = np.sqrt(np.dot(weights.T, (ratio - scale) ** 2) / weights.sum(axis=0))[0]
    else:
        scale, e_scale = np.zeros(syn.shape[1]), np.zeros(syn.shape[1])

    # -- we don't need to scale the colors, only the absolute fluxes
    chisq = np.where(colors.reshape(-1, 1), (syn - meas) ** 2 / e_meas ** 2, (syn * scale - meas) ** 2 / e_meas ** 2)

    # Then add Chi2 of derived properties as distance, mass ratio, ...
    # =================================================================
    derived_properties = kwargs.get('derived_properties', {})
    constraints = kwargs.get('constraints', {})

    for con, (c, c_m, c_p) in list(constraints.items()):
        syn_c = None
        if con == 'distance' and sum(~colors) > 0:
            # distance can only be constrained if there is at least 1 absolute measurement
            c, c_m, c_p = 1. / c ** 2, 2. * c_m / c ** 3, 2. * c_p / c ** 3
            syn_c = scale

        elif con in pars:
            # constrained on a fitted parameter
            syn_c = pars[con]

        elif con in derived_properties:
            # constrained on a derived parameter
            syn_c = derived_properties[con]

        if not syn_c is None:
            chi2_c = np.where(syn_c < c, (syn_c - c) ** 2 / c_m ** 2, (syn_c - c) ** 2 / c_p ** 2)

            # append to chisq array
            chisq = np.vstack([chisq, chi2_c])

    return chisq.sum(axis=0), scale, e_scale
