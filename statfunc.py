# -*- coding: utf-8 -*-
"""
Fit SED models to observed data using various approaches.
"""
import logging
import sys
import pyfits
import itertools
import re
import copy

import pylab as pl
import numpy as np
from numpy import inf
import scipy
from scipy.interpolate import Rbf
from scipy.optimize import fmin,fmin_powell
from ivs.statistics import pca
from ivs.sed import model
from ivs.sed import filters
from ivs.sed.decorators import iterate_gridsearch,parallel_gridsearch
from ivs.sigproc import fit as sfit
from ivs.aux import numpy_ext
from ivs.aux import progressMeter
from ivs.aux.decorators import make_parallel
from ivs.units import constants
from ivs.units import conversions as cv

logger = logging.getLogger("SED.FIT")

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
   
   #-- derive luminosity ratio
   if 'rad' in pnames and 'teff' in pnames and 'rad2' in pnames and 'teff2' in pnames:
      l1 = ( theta[pnames.index('rad')]**2 * theta[pnames.index('teff')]**4 )
      l2 = ( theta[pnames.index('rad2')]**2 * theta[pnames.index('teff2')]**4 )
      derived_properties['lr'] = l1 / l2
      
   if 'rad1' in pnames and 'teff1' in pnames and 'rad2' in pnames and 'teff2' in pnames:
      l1 = ( theta[pnames.index('rad1')]**2 * theta[pnames.index('teff1')]**4 )
      l2 = ( theta[pnames.index('rad2')]**2 * theta[pnames.index('teff2')]**4 )
      derived_properties['lr'] = l1 / l2
   
   return derived_properties



def stat_chi2(meas, e_meas, colors, syn, **kwargs):
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
   #=========================================
   #-- if syn represents only one measurement
   if sum(-colors) > 0:
      ratio = (meas/syn)[-colors]
      weights = (meas/e_meas)[-colors]
      #-- weighted average and standard deviation
      scale = np.average(ratio,weights=weights)
      e_scale = np.sqrt(np.dot(weights, (ratio-scale)**2)/weights.sum())
   else:
      scale,e_scale = 0,0
      
   #-- we don't need to scale the colors, only the absolute fluxes
   chisq = np.where(colors, (syn-meas)**2/e_meas**2, (syn*scale-meas)**2/e_meas**2)
   
   
   # Then add Chi2 of derived properties as distance, mass ratio, ...
   #=================================================================
   derived_properties = kwargs.get('constraints_syn', {})
   
   #-- Chi2 of the distance measurement
   #   this can only be done if absolute fluxes are included in the fit.
   if 'distance' in kwargs and sum(-colors) > 0:
      syn_scale = 1/kwargs['distance'][0]**2
      syn_scale_e = 2. / kwargs['distance'][0]**3 * kwargs['distance'][1]
      
      chi2_d = (scale - syn_scale)**2 / syn_scale_e**2
      
      # append to chisq array
      chisq = np.append(chisq, chi2_d)
      
      ##-- add distance to derived properties
      #derived_properties['distance'] = 1 / np.sqrt(scale)
   
   #-- Chi2 of the Mass-Ratio
   #   q needs to be computed elsewhere and provided in the kwargs
   
   if 'q' in kwargs and 'q' in derived_properties:
      q, q_e = kwargs['q'][0], kwargs['q'][1]
      syn_q = derived_properties['q']
      
      chi2_q = (q - syn_q)**2 / q_e**2
      
      # append to chisq array
      chisq = np.append(chisq, chi2_q)
   
   
   return chisq.sum(axis=0), scale, e_scale
  

