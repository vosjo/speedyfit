 
import numpy as np

import  unittest

from speedyfit.speedyfit import mcmc, statfunc


class TestMCMC(unittest.TestCase):
   
   #def setUp(self):
      
      #self.evolution_model = 'mist'
      #self.variables = ['log_L', 'log_Teff', 'log_g', 'M_H']
      #self.limits = [(0.2, 1.1), (-1.0, 0.25), (5.0, 9.0)]
      #self.obs = np.array([-0.55, 3.67, 4.50, -0.35])
      #self.obs_err = np.array([0.15, 0.05, 0.50, 0.40])
      
   
   def test_prior(self):
      
      theta = (26000, 5.8, 0.1429, 5771, 4.438, 1, 0.02)
      pnames = ('teff', 'logg', 'rad', 'teff2', 'logg2', 'rad2', 'ebv')
      
      limits = np.array([(20000, 40000), (5.0, 6.5), (0.05, 0.25), 
                         (4000, 8000),   (4.0, 5.0), (0.50, 2.10),
                         (0.00, 0.05)])
      
      derived_limits = {'q':(0.1, 1.0), 'mass':(0.40, 0.55), 'mass2':(0.70, 2.0)}
      
      derived_properties = statfunc.get_derived_properties(theta, pnames)
      
      prior = mcmc.lnprior(theta, derived_properties, limits,
                           pnames = pnames, derived_limits=derived_limits)
      
      self.assertEqual(prior, 0,
                       msg="Prior should have been 0 for accepted limits, was {}".format(prior))
      
      
      
      theta = (15000, 5.8, 0.1429, 5771, 4.438, 1, 0.02)
      derived_properties = statfunc.get_derived_properties(theta, pnames)
      
      prior = mcmc.lnprior(theta, derived_properties, limits,
                           pnames = pnames, derived_limits=derived_limits)
      
      self.assertEqual(prior, -np.inf,
                       msg="theta out of limits, expected prior = -inf, was {}".format(prior))
      
      
      
      theta = (26000, 5.8, 0.10, 5771, 4.438, 1, 0.02)
      derived_properties = statfunc.get_derived_properties(theta, pnames)
      
      prior = mcmc.lnprior(theta, derived_properties, limits,
                           pnames = pnames, derived_limits=derived_limits)
      
      self.assertEqual(prior, -np.inf,
                       msg="mass out of limits, expected prior = -inf, was {}".format(prior))
      
      
   #def test_likelihood(self):
      
      #photbands = ['STROMGREN.U', 'STROMGREN.B']
      #gridfilename = '/home/joris/Python/ivsdata/sedtables/modelgrids/ikurucz93_z0.0_k2odfnew_sed_lawfitzpatrick2004_Rv3.10.fits'
      
      #axis_values, grid_pars, pixelgrid, grid_names = model.prepare_grid(photbands, gridfilename,
                 #teffrange=(4000, 5000),loggrange=(3.5, 4.0),
                 #ebvrange=(0.0, 0.02),
                 #variables=['teff','logg','ebv'])
      
      #grids = [[axis_values, pixelgrid],]
      
      
      #theta = (4700, 3.8, 0.01)
      #pnames = ['teff', 'logg', 'ebv']
      #derived_properties = statfunc.get_derived_properties(theta, pnames)
      
      #y = np.array([5.0, 6.0])
      #yerr = np.array([0.5, 0.75])
      
      #kwargs = dict(pnames = pnames,
                    #photbands = photbands, 
                    #grid=grids)
      
      #print mcmc.lnlike(theta, derived_properties, y, yerr, **kwargs)
      
      
      #self.assertTrue(True)
      