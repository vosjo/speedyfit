 
import numpy as np

import  unittest

from speedyfit import mcmc

class TestMCMC(unittest.TestCase):
   
   #def setUp(self):
      
      #self.evolution_model = 'mist'
      #self.variables = ['log_L', 'log_Teff', 'log_g', 'M_H']
      #self.limits = [(0.2, 1.1), (-1.0, 0.25), (5.0, 9.0)]
      #self.obs = np.array([-0.55, 3.67, 4.50, -0.35])
      #self.obs_err = np.array([0.15, 0.05, 0.50, 0.40])
      
   def test_derived_properties(self):
      
      theta = (26000, 5.8, 0.1429, 5771, 4.438, 1, 0.02)
      pnames = ('teff1', 'logg1', 'rad1', 'teff2', 'logg2', 'rad2', 'ebv')
      
      derived_properties = mcmc.get_derived_properties(theta, pnames)
      
      self.assertAlmostEqual(derived_properties['mass1'], 0.47, places=2, 
                           msg="sdB mass wrongly calculated: {}".format(derived_properties['mass1']))
      
      self.assertAlmostEqual(derived_properties['mass2'], 1.00, places=2, 
                           msg="Solar mass wrongly calculated: {}".format(derived_properties['mass2']))
      
      self.assertAlmostEqual(derived_properties['q'], 0.47, places=2, 
                           msg="q wrongly calculated: {}".format(derived_properties['q']))
      
   
   def test_prior(self):
      
      theta = (26000, 5.8, 0.1429, 5771, 4.438, 1, 0.02)
      pnames = ('teff1', 'logg1', 'rad1', 'teff2', 'logg2', 'rad2', 'ebv')
      
      limits = np.array([(20000, 40000), (5.0, 6.5), (0.05, 0.25), 
                         (4000, 8000),   (4.0, 5.0), (0.50, 2.10),
                         (0.00, 0.05)])
      
      derived_limits = {'q':(0.1, 1.0), 'mass1':(0.40, 0.55), 'mass2':(0.70, 2.0)}
      
      derived_properties = mcmc.get_derived_properties(theta, pnames)
      
      prior = mcmc.lnprior(theta, derived_properties, limits, 
                           pnames = pnames, derived_limits=derived_limits)
      
      self.assertEqual(prior, 0,
                       msg="Prior should have been 0 for accepted limits, was {}".format(prior))
      
      
      
      theta = (15000, 5.8, 0.1429, 5771, 4.438, 1, 0.02)
      derived_properties = mcmc.get_derived_properties(theta, pnames)
      
      prior = mcmc.lnprior(theta, derived_properties, limits, 
                           pnames = pnames, derived_limits=derived_limits)
      
      self.assertEqual(prior, -np.inf,
                       msg="theta out of limits, expected prior = -inf, was {}".format(prior))
      
      
      
      theta = (26000, 5.8, 0.10, 5771, 4.438, 1, 0.02)
      derived_properties = mcmc.get_derived_properties(theta, pnames)
      
      prior = mcmc.lnprior(theta, derived_properties, limits, 
                           pnames = pnames, derived_limits=derived_limits)
      
      self.assertEqual(prior, -np.inf,
                       msg="mass1 out of limits, expected prior = -inf, was {}".format(prior))