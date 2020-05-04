 
import numpy as np

import  unittest

from speedyfit.speedyfit import statfunc

class TestMCMC(unittest.TestCase):
      
   def test_derived_properties(self):
      
      theta = (26000, 5.8, 0.1429, 5771, 4.438, 1, 0.02)
      pnames = ('teff', 'logg', 'rad', 'teff2', 'logg2', 'rad2', 'ebv')
      
      derived_properties = statfunc.get_derived_properties(theta, pnames)
      
      self.assertAlmostEqual(derived_properties['mass'], 0.47, places=2, 
                           msg="sdB mass wrongly calculated: {}".format(derived_properties['mass']))
      
      self.assertAlmostEqual(derived_properties['mass2'], 1.00, places=2, 
                           msg="Solar mass wrongly calculated: {}".format(derived_properties['mass2']))
      
      self.assertAlmostEqual(derived_properties['q'], 0.47, places=2, 
                           msg="q wrongly calculated: {}".format(derived_properties['q']))
      
      
   def test_stat_chi2(self):
      
      meas = np.array([1.0,2.0,3.0])
      e_meas = np.array([0.25, 0.30, 0.15])
      syn = np.array([1.5, 1.89, 3.6])
      
      #-- colors (none scaled)
      colors = np.array([True,True,True])
      chi2, scale, e_scale = statfunc.stat_chi2(meas, e_meas, colors, syn)
      
      self.assertAlmostEqual(chi2, 20.1344, 2,
                             msg="Chi2 for colours not correct ({} != {})".format(chi2, 20.1344))
      
      #-- absolute magnitudes (scaled)
      colors = np.array([False,False,False])
      chi2, scale, e_scale = statfunc.stat_chi2(meas, e_meas, colors, syn)
      
      self.assertAlmostEqual(chi2, 3.32834, 2,
                             msg="Chi2 for colours not correct ({} != {})".format(chi2, 3.3283))
      
      #-- added distance
      constraints = {'distance': (1.5, 0.15)}
      colors = np.array([False,False,False])
      chi2, scale, e_scale = statfunc.stat_chi2(meas, e_meas, colors, syn, **constraints)
      
      self.assertAlmostEqual(chi2, 25.2343, 2,
                             msg="Chi2 for distance not correct ({} != {})".format(chi2, 25.2343))
      
      #self.assertTrue(False)