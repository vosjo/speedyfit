 
import numpy as np

import  unittest

from speedyfit import mcmc, model

class TestModel(unittest.TestCase):
   
   #def setUp(self):
      
      #self.evolution_model = 'mist'
      #self.variables = ['log_L', 'log_Teff', 'log_g', 'M_H']
      #self.limits = [(0.2, 1.1), (-1.0, 0.25), (5.0, 9.0)]
      #self.obs = np.array([-0.55, 3.67, 4.50, -0.35])
      #self.obs_err = np.array([0.15, 0.05, 0.50, 0.40])
      
   def test_itable_single(self):
      
      photbands = ['STROMGREN.U', 'STROMGREN.B']
      gridfilename = '/home/joris/Python/ivsdata/sedtables/modelgrids/ikurucz93_z0.0_k2odfnew_sed_lawfitzpatrick2004_Rv3.10.fits'
      
      axis_values, grid_pars, pixelgrid, grid_names = model.prepare_grid(photbands, gridfilename,
                 teffrange=(4000, 5000),loggrange=(3.5, 4.0),
                 ebvrange=(0.0, 0.02),
                 variables=['teff','logg','ebv'])
      
      grid = [axis_values, pixelgrid]
      
      #-- check if interpolated model is correct
      flux1, Labs1 = model.get_itable_single(teff=4750, logg=3.75, ebv=0.01, grid=grid)
      
      self.assertAlmostEqual(Labs1, 0.45605564, 3,
                             msg='Basic model: Luminosity is not correct ({} != 0.456)'.format(Labs1))
      
      #-- check if radius is treated correctly
      flux1, Labs1 = model.get_itable_single(teff=4750, logg=3.75, ebv=0.01, rad=1, grid=grid)
      self.assertAlmostEqual(Labs1, 0.45605564, 3,
                             msg='Solar Radius: Luminosity is not correct ({} != 0.456)'.format(Labs1))
      
      flux1, Labs1 = model.get_itable_single(teff=4750, logg=3.75, ebv=0.01, rad=2, grid=grid)
      self.assertAlmostEqual(Labs1, 2**2*0.45605564, 3,
                             msg='Double Radius: Luminosity is not correct ({} != 1.824)'.format(Labs1))
      
      #-- compare with get_itable
      flux1, Labs1 = model.get_itable_single(teff=4750, logg=3.75, ebv=0.01, grid=grid)
      flux2, Labs2 = model.get_itable(teff=4750, logg=3.75, ebv=0.01, grid=grid)
      
      
      for f1, f2 in zip(flux1, flux2):
         self.assertEqual(f1, f2, 
                          msg="output of get_itable should be same as get_itable_single, ({} != {})".format(flux1, flux2))
      self.assertEqual(Labs1, Labs2, 
                          msg="output of get_itable should be same as get_itable_single, ({} != {})".format(flux1, flux2))
      
      ##print axis_values
      ##print pixelgrid
      #print pixelgrid.shape
      
      #self.assertTrue(False)
      
   def test_itable(self):
      
      photbands = ['STROMGREN.U', 'STROMGREN.B']
      gridfilename = '/home/joris/Python/ivsdata/sedtables/modelgrids/ikurucz93_z0.0_k2odfnew_sed_lawfitzpatrick2004_Rv3.10.fits'
      
      
      axis_values, grid_pars, pixelgrid, grid_names = model.prepare_grid(photbands, gridfilename,
                 teffrange=(4000, 5000),loggrange=(3.5, 4.0),
                 ebvrange=(0.0, 0.02),
                 variables=['teff','logg','ebv'])
      
      grid1 = [axis_values, pixelgrid]
      
      gridfilename = '/home/joris/Python/ivsdata/sedtables/modelgrids/iTMAP2012_sdB_extended_lawfitzpatrick2004_Rv3.10.fits'
      
      axis_values, grid_pars, pixelgrid, grid_names = model.prepare_grid(photbands, gridfilename,
                 teffrange=(20000, 40000),loggrange=(5.0, 6.5),
                 ebvrange=(0.0, 0.02),
                 variables=['teff','logg','ebv'])
      
      grid2 = [axis_values, pixelgrid]
      
      flux1, Labs1 = model.get_itable(teff=4750, logg=3.75, ebv=0.01, rad=1, grid=grid1)
      flux2, Labs2 = model.get_itable(teff=4250, logg=3.50, ebv=0.01, rad=1, grid=grid1)
      
      fluxt, Labst = model.get_itable(teff=4750, logg=3.75, ebv=0.01, rad=1, 
                                      teff2=4250, logg2=3.50, ebv2=0.01, rad2=1, grid=[grid1, grid1])
      
      for f1, f2, ft in zip(flux1, flux2, fluxt):
         self.assertEqual(f1+f2, ft, 
                          msg="binary SED should have sum of flux of single SED, ({} + {} != {})".format(flux1, flux2, fluxt))
      
      
      flux1, Labs1 = model.get_itable(teff=4800, logg=3.8, ebv=0.01, rad=2, grid=grid1)
      flux2, Labs2 = model.get_itable(teff=25000, logg=5.8, ebv=0.01, rad=0.14, grid=grid2)
      
      fluxt, Labst = model.get_itable(teff=4800, logg=3.8, ebv=0.01, rad=2, 
                                      teff2=25000, logg2=5.8, ebv2=0.01, rad2=0.14, 
                                      grid=[grid1, grid2])
      
      for f1, f2, ft in zip(flux1, flux2, fluxt):
         self.assertEqual(f1+f2, ft, 
                          msg="binary SED should have sum of flux of single SED, ({} + {} != {})".format(flux1, flux2, fluxt))
      