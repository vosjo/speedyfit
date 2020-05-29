
import os

import numpy as np
 
from astropy.io import ascii
basedir = os.path.dirname(__file__)

def is_color(photband):
   """
   Return true if the photometric passband is actually a color.
   
   @param photband: name of the photometric passband
   @type photband: string
   @return: True or False
   @rtype: bool
   """
   if '-' in photband.split('.')[1]:
      return True
   elif photband.split('.')[1].upper() in ['M1','C1']:
      return True
   else:
      return False


def eff_wave(photband):
   """
   Returns the effective wavelength of the pass band in angstrom
   
   @param photband: name of the photometric passband
   @type photband: string
   @return: effective wavelength
   @rtype: float
   """
   
   data = ascii.read(os.path.join(basedir, 'zeropoints.dat'), comment="\s*#")
   
   s = np.where(data['photband'] == photband)
   
   return data['eff_wave'][s][0]


def mag2flux(mag, error, photband):
   """
   Converts a magnitude in a given photband to a flux
   
   Flux is returned in units of erg/s/cm2/AA
   """
   data = ascii.read(os.path.join(basedir, 'zeropoints.dat'), comment="\s*#")
   
   if not isinstance(photband, str):
      flux, err = np.zeros_like(mag), np.zeros_like(mag)
      for i, (m, e, b) in enumerate(zip(mag, error, photband)):
         
         s = np.where(data['photband'] == b)
         
         F0 = data['Flam0'][s][0]
         zpcor = data['zp_corr'][s][0]
         
         f = 10**(-(m-zpcor)/2.5)*F0
         fe = np.log(10) * e / 2.5 * f
         
         flux[i] = f
         err[i] = fe
   
   else:
      s = np.where(data['photband'] == photband)
      data = data[s]

      F0 = data['Flam0'][0]
      zpcor = data['zp_corr'][0]
      
      flux = 10**(-(mag-zpcor)/2.5)*F0
      err = np.log(10) * error / 2.5 * flux

      if hasattr(err, '__iter__'):
         err = np.where(err == 0, flux / 10, err)
      elif err == 0:
         err = flux / 10
   
   return flux, err


def flux2mag(flux, error, photband):
   """
   Converts a flux in a given photband to a magnitude
   
   Flux has to be provided in units of erg/s/cm2/AA
   """
   data = ascii.read(os.path.join(basedir, 'zeropoints.dat'), comment="\s*#")

   if not isinstance(photband, str):
      mag, err = np.zeros_like(flux), np.zeros_like(flux)
      for i, (f, e, b) in enumerate(zip(flux, error, photband)):
         
         s = np.where(data['photband'] == b)
         
         F0 = data['Flam0'][s][0]
         zpcor = data['zp_corr'][s][0]
         
         m = -2.5*np.log10(f/F0)+zpcor
         me = 2.5 / np.log(10) * e / f
         
         mag[i] = m
         err[i] = me
   
   else:
      s = np.where(data['photband'] == photband)
      data = data[s]
      
      F0 = data['Flam0'][0]
      zpcor = data['zp_corr'][0]
      
      mag = -2.5*np.log10(flux/F0)+zpcor
      err = 2.5 / np.log(10) * error / flux
   
   return  mag, err


if __name__=="__main__":
   
   print((mag2flux(13.0, 0.01, '2MASS.H')))
