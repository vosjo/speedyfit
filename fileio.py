import os
import h5py
import yaml
import datetime

from astropy.io import fits

def write2hdf5(samples, blobs, filename):
   """
   
   """
   
   hdf = h5py.File(filename)
   
   hdf.create_dataset('samples', data=samples)
   hdf.create_dataset('blobs', data=blobs)
   
   hdf.close()
   
   
def read_hdf5(filename):
   
   hdf = h5py.File(filename, 'r')
   
   samples = hdf['samples'].value
   blobs = hdf['blobs'].value
   
   hdf.close()
   
   return samples, blobs
   
   
def write2fits(samples, filename, setup=None):
   """
   writes recarray containing all samples to fits file
   """
   
   #-- create Table data for fits file
   cols = []
   
   for (name, fmt) in samples.dtype.descr:
      col = fits.Column(name=name, format='D', array=samples[name])
      cols.append(col)
      
   cols = fits.ColDefs(cols)
   
   tbhdu = fits.BinTableHDU.from_columns(cols)
   
   #-- add date to header
   header = tbhdu.header
   header['date'] = str(datetime.datetime.now())
    
   #-- add complete setup as comment if given as str
   for line in setup.split('\n'):
      header['comment'] = line
   
   #-- delete file if it already exists
   if os.path.isfile(filename):
      os.remove(filename)
   
   tbhdu.writeto(filename)
   
def read_fits(filename):
   """
   Reads fits file and returns samples and the settings if they exist
   """
   
   hdu = fits.open(filename)
   
   samples = hdu[1].data
   
   setup = hdu[1].header['comment']
   
   return samples, str(setup)
   