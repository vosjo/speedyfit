import os
import h5py
import yaml
import datetime

import numpy as np

try:
   import gzip
   has_zip = True
except Exception, d:
   has_zip = False

from astropy.io import fits

def read2list(filename, commentchar='#', splitchar=None, skip_empty=True, skip_lines=0, **kwargs):
   """
   Load an ASCII file to list of lists.
   
   The comments and data go to two different lists.
   
   Also opens gzipped files.
   
   @param filename: name of file with the data
   @type filename: string
   @keyword commentchar: character(s) denoting comment rules
   @type commentchar: list of str
   @keyword splitchar: character seperating entries in a row (default: whitespace)
   @type splitchar: str or None
   @keyword skip_empty: skip empty lines
   @type skip_empty: bool
   @keyword skip_lines: skip nr of lines (including comment and empty lines)
   @type skip_lines: integer
   @return: list of lists (data rows)
            list of lists (comments lines without commentchar),
   @rtype: (list,list)
   """
   
   if os.path.splitext(filename)[1] == '.gz' and has_zip:
      ff = gzip.open(filename)
   elif os.path.splitext(filename)[1] == '.gz':
      print 'gzip not installed'
      return [], []
   else:
      ff = open(filename)
      
   data = []  # data
   comm = []  # comments
   
   line_nr = -1
   
   #-- fixwidth split or character split?
   if splitchar is None or isinstance(splitchar,str):
      fixwidth = False
   else:
      fixwidth = True
   
   while 1:  # might call read several times for a file
      line = ff.readline()
      if not line: break  # end of file
      line_nr += 1
      if line_nr<skip_lines:        
            continue
      
      #-- strip return character from line
      if skip_empty and line.isspace():
            continue # empty line
      
      #-- remove return characters
      line = line.replace('\n','')
      #-- when reading a comment line
      if line[0] in commentchar:
            comm.append(line[1:])
            continue # treat next line
      
      #-- when reading data, split the line
      if not fixwidth:
            data.append(line.split(splitchar))
      else:
            data.append(fw2python(line,splitchar))
   ff.close()
   
   #-- and return the contents
   return data,comm

def read2array(filename,**kwargs):
   """
   Load ASCII file to a numpy array.
   
   For a list of extra keyword arguments, see <read2list>.
   
   If you want to return a list of the columns instead of rows, just do
   
   >>> col1,col2,col3 = ascii.read2array(myfile).T
   
   @param filename: name of file with the data
   @type filename: string
   @keyword dtype: type of numpy array (default: float)
   @type dtype: numpy dtype
   @keyword return_comments: flag to return comments (default: False)
   @type return_comments: bool
   @return: data array (, list of comments)
   @rtype: ndarray (, list)
   """
   dtype = kwargs.get('dtype',np.float)
   return_comments = kwargs.get('return_comments',False)
   data,comm = read2list(filename,**kwargs)
   data = np.array(data,dtype=dtype)
   return return_comments and (data,comm) or data

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
   