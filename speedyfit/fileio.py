import os
import h5py
import datetime

import numpy as np

try:
   import gzip
   has_zip = True
except Exception as d:
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
      print('gzip not installed')
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

def write_array(data, filename, **kwargs):
    """
    Save a numpy array to an ASCII file.
    
    Add comments via keyword comments (a list of strings denoting every comment
    line). By default, the comment lines will be preceded by the C{commentchar}.
    If you want to override this behaviour, set C{commentchar=''}.
    
    If you give a record array, you can simply set C{header} to C{True} to write
    the header, instead of specifying a list of strings.
    
    @keyword header: optional header for column names
    @type header: list of str (or boolean for record arrays)
    @keyword comments: comment lines
    @type comments: list of str
    @keyword commentchar: comment character
    @type commentchar: str
    @keyword sep: separator for the columns and header names
    @type sep: str
    @keyword axis0: string denoting the orientation of the matrix. If you gave
    a list of columns, set C{axis0='cols'}, otherwise C{axis='rows'} (default).
    @type axis0: str, one of C{cols}, C{rows}.
    @keyword mode: file mode (a for appending, w for (over)writing...)
    @type mode: char (one of 'a','w'...)
    @keyword auto_width: automatically determine the width of the columns
    @type auto_width: bool
    @keyword formats: formats to use to write each column
    @type formats: list of string formatters
    """
    header = kwargs.get('header',None)
    comments = kwargs.get('comments',None)
    commentchar = kwargs.get('commentchar','#')
    sep = kwargs.get('sep',' ')
    axis0 = kwargs.get('axis0','rows')
    mode = kwargs.get('mode','w')
    auto_width = kwargs.get('auto_width',False)
    formats = kwargs.get('formats',None)
    # use '%g' or '%f' or '%e' for writing floats automatically from record arrays with auto width
    use_float = kwargs.get('use_float','%f') 
    
    #-- switch to rows first if a list of columns is given
    if not isinstance(data,np.ndarray):
        data = np.array(data)
    if not 'row' in axis0.lower():
        data = data.T
    
    if formats is None:
        try:
            formats = [('S' in str(data[col].dtype) and '%s' or use_float) for col in data.dtype.names]
        except TypeError:
            formats = [('S' in str(col.dtype) and '%s' or '%s') for col in data.T]
    #-- determine width of columns: also take the header label into account
    col_widths = []
    #-- for record arrays
    if auto_width is True and header==True:
        for fmt,head in zip(formats,data.dtype.names):
            col_widths.append(max([len('%s'%(fmt)%(el)) for el in data[head]]+[len(head)]))
    #-- for normal arrays and specified header
    elif auto_width is True and header is not None:
        for i,head in enumerate(header):
            col_widths.append(max([len('%s'%(formats[i])%(el)) for el in data[:,i]]+[len(head)]))
    #-- for normal arrays without header
    elif auto_width is True and header is not None:
        for i in range(data.shape[1]):
            col_widths.append(max([len('%s'%(formats[i])%(el)) for el in data[:,i]]))
    
    if header is True:
        col_fmts = [str(data.dtype[i]) for i in range(len(data.dtype))]
        header = data.dtype.names
    else:
        col_fmts = None
    
    ff = open(filename,mode)
    if comments is not None:
        ff.write('\n'.join(comments)+'\n')
    
    #-- WRITE HEADER
    #-- when header is desired and automatic width
    if header is not None and col_widths:
        ff.write('#'+sep.join(['%%%s%ss'%(('s' in fmt and '-' or ''),cw)%(head) for fmt,head,cw in zip(formats,header,col_widths)])+'\n')
    #-- when header is desired
    elif header is not None:
        ff.write('#'+sep.join(header)+'\n')
    
    #-- WRITE COLUMN FORMATS
    if col_fmts is not None and col_widths:
        ff.write('#'+sep.join(['%%%s%ss'%(('s' in fmt and '-' or ''),cw)%(colfmt) for fmt,colfmt,cw in zip(formats,col_fmts,col_widths)])+'\n')
    elif col_fmts is not None:
        ff.write('#'+sep.join(['%%%ss'%('s' in fmt and '-' or '')%(colfmt) for fmt,colfmt in zip(formats,col_fmts)])+'\n')
    
    #-- WRITE DATA
    #-- with automatic width
    if col_widths:
        for row in data:
            ff.write(' '+sep.join(['%%%s%s%s'%(('s' in fmt and '-' or ''),cw,fmt[1:])%(col) for col,cw,fmt in zip(row,col_widths,formats)])+'\n')
    #-- without automatic width
    else:
        for row in data:
            ff.write(sep.join(['%s'%(col) for col in row])+'\n')
    ff.close()

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
   

def write_summary2hdf5(objectname, samples, obs, obs_err, photbands, pars={}, grids=[], filename=None):
   """
   Writes the obeservations and results to hdf5 format readable by AOTS
   """
   from astropy.table import Table
   from . import model
   from . import filters
   from .photometry_query import get_coordinate
   
   if filename is None: filename = objectname + '_sedfit.hdf5'
   f = h5py.File(filename, "w")
   
   ra, dec = get_coordinate(objectname)
      
   f.attrs['systemname'] = objectname
   f.attrs['ra'] = ra
   f.attrs['dec'] = dec
   
   f.attrs['name'] = 'SED fit'
   f.attrs['note'] = ''
   f.attrs['type'] = 'sedfit'
   
   pars = pars.copy()
   
   # use model from 'best' results or 'pc' results
   result = 'best'
   resi = 0 if result == 'best' else 1
   
   
   colors = np.array([filters.is_color(p) for p in photbands])
   waves = np.array([filters.eff_wave(p) for p in photbands[~colors]])
   pbands = np.array(photbands[~colors], dtype='<S') # needed to convert to S type for compatibility with h5 format.
   obsdata = Table([waves, pbands, obs[~colors], obs_err[~colors]], names=['wave', 'photband', 'flux', 'flux_err'])
   
   #-- create the data group
   data = f.create_group('DATA')
   data.attrs['xlabel'] = 'Wavelength (AA)'
   data.attrs['ylabel'] = 'flux'
   data.attrs['xmin'] = np.min(waves)
   data.attrs['xmax'] = np.max(waves)
   data.attrs['xscale'] = 'log'
   data.attrs['yscale'] = 'log'
   
   obsgr = data.create_dataset('Obs', data=obsdata)
   obsgr.attrs['label'] = 'OBS'
   obsgr.attrs['datatype'] = 'discrete'
   obsgr.attrs['xpar'] = 'wave'
   obsgr.attrs['ypar'] = 'flux'
   
   
   def get_unit(par):
      
      if 'teff' in par: return 'K'
      if 'logg' in par: return 'dex'
      if 'L' in par: return 'Lsol'
      if 'rad' in par: return 'Rsol'
      if 'd' in par: return 'pc'
      if 'ebv' in par: return ''
   
      return ''
   
   if not len(list(pars.keys())) == 0:
      
      
      #-- create parameters
      group = f.create_group('PARAMETERS')
      ipars = pars.copy()
      for key, value in list(pars.items()):
         
         ipars[key] = [value[resi]]
         pars[key] = value[resi]
         
         # if key == 'scale' or key == 'chi2' or key == 'mass' or key == 'logg': continue
         if key in ['scale', 'chi2', 'mass', 'logg']: continue
         if key not in samples.dtype.names: continue
         
         pg = group.create_group(key)
         pg.attrs['unit'] = get_unit(key)
         pg.attrs['value'] = value[resi]
         pg.attrs['err'] = value[2]
         if key == 'g' or key == 'logg':
            pg.attrs['valid'] = False
         else:
            pg.attrs['valid'] = True
         
         y, x = np.histogram(samples[key], 20, density=False)
         x = (x[0:-1] + x[1:])/2.
         y = y / np.float(np.max(y))
         
         dis = Table([x, y], names=['x', 'y'])
         
         
         disgr = pg.create_dataset('DISTRIBUTION', data=dis)
         disgr.attrs['datatype'] = 'bars'
         disgr.attrs['xpar'] = 'x'
         disgr.attrs['ypar'] = 'y'
         
      _ = ipars.pop('d')
      _ = pars.pop('d')
      
      
       #-- create the model group
      group = f.create_group('MODEL')
      group.attrs['xlabel'] = 'Wavelength (AA)'
      group.attrs['ylabel'] = 'Flux'
      #group.attrs['xmin'] = 0.
      #group.attrs['xmax'] = 1.
      
      
      syn, Labs = model.get_itable(grid=grids, photbands=photbands, **ipars)
      syn = syn[:,0]
      
      scale = pars['scale']
      
      
      itable = Table([waves, scale*syn], names=['wave', 'flux'])
      
      itb = group.create_dataset('Iflux', data=itable)
      itb.attrs['label'] = 'Synth. photometry'
      itb.attrs['datatype'] = 'discrete'
      itb.attrs['xpar'] = 'wave'
      itb.attrs['ypar'] = 'flux'
      
      
      #-- if possible, plot a non integrated model
      #   can only be done if provided gridnames are not paths to integrated files.
      if len(grids) > 0 and not os.path.isfile(grids[0]):
         #-- synthetic model
         wave, flux = model.get_table(grid=grids, **pars)
         s = np.where((wave > np.min(waves)-1000) & (wave < np.max(waves) + 1000))
         table = Table([wave[s], scale*flux[s]], names=['wave', 'flux'])
         
         name = 'binary' if 'teff2' in pars else grids[0]
         tb = group.create_dataset(name, data=table)
         tb.attrs['label'] = name
         tb.attrs['datatype'] = 'continuous'
         tb.attrs['xpar'] = 'wave'
         tb.attrs['ypar'] = 'flux'
         
         
         #-- plot components
         if 'teff2' in pars:
            wave, flux = model.get_table(grid=grids[0], teff=pars['teff'], logg=pars['logg'], rad=pars['rad'], ebv=pars['ebv'])
            s = np.where((wave > np.min(waves)-1000) & (wave < np.max(waves) + 1000))
            table1 = Table([wave[s], scale*flux[s]], names=['wave', 'flux'])
            
            tb1 = group.create_dataset(grids[0], data=table1)
            tb1.attrs['label'] = grids[0]
            tb1.attrs['datatype'] = 'continuous'
            tb1.attrs['xpar'] = 'wave'
            tb1.attrs['ypar'] = 'flux'
            
            wave, flux = model.get_table(grid=grids[1], teff=pars['teff2'], logg=pars['logg2'], rad=pars['rad2'], ebv=pars['ebv'])
            s = np.where((wave > np.min(waves)-1000) & (wave < np.max(waves) + 1000))
            table2 = Table([wave[s], scale*flux[s]], names=['wave', 'flux'])
            
            tb2 = group.create_dataset(grids[1], data=table2)
            tb2.attrs['label'] = grids[1]
            tb2.attrs['datatype'] = 'continuous'
            tb2.attrs['xpar'] = 'wave'
            tb2.attrs['ypar'] = 'flux'
            
            
         #-- create the O-C data
         group = f.create_group('O-C')
         group.attrs['xlabel'] = 'Wavelength (AA)'
         group.attrs['ylabel'] = 'O-C (mag)'
         group.attrs['xscale'] = 'log'
         group.attrs['xmin'] = np.min(waves)-1000
         group.attrs['xmax'] = np.max(waves)+1000
         
         
         y = 2.5*np.log10(syn*scale / obs[~colors])
         yerr = 2.5 / np.log(10.0) * obs_err[~colors] / obs[~colors]
         obsdata = Table([waves, y, yerr], names=['wave', 'o-c', 'o-c_err'])
         
         octb = group.create_dataset('Obs', data=obsdata)
         octb.attrs['label'] = 'Obs - Calc'
         octb.attrs['datatype'] = 'discrete'
         octb.attrs['xpar'] = 'wave'
         octb.attrs['ypar'] = 'o-c'
               
   f.close()
