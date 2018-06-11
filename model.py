import re
import numpy as np
 
from astropy.io import fits 

import interpol

from ivs.aux import loggers
logger = loggers.get_basic_logger()

def load_grids(gridnames, pnames, limits):
   """
   prepares the integrated photometry grid by loading the grid and cutting it to the size
   given in limits.
   """
   grids = []
   for i, name in enumerate(gridnames):
      
      ind = '' if i == 0 else str(i + 1)
      
      axis_values, grid_pars, pixelgrid, grid_names = model.prepare_grid(photbands, name,
            teffrange=limits[pnames.index('teff'+ind)],
            loggrange=limits[pnames.index('logg'+ind)],
            ebvrange =limits[pnames.index('ebv')],
            variables=['teff','logg','ebv'])

      grids.append([axis_values, pixelgrid])
      
   return grids

def prepare_grid(photbands, gridfilename,
                 teffrange=(-np.inf,np.inf),loggrange=(-np.inf,np.inf),
                 ebvrange=(-np.inf,np.inf),
                 variables=['teff','logg','ebv'],**kwargs):
   
   flux = []
   grid_pars = []
   grid_names = np.array(variables)
   
   with fits.open(gridfilename) as ff:
      #-- make an alias for further reference
      ext = ff[1]
      
      #-- the grid is already cut here to limit memory usage
      keep = np.ones(len(ext.data),bool)
      for name in variables:
         #-- first find the closest actuall grid points
         low,high = locals()[name+'range']
         lidx = np.abs(ext.data.field(name)[ext.data.field(name) <= low]-low).argmin()
         hidx = np.abs(ext.data.field(name)[ext.data.field(name) >= high]-high).argmin()
         low = ext.data.field(name)[ext.data.field(name) <= low][lidx]
         high = ext.data.field(name)[ext.data.field(name) >= high][hidx]
         
         #-- we need to be carefull for rounding errors
         in_range = (low<=ext.data.field(name)) & (ext.data.field(name)<=high)
         on_edge  = np.allclose(ext.data.field(name),low) | np.allclose(ext.data.field(name),high)
         
         keep = keep & (in_range | on_edge)
         
      grid_pars = np.vstack([ext.data.field(name)[keep] for name in variables])
      flux = _get_flux_from_table(ext,photbands,include_Labs=True)[keep]
      
   flux = np.log10(flux)
   
   #-- create the pixeltype grid
   axis_values, pixelgrid = interpol.create_pixeltypegrid(grid_pars,flux.T)
   return axis_values,grid_pars.T,pixelgrid,grid_names

def get_itable_single(teff=None, logg=None, g=None, ebv=0.0, **kwargs):
   
   #-- get the grid
   axis_values, pixelgrid = kwargs['grid']
   
   #-- check if logg or g is given
   if logg is None and not g is None:
      logg = np.log10(g)
   
   multiple = False
   if not hasattr(teff, '__iter__') or not hasattr(logg, '__iter__') or \
      not hasattr(ebv, '__iter__'):
      teff = np.array(teff)
      logg = np.array(logg)
      ebv = np.array(ebv)
      multiple = True
   
   p = np.vstack([teff, logg, ebv])
   
   values = interpol.interpolate(p, axis_values, pixelgrid)
      
   #-- switch logarithm to normal
   values = 10**values
   flux,Labs = values[:-1],values[-1]
   
   #-- Take radius into account when provided
   if 'rad' in kwargs:
      rad = np.array(kwargs['rad'])
      flux,Labs = flux*rad**2, Labs*rad**2
      
   if multiple:
      flux = flux.flatten()
      Labs = Labs.flatten()
   
   return flux, Labs
   

def get_itable(grid=[], **kwargs):
   
   values, parameters, components = {}, set(), set()
   for key in kwargs.keys():
      if re.search("^(teff|logg|g|ebv|rad)\d?$", key):
         par, comp = re.findall("^(teff|logg|g|ebv|rad)(\d?)$", key)[0]
         values[key] = kwargs.pop(key)
         parameters.add(par)
         components.add(comp)
   
   # need to check here that grid is same length as components
   if hasattr(grid, '__iter__') and  hasattr(grid[0], '__iter__') and len(grid[0]) == 2:
      grids = grid
   else:
      grids = [grid]
   
   #-- If there is only one component, we can directly return the result
   if len(components) == 1:
      kwargs.update(values)
      fluxes, Labs = get_itable_single(grid=grids[0], **kwargs)
      return fluxes, {'L':Labs}
   
   fluxes, Labs = [],{}                               
   for i, (comp, grid) in enumerate(zip(components,grids)):
      kwargs_ = kwargs.copy()
      for par in parameters:
         kwargs_[par] = values[par+comp] if par+comp in values else values[par]
      
      f,L = get_itable_single(grid=grid, **kwargs_)
                                    
      fluxes.append(f)
      Labs['L'+comp] = np.sum(L,axis=0)
   
   fluxes = np.sum(fluxes,axis=0)
   return fluxes, Labs


def _get_flux_from_table(fits_ext,photbands,index=None,include_Labs=True):
   """
   Retrieve flux and flux ratios from an integrated SED table.
   
   @param fits_ext: fits extension containing integrated flux
   @type fits_ext: FITS extension
   @param photbands: list of photometric passbands
   @type photbands: list of str
   @param index: slice or index of rows to retrieve
   @type index: slice or integer
   @return: fluxes or flux ratios
   #@rtype: list
   """
   if index is None:
      index = slice(None) #-- full range
   fluxes = []
   for photband in photbands:
      try:
         if not is_color(photband):
            fluxes.append(fits_ext.data.field(photband)[index])
         else:
            system,color = photband.split('.')
            if '-' in color:
               band0,band1 = color.split('-')
               fluxes.append(fits_ext.data.field('%s.%s'%(system,band0))[index]/fits_ext.data.field('%s.%s'%(system,band1))[index])
            elif color=='M1':
               fv = fits_ext.data.field('STROMGREN.V')[index]
               fy = fits_ext.data.field('STROMGREN.Y')[index]
               fb = fits_ext.data.field('STROMGREN.B')[index]
               fluxes.append(fv*fy/fb**2)
            elif color=='C1':
               fu = fits_ext.data.field('STROMGREN.U')[index]
               fv = fits_ext.data.field('STROMGREN.V')[index]
               fb = fits_ext.data.field('STROMGREN.B')[index]
               fluxes.append(fu*fb/fv**2)
      except KeyError:
         logger.warning('Passband %s missing from table'%(photband))
         fluxes.append(np.nan*np.ones(len(fits_ext.data)))
   #-- possibly include absolute luminosity
   if include_Labs:
      fluxes.append(fits_ext.data.field("Labs")[index])
   fluxes = np.array(fluxes).T
   if index is not None:
      fluxes = fluxes
   return fluxes

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