import numpy as np
from scipy import ndimage

def create_pixeltypegrid(grid_pars,grid_data):
   """
   Creates pixelgrid and arrays of axis values.
   
   Starting from:
      * grid_pars: 2D numpy array, 1 column per parameter, unlimited number of cols
      * grid_data: 2D numpy array, 1 column per variable, data corresponding to the rows in grid_pars
   
   example: interpolation in a 3D grid containing stellar evolution models. Say we have as
   input parameters mass, age and metalicity, and want to obtain teff and logg as variables.
   
   grid_pars =
      +------+-----+------+
      | mass | age | Fe/H |
      +------+-----+------+
      | 1.0  | 1.0 | -0.5 |
      +------+-----+------+
      | 2.0  | 1.0 | -0.5 |
      +------+-----+------+
      | 1.0  | 2.0 | -0.5 |
      +------+-----+------+
      | 2.0  | 2.0 | -0.5 |
      +------+-----+------+
      | 1.0  | 1.0 |  0.0 |
      +------+-----+------+
      | 2.0  | 1.0 |  0.0 |
      +------+-----+------+
      |...   |...  |...   |
      +------+-----+------+
      
   grid_data = 
      +------+------+
      | teff | logg |
      +------+------+
      | 5000 | 4.45 |
      +------+------+
      | 6000 | 4.48 |
      +------+------+
      |...   |...   |
      +------+------+
   
   The resulting grid will be rectangular and complete. This means that every
   combination of unique values in grid_pars should exist. If this is not the
   case, a +inf value will be inserted in grid_data at all locations that are 
   missing!

   
   :param grid_pars: Npar x Ngrid array of parameters
   :type grid_pars: array
   :param grid_data: Ndata x Ngrid array of data
   :type grid_data: array
   
   :return: axis values and pixelgrid
   :rtype: array, array
   """

   uniques = [np.unique(column, return_inverse=True) for column in grid_pars]
   #[0] are the unique values, [1] the indices for these to recreate the original array

   axis_values = [uniques_[0] for uniques_ in uniques]
   unique_val_indices = [uniques_[1] for uniques_ in uniques]
   
   data_dim = np.shape(grid_data)[0]

   par_dims   = [len(uv[0]) for uv in uniques]

   par_dims.append(data_dim)
   pixelgrid = np.ones(par_dims)
   
   # We put np.inf as default value. If we get an inf, that means we tried to access
   # a region of the pixelgrid that is not populated by the data table
   pixelgrid[pixelgrid==1] = np.inf
   
   # now populate the multiDgrid
   indices = [uv[1] for uv in uniques]
   pixelgrid[indices] = grid_data.T
   
   return axis_values, pixelgrid

def interpolate(p, axis_values, pixelgrid):
   """
   Interpolates in a grid prepared by create_pixeltypegrid().
   
   p is an array of parameter arrays
   each collumn contains the value for the corresponding parameter in grid_pars
   each row contains a set of model parameters for wich the interpolated values
   in grid_data are requested.
   
   example: continue with stellar evolution models used in create_pixeltypegrid
   
   p = 
      +------+-----+-------+
      | mass | age | Fe/H  | 
      +------+-----+-------+
      | 1.21 | 1.3 | 0.24  |
      +------+-----+-------+
      | 1.57 | 2.4 | -0.15 |
      +------+-----+-------+
      |...   |...  |...    |
      +------+-----+-------+
      
   >>> p = np.array([[1.21, 1.3, 0.24], [1.57, 2.4, -0.15]])
   >>> interpolate(p, axis_values, pixelgrid)
   >>> some output
   
   :param p: Npar x Ninterpolate array containing the points which to
             interpolate in axis_values
   :type p: array
   :param axis_values: output from create_pixeltypegrid
   :type axis_values: array
   :param pixelgrid: output from create_pixeltypegrid
   :type pixelgrid: array
   
   :return: Ndata x Ninterpolate array containing the interpolated values
            in pixelgrid
   :rtype: array
   
   """
   # convert requested parameter combination into a coordinate
   #p_ = [np.searchsorted(av_,val) for av_, val in zip(axis_values,p)]
   # we force the values to be inside the grid, to avoid edge-effect rounding
   # (e.g. 3.099999 is edge, while actually it is 3.1). For values below the
   # lowest value, this is automatically done via searchsorted (it return 0)
   # for values higher up, we need to force it
   
   p_ = []
   for av_,val in zip(axis_values,p):
      indices = np.searchsorted(av_,val)
      indices[indices==len(av_)] = len(av_)-1
      p_.append(indices)
   

   #-- The type of p is changes to the same type as in axis_values to catch possible rounding errors
   #   when comparing float64 to float32.
   for i, ax in enumerate(axis_values):
      p[i] = np.array(p[i], dtype = ax.dtype)
   
   #-- Convert requested parameter combination into a coordinate
   p_ = np.array([np.searchsorted(av_,val) for av_, val in zip(axis_values,p)])
   lowervals_stepsize = np.array([[av_[p__-1], av_[p__]-av_[p__-1]] \
                           for av_, p__ in zip(axis_values,p_)])
   p_coord = (p-lowervals_stepsize[:,0])/lowervals_stepsize[:,1] + np.array(p_)-1

   # interpolate
   return np.array([ndimage.map_coordinates(pixelgrid[...,i],p_coord, order=1, prefilter=False) \
               for i in range(np.shape(pixelgrid)[-1])])

