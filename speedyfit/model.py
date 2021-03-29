import re
import os
import yaml
import numpy as np

from astropy.io import fits

from speedyfit import interpol, filters

__defaults__ = dict(grid='kurucz',
                    directory=os.environ.get('SPEEDYFIT_MODELS', None), )
defaults = __defaults__.copy()

# load a list of all available integrated grids
try:
    ifile = open("{}/grid_description.yaml".format(defaults['directory']))
    grid_description = yaml.safe_load(ifile)
    ifile.close()
except:
    grid_description = {}


def check_grids(print_bands=False):
    """
    Check which atmospheric model grids are installed and can be used by speedyfit.

    :param print_bands: If True, also print the photometric pass bands that are available in the integrated grid.
    :type print_bands: bool
    """
    print("Checking which atmosphere models are available...")

    if defaults['directory'] is not None:
        print("Checking for models in {}".format(defaults['directory']))
    else:
        print("SPEEDYFIT_MODELS environmental variable not set. CAN NOT find models!")
        print("Please point the SPEEDYFIT_MODELS variable to the directory where you stored models. On bash use:")
        print("\texport SPEEDYFIT_MODELS='<path to extracted atmosphere models>'")
        return

    if len(grid_description.keys()) == 0:
        print("grid_description.yaml file not found or no models included in the description file.")
        print("Please add a grid_description.yaml file in the SPEEDYFIT_MODELS directory.")
        return

    for grid in grid_description.keys():
        print(grid)
        gridpath = get_grid_file(integrated=False, grid=grid)
        if os.path.isfile(gridpath):
            print("\t raw: available")
        else:
            print("\t raw: NOT FOUND")

        gridpath = get_grid_file(integrated=True, grid=grid)
        if os.path.isfile(gridpath):
            print("\t integrated: available")
        else:
            print("\t integrated: NOT FOUND")

        if print_bands:
            hdu = fits.open(gridpath)
            bands = set([b.split('.')[0] for b in hdu[1].data.dtype.names if '.' in b])
            for b in bands:
                print("\t - " + b)


        if 'info' in grid_description[grid]:
            print('\t info: ' + grid_description[grid]['info'])


def get_grid_file(integrated=False, **kwargs):
    grid = kwargs.get('grid', defaults['grid'])

    if os.path.isfile(grid):
        return grid

    if grid in grid_description:
        filename = grid_description[grid]['filename']
    else:
        raise ValueError('Grid name ({}) not recognized!'.format(grid))

    if integrated:
        filename = 'i' + filename + '_lawfitzpatrick2004_Rv3.10'

    directory = kwargs.get('directory', defaults['directory'])

    return directory + filename + '.fits'


def get_grid_ranges(**kwargs):
    grid = kwargs.get('grid', defaults['grid'])

    if grid == 'kurucz' or grid == 'kurucz2': # kurucz2 is necessary for backwards compatibility
        teff = (3500, 10000)
        logg = (4.32, 4.32)
        rad = (0.05, 2.5)

    elif grid == 'munari':
        teff = (3500, 20000)
        logg = (4.32, 4.32)
        rad = (0.05, 2.5)

    elif grid == 'tmap':
        teff = (20000, 80000)
        logg = (5.88, 5.88)
        rad = (0.01, 0.5)

    elif grid == 'koester':
        teff = (20000, 80000)
        logg = (5.88, 5.88)
        rad = (0.01, 0.5)

    elif grid == 'blackbody':
        teff = (20000, 80000)
        logg = (5.88, 5.88)
        rad = (0.01, 2.5)

    else:
        raise ValueError('Grid name ({}) not recognized!'.format(grid))

    return {'teff': teff, 'logg': logg, 'rad': rad}


def get_grid_dimensions(**kwargs):
    """
    Retrieve possible effective temperatures and gravities from a grid.

    E.g. kurucz, sdB, fastwind...

    :rtype: (ndarray,ndarray)
    :return: effective temperatures, gravities
    """
    gridfile = get_grid_file(**kwargs)
    ff = fits.open(gridfile)
    teffs = []
    loggs = []
    for hdu in ff[1:]:
        teffs.append(float(hdu.header['TEFF']))
        loggs.append(float(hdu.header['LOGG']))
    ff.close()

    # # -- maybe the fits extensions are not in right order...
    # matrix = np.vstack([np.array(teffs), np.array(loggs)]).T
    # matrix = numpy_ext.sort_order(matrix, order=[0, 1])
    # teffs, loggs = matrix.T

    return teffs, loggs


def load_grids(gridnames, pnames, limits, photbands):
    """
    prepares the integrated photometry grid by loading the grid and cutting it to the size
    given in limits.
    """
    grids = []
    for i, name in enumerate(gridnames):
        ind = '' if i == 0 else str(i + 1)

        axis_values, grid_pars, pixelgrid, grid_names = prepare_grid(photbands, name,
                                                                     teffrange=limits[pnames.index('teff' + ind)],
                                                                     loggrange=limits[pnames.index('logg' + ind)],
                                                                     ebvrange=limits[pnames.index('ebv')],
                                                                     variables=['teff', 'logg', 'ebv'])

        grids.append([axis_values, pixelgrid])

    return grids


def prepare_grid(photbands, gridname,
                 teffrange=(-np.inf, np.inf), loggrange=(-np.inf, np.inf),
                 ebvrange=(-np.inf, np.inf),
                 variables=['teff', 'logg', 'ebv'], **kwargs):
    flux = []
    grid_pars = []
    grid_names = np.array(variables)

    gridfilename = get_grid_file(integrated=True, grid=gridname)

    with fits.open(gridfilename) as ff:
        # -- make an alias for further reference
        ext = ff[1]

        # -- the grid is already cut here to limit memory usage
        keep = np.ones(len(ext.data), bool)
        for name in variables:
            # -- first find the closest actuall grid points
            low, high = locals()[name + 'range']
            lidx = np.abs(ext.data.field(name)[ext.data.field(name) <= low] - low).argmin()
            hidx = np.abs(ext.data.field(name)[ext.data.field(name) >= high] - high).argmin()
            low = ext.data.field(name)[ext.data.field(name) <= low][lidx]
            high = ext.data.field(name)[ext.data.field(name) >= high][hidx]

            # -- we need to be carefull for rounding errors
            in_range = (low <= ext.data.field(name)) & (ext.data.field(name) <= high)
            on_edge = np.allclose(ext.data.field(name), low) | np.allclose(ext.data.field(name), high)

            keep = keep & (in_range | on_edge)

        grid_pars = np.vstack([ext.data.field(name)[keep] for name in variables])
        flux = _get_flux_from_table(ext, photbands, include_Labs=True)[keep]

    flux = np.log10(flux)

    # -- create the pixeltype grid
    axis_values, pixelgrid = interpol.create_pixeltypegrid(grid_pars, flux.T)
    return axis_values, grid_pars.T, pixelgrid, grid_names


def get_itable_single(teff=None, logg=None, g=None, ebv=0.0, **kwargs):
    # -- check if logg or g is given
    if logg is None and not g is None:
        logg = np.log10(g)

    multiple = False
    if not hasattr(teff, '__iter__') or not hasattr(logg, '__iter__') or \
            not hasattr(ebv, '__iter__'):
        teff = np.array(teff)
        logg = np.array(logg)
        ebv = np.array(ebv)
        multiple = True

    # -- get the grid from the keyword, or prepare it if a grid name is given
    if isinstance(kwargs['grid'], str):
        axis_values, grid_pars, pixelgrid, grid_names = prepare_grid(kwargs['photbands'], kwargs['grid'],
                                                                     teffrange=(np.min(teff), np.max(teff)),
                                                                     loggrange=(np.min(logg), np.max(logg)),
                                                                     ebvrange=(np.min(ebv), np.max(ebv)),
                                                                     variables=['teff', 'logg', 'ebv'])
    else:
        axis_values, pixelgrid = kwargs['grid']

    p = np.vstack([teff, logg, ebv])

    values = interpol.interpolate(p, axis_values, pixelgrid)

    # -- switch logarithm to normal
    values = 10 ** values
    flux, Labs = values[:-1], values[-1]

    # -- Take radius into account when provided
    if 'rad' in kwargs:
        rad = np.array(kwargs['rad'])
        flux, Labs = flux * rad ** 2, Labs * rad ** 2

    if 'd' in kwargs:
        d = np.array(kwargs['d'])
        flux, Labs = flux / d ** 2, Labs * d ** 2

    if multiple:
        flux = flux.flatten()
        Labs = Labs.flatten()

    return flux, Labs


def get_itable(grid=[], **kwargs):
    values, parameters, components = {}, set(), set()
    for key in list(kwargs.keys()):
        if re.search("^(teff|logg|g|ebv|rad)\d?$", key):
            par, comp = re.findall("^(teff|logg|g|ebv|rad)(\d?)$", key)[0]
            values[key] = kwargs.pop(key)
            parameters.add(par)
            components.add(comp)

    # need to check here that grid is same length as components
    if isinstance(grid, str):
        grids = [grid]
    elif hasattr(grid, '__iter__') and isinstance(grid[0], str):
        grids = grid
    elif hasattr(grid, '__iter__') and hasattr(grid[0], '__iter__') and len(grid[0]) == 2:
        grids = grid
    else:
        grids = [grid]

    # -- If there is only one component, we can directly return the result
    if len(components) == 1:
        kwargs.update(values)
        fluxes, Labs = get_itable_single(grid=grids[0], **kwargs)
        return fluxes, {'L': Labs}

    fluxes, Labs = [], {}
    for i, (comp, grid) in enumerate(zip(components, grids)):
        kwargs_ = kwargs.copy()
        for par in parameters:
            kwargs_[par] = values[par + comp] if par + comp in values else values[par]

        f, L = get_itable_single(grid=grid, **kwargs_)

        fluxes.append(f)
        Labs['L' + comp] = np.sum(L, axis=0)

    fluxes = np.sum(fluxes, axis=0)
    return fluxes, Labs


def get_table_single(teff=None, logg=None, ebv=0.0, **kwargs):
    """
   No interpolating, just returns the closest gridpoint
   """

    # -- get the grid
    gridname = kwargs['grid']
    gridfilename = get_grid_file(integrated=False, grid=gridname)

    hdu = fits.open(gridfilename)

    teffs, loggs = np.zeros(len(hdu) - 1), np.zeros(len(hdu) - 1)
    for i in range(1, len(hdu)):
        teffs[i - 1] = hdu[i].header['TEFF']
        loggs[i - 1] = hdu[i].header['LOGG']

    dteff = abs(teffs - teff) / teff
    dlogg = abs(loggs - logg) / logg

    s = np.where(np.sqrt(dteff ** 2 + dlogg ** 2) == np.min(np.sqrt(dteff ** 2 + dlogg ** 2)))

    model = hdu[s[0][0] + 1].data
    wave, flux = model['wavelength'], model['flux']

    # -- Take radius into account when provided
    if 'rad' in kwargs:
        rad = np.array(kwargs['rad'])
        flux = flux * rad ** 2

    if 'd' in kwargs:
        d = np.array(kwargs['d'])
        flux = flux / d ** 2

    return wave, flux


def get_table(grid=[], **kwargs):
    """
   Returns the closest model atmosphere available in the grid. No interpolation is done!
   """
    values, parameters, components = {}, set(), set()
    for key in list(kwargs.keys()):
        if re.search("^(teff|logg|g|ebv|rad)\d?$", key):
            par, comp = re.findall("^(teff|logg|g|ebv|rad)(\d?)$", key)[0]
            values[key] = kwargs.pop(key)
            parameters.add(par)
            components.add(comp)

    # need to check here that grid is same length as components
    if hasattr(grid, '__iter__') and not isinstance(grid, str):
        grids = grid
    else:
        grids = [grid]

    # -- If there is only one component, we can directly return the result
    if len(components) == 1:
        kwargs.update(values)
        wave, flux = get_table_single(grid=grids[0], **kwargs)
        return wave, flux

    waves, fluxes = [], []
    for i, (comp, grid) in enumerate(zip(components, grids)):
        kwargs_ = kwargs.copy()
        for par in parameters:
            kwargs_[par] = values[par + comp] if par + comp in values else values[par]

        w, f = get_table_single(grid=grid, **kwargs_)

        waves.append(w)
        fluxes.append(f)

    # -- interpolate and combine the models
    wave = waves[0]
    flux = np.zeros_like(wave)
    for w, f in zip(waves, fluxes):
        flux += np.interp(wave, w, f)

    return wave, flux


def _get_flux_from_table(fits_ext, photbands, index=None, include_Labs=True):
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
        index = slice(None)  # -- full range
    fluxes = []
    for photband in photbands:
        try:
            if not filters.is_color(photband):
                fluxes.append(fits_ext.data.field(photband)[index])
            else:
                system, color = photband.split('.')
                if '-' in color:
                    band0, band1 = color.split('-')
                    fluxes.append(fits_ext.data.field('%s.%s' % (system, band0))[index] /
                                  fits_ext.data.field('%s.%s' % (system, band1))[index])
                elif color == 'M1':
                    fv = fits_ext.data.field('STROMGREN.V')[index]
                    fy = fits_ext.data.field('STROMGREN.Y')[index]
                    fb = fits_ext.data.field('STROMGREN.B')[index]
                    fluxes.append(fv * fy / fb ** 2)
                elif color == 'C1':
                    fu = fits_ext.data.field('STROMGREN.U')[index]
                    fv = fits_ext.data.field('STROMGREN.V')[index]
                    fb = fits_ext.data.field('STROMGREN.B')[index]
                    fluxes.append(fu * fb / fv ** 2)
        except KeyError:
            print('Passband %s missing from table' % (photband))
            fluxes.append(np.nan * np.ones(len(fits_ext.data)))
    # -- possibly include absolute luminosity
    if include_Labs:
        fluxes.append(fits_ext.data.field("Labs")[index])
    fluxes = np.array(fluxes).T
    if index is not None:
        fluxes = fluxes
    return fluxes


def luminosity(wave, flux, radius=1.):
    """
    Calculate the bolometric luminosity of a model SED.

    Flux should be in cgs per unit wavelength (same unit as wave).
    The latter is integrated out, so it is of no importance. After integration,
    flux, should have units erg/s/cm2.

    Returned luminosity is in solar units.

    If you give radius=1 and want to correct afterwards, multiply the obtained
    Labs with radius**2.

    :param wave: model wavelengths
    :type wave: ndarray
    :param flux: model fluxes (Flam)
    :type flux: ndarray
    :param radius: stellar radius in solar units
    :type radius: float
    :return: total bolometric luminosity
    :rtype: float
    """
    Lsol_cgs = 3.846e33
    Rsol_cgs = 6.95508e10
    Lint = np.trapz(flux, x=wave)
    Labs = Lint * 4 * np.pi / Lsol_cgs * (radius * Rsol_cgs) ** 2
    return Labs