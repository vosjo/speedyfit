
import pylab as pl
import numpy as np
from astropy.io import fits
from speedyfit.model import get_grid_file, get_table_single

# grid = 'blackbody'
#
# hdu = fits.open(get_grid_file(grid=grid, integrated=True))
#
# teff = hdu[1].data['teff']
# logg = hdu[1].data['logg']
#
# pl.figure()
# pl.plot(teff, logg, 'oC0')
# pl.xlabel('Effective temperature [K]')
# pl.ylabel('log(surface gravity) [cgs]')
# pl.title('grid = ' + grid)
# pl.tight_layout()
#
# pl.savefig('figures/' + grid + '.png')
#
# pl.show()

grids = ['kurucz', 'munari', 'tmap', 'blackbody']

pl.figure(figsize=(9,4))

for grid in grids:

    wave, flux = get_table_single(teff=20000, logg=5.0, ebv=0.0, grid=grid)
    s = np.where((wave > 3000) & (wave<10000))
    pl.plot(wave[s], flux[s], label=grid)

pl.xlabel('Wavelength (A)')
pl.ylabel('Flux (erg/s/cm2)')
pl.legend(loc='best')
pl.tight_layout()
pl.savefig('figures/model_comparison.png')
pl.show()
