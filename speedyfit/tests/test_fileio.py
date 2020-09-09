import pytest

import os

import numpy as np

from speedyfit import fileio

default = """
# photometry file and column indices
photometryfile: path/to/file.dat
photband_index: 0
obs_index: 1
err_index: 2
# parameters to fit and the limits on them in same order as parameters
pnames: [teff, logg, rad, teff2, logg2, rad2, ebv]
limits:
- [3500, 6000]
- [3.5, 5.0]
- [0.7, 1.5]
- [20000, 40000]
- [4.5, 6.5]
- [0.05, 0.3]
- [0, 0.02]
# constraints on distance and mass ratio if known
constraints:
  q: [3.03, 0.2]
  distance: [600, 50] # in parsec
# constraints on derived properties
derived_limits:
  mass: [0.5, 1.0]
  mass2: [0.1, 1.0]
# path to the model grids with integrated photometry
grids:
- path/to/grid/1.fits
- path/to/grid/2.fits
# setup for the MCMC algorithm
nwalkers: 100    # total number of walkers
nsteps: 2000     # steps taken by each walker (not including burn-in)
nrelax: 500      # burn-in steps taken by each walker
a: 10            # relative size of the steps taken
# output options
datafile: none   # filepath to write results of all walkers
"""


class TestFileIO:

    def test_write2fits(self):

        data = np.array([(1, 2), (2, 3)], dtype=[('a', 'f8'), ('b', 'f8')])

        try:
            fileio.write2fits(data, 'testfile.fits', setup=default)

            samples, setup = fileio.read_fits('testfile.fits')

            assert setup == default

            for name in data.dtype.names:
                for v1, v2 in zip(data[name], samples[name]):
                    assert v1 == v2

            os.remove('testfile.fits')

        finally:
            if os.path.isfile('testfile.fits'):
                os.remove('testfile.fits')
