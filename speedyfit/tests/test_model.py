import pytest

import numpy as np

from speedyfit import model

import os


class TestModel:

    def test_itable_single(self):

        print('OS: ' , os.environ.get('SPEEDYFIT_MODELS', None))

        photbands = ['STROMGREN.U', 'STROMGREN.B']
        gridfilename = 'kurucz2'

        axis_values, grid_pars, pixelgrid, grid_names = model.prepare_grid(photbands, gridfilename,
                                                                           teffrange=(4000, 5000), loggrange=(3.5, 4.0),
                                                                           ebvrange=(0.0, 0.02),
                                                                           variables=['teff', 'logg', 'ebv'])

        grid = [axis_values, pixelgrid]

        # -- check if interpolated model is correct
        flux1, Labs1 = model.get_itable_single(teff=4750, logg=3.75, ebv=0.01, grid=grid)

        assert Labs1 == pytest.approx(0.45605564, 0.001), \
            'Basic model: Luminosity is not correct ({} != 0.456)'.format(Labs1)

        # -- check if radius is treated correctly
        flux1, Labs1 = model.get_itable_single(teff=4750, logg=3.75, ebv=0.01, rad=1, grid=grid)
        assert Labs1 == pytest.approx(0.45605564, 0.001), \
            'Solar Radius: Luminosity is not correct ({} != 0.456)'.format(Labs1)

        flux1, Labs1 = model.get_itable_single(teff=4750, logg=3.75, ebv=0.01, rad=2, grid=grid)
        assert Labs1 == pytest.approx(2 ** 2 * 0.45605564, 0.001), \
            'Double Radius: Luminosity is not correct ({} != 1.824)'.format(Labs1)

        # -- compare with get_itable
        flux1, Labs1 = model.get_itable_single(teff=4750, logg=3.75, ebv=0.01, grid=grid)
        flux2, Labs2 = model.get_itable(teff=4750, logg=3.75, ebv=0.01, grid=grid)

        for f1, f2 in zip(flux1, flux2):
            assert f1 == f2, \
                "flux output of get_itable should be same as get_itable_single, ({} != {})".format(flux1, flux2)

        assert Labs1[0] == Labs2['L'], \
            "Labs output of get_itable should be same as get_itable_single, ({} != {})".format(Labs1[0], Labs2['L'])

    def test_itable(self):

        photbands = ['STROMGREN.U', 'STROMGREN.B']
        gridfilename = 'kurucz2'

        axis_values, grid_pars, pixelgrid, grid_names = model.prepare_grid(photbands, gridfilename,
                                                                           teffrange=(4000, 5000), loggrange=(3.5, 4.0),
                                                                           ebvrange=(0.0, 0.02),
                                                                           variables=['teff', 'logg', 'ebv'])

        grid1 = [axis_values, pixelgrid]

        gridfilename = 'tmap'

        axis_values, grid_pars, pixelgrid, grid_names = model.prepare_grid(photbands, gridfilename,
                                                                           teffrange=(20000, 40000),
                                                                           loggrange=(5.0, 6.5),
                                                                           ebvrange=(0.0, 0.02),
                                                                           variables=['teff', 'logg', 'ebv'])

        grid2 = [axis_values, pixelgrid]

        flux1, Labs1 = model.get_itable(teff=4750, logg=3.75, ebv=0.01, rad=1, grid=grid1)
        flux2, Labs2 = model.get_itable(teff=4250, logg=3.50, ebv=0.01, rad=1, grid=grid1)

        fluxt, Labst = model.get_itable(teff=4750, logg=3.75, ebv=0.01, rad=1,
                                        teff2=4250, logg2=3.50, ebv2=0.01, rad2=1, grid=[grid1, grid1])

        for f1, f2, ft in zip(flux1, flux2, fluxt):
            assert f1 + f2 == ft, \
                "binary SED should have sum of flux of single SED, ({} + {} != {})".format(flux1, flux2, fluxt)

        flux1, Labs1 = model.get_itable(teff=4800, logg=3.8, ebv=0.01, rad=2, grid=grid1)
        flux2, Labs2 = model.get_itable(teff=25000, logg=5.8, ebv=0.01, rad=0.14, grid=grid2)

        fluxt, Labst = model.get_itable(teff=4800, logg=3.8, ebv=0.01, rad=2,
                                        teff2=25000, logg2=5.8, ebv2=0.01, rad2=0.14,
                                        grid=[grid1, grid2])

        for f1, f2, ft in zip(flux1, flux2, fluxt):
            assert f1 + f2 == ft, \
                "binary SED should have sum of flux of single SED, ({} + {} != {})".format(flux1, flux2, fluxt)

    def test_itable_vega(self):
        """
        Test the returned flux for Vega
        Teff = 9602 +- 180 K
        logg = 4.10 +- 0.1
        Rad = 2.362 Rsol
        d = 7.68 +- 0.02 pc
        """

        photbands = ['STROMGREN.U', 'STROMGREN.B', 'STROMGREN.V', 'STROMGREN.Y', '2MASS.KS']
        gridfilename = 'kurucz2'

        axis_values, grid_pars, pixelgrid, grid_names = model.prepare_grid(photbands, gridfilename,
                                                                           teffrange=(8000, 10000),
                                                                           loggrange=(4.0, 5.0),
                                                                           ebvrange=(0.0, 0.02),
                                                                           variables=['teff', 'logg', 'ebv'])

        grid = [axis_values, pixelgrid]

        syn, Labs = model.get_itable_single(teff=9602, logg=4.10, ebv=0.0, rad=2.362, grid=grid)

        syn = syn / 340745610.888946 ** 2  # correct for distance
        Labs = Labs[0]

        obs = np.array([8.528e-10, 5.664e-09, 6.190e-09, 3.547e-09, 3.803e-11])
        obs_e = np.array([3.458e-11, 1.484e-10, 1.788e-10, 9.148e-11, 6.515e-12])

        ratio = (obs / syn)
        weights = (obs / obs_e)
        scale = np.average(ratio, weights=weights)
        e_scale = np.sqrt(np.dot(weights, (ratio - scale) ** 2) / weights.sum())

        assert abs(scale - 1.0) < e_scale, \
            "Cannot reproduce VEGA at correct distance! scale = {:0.1f}+-{:0.1f} instead of 1".format(scale, e_scale)

        assert abs(Labs - 40) < 2.5, "Vega should have luminosity of 40.12 Lsol, model gives {:0.2f} ".format(Labs)

    def test_itable_EG50(self):
        """
        Test the returned flux for a known white dwarf star EG 50
        Teff = 22290 +- 350 K
        logg = 8.10 +- 0.05
        Rad = 0.012 +- 0.001 Rsol
        d = 19 +- 3 pc
        """

        photbands = ['STROMGREN.U', 'STROMGREN.V', 'STROMGREN.B', 'STROMGREN.Y', '2MASS.H']
        gridfilename = 'tmap'

        axis_values, grid_pars, pixelgrid, grid_names = model.prepare_grid(photbands, gridfilename,
                                                                           teffrange=(20000, 40000),
                                                                           loggrange=(5.0, 6.5),
                                                                           ebvrange=(0.0, 0.02),
                                                                           variables=['teff', 'logg', 'ebv'])

        grid = [axis_values, pixelgrid]

        syn, Labs = model.get_itable_single(teff=22290, logg=6.5, ebv=0.0, rad=0.012, grid=grid)

        syn = syn / 842950390.9165244 ** 2  # correct for distance

        obs = np.array([4.89e-14, 1.02e-13, 9.31e-14, 5.34e-14, 9.76e-16])
        obs_e = np.array([2.53e-15, 3.99e-15, 3.08e-15, 1.72e-15, 5.48e-17])

        ratio = (obs / syn)
        weights = (obs / obs_e)
        scale = np.average(ratio, weights=weights)
        e_scale = np.sqrt(np.dot(weights, (ratio - scale) ** 2) / weights.sum())

        assert abs(scale - 1.0) < e_scale, \
            "Cannot reproduce EG 50 at correct distance! scale = {:0.1f}+-{:0.1f} instead of 1".format(scale, e_scale)
