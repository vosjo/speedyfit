import pytest

import numpy as np

from speedyfit import statfunc


class TestMCMC:

    def test_derived_properties(self):
        theta = (26000, 5.8, 0.1429, 5771, 4.438, 1, 0.02)
        pnames = ('teff', 'logg', 'rad', 'teff2', 'logg2', 'rad2', 'ebv')

        derived_properties = statfunc.get_derived_properties_binary(theta, pnames)

        assert derived_properties['mass'] == pytest.approx(0.47, 0.01), \
            "sdB mass wrongly calculated: {}".format(derived_properties['mass'])

        assert derived_properties['mass2'] == pytest.approx(1.00, 0.01), \
            "Solar mass wrongly calculated: {}".format(derived_properties['mass2'])

        assert derived_properties['q'] == pytest.approx(0.47, 0.001), \
            "q wrongly calculated: {}".format(derived_properties['q'])

    def test_stat_chi2(self):
        meas = np.array([1.0, 2.0, 3.0])
        e_meas = np.array([0.25, 0.30, 0.15])
        syn = np.array([1.5, 1.89, 3.6])

        # -- colors (none scaled)
        colors = np.array([True, True, True])
        chi2, scale, e_scale = statfunc.stat_chi2(meas, e_meas, colors, syn)

        assert chi2 == pytest.approx(20.1344, 0.0001), \
            "Chi2 for colours not correct ({} != {})".format(chi2, 20.1344)

        # -- absolute magnitudes (scaled)
        colors = np.array([False, False, False])
        chi2, scale, e_scale = statfunc.stat_chi2(meas, e_meas, colors, syn)

        assert chi2 == pytest.approx(3.32834, 0.0001), \
            "Chi2 for colours not correct ({} != {})".format(chi2, 3.3283)

        # -- added distance
        constraints = {'distance': (1.5, 0.15)}
        colors = np.array([False, False, False])
        chi2, scale, e_scale = statfunc.stat_chi2(meas, e_meas, colors, syn, **constraints)

        assert chi2 == pytest.approx(25.2343, 0.0001), "Chi2 for distance not correct ({} != {})".format(chi2, 25.2343)
