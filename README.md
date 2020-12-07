![PyPI - Python Version](https://img.shields.io/pypi/pyversions/speedyfit)
[![PyPI](https://img.shields.io/pypi/v/speedyfit?color=blue)](https://pypi.org/project/speedyfit/)
[![Documentation Status](https://readthedocs.org/projects/speedyfit/badge/?version=latest)](https://speedyfit.readthedocs.io/en/latest/?badge=latest)
[![GitHub](https://img.shields.io/github/license/vosjo/speedyfit)](https://github.com/vosjo/speedyfit/blob/master/LICENSE)

# Speedyfit

A python package to fit the photometric spectral energy distribution of stars. Uses a Markov chain Monte Carlo approach 
to determine the errors on the derived parameters.

Speedyfit is a command line tool writen in Python 3 that allows you to search the most common online databases for 
photometric observations of your target, and fit theoretical atmosphere models to the obtained photometry. Speedyfit can
deal with both single and binary stars, and allows for the inclusion of constraints from other sources, as for example
the distance or reddening. 

## Documentation

The detailed documentation is a work in progress, but the command line documentation is complete. You can read the
 current version at:
[speedyfit.readthedocs.io](https://speedyfit.readthedocs.io/en/latest/)
