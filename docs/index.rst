.. speedyfit documentation master file, created by
   sphinx-quickstart on Mon Sep 14 09:27:51 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Speedyfit
=========

Speedyfit is a python package to fit the photometric spectral energy distribution of stars. It uses Bayesian inference
to correctly take known parameters with errorbars into account as priors, and uses a a Markov chain Monte Carlo approach
to fit the requested parameters and determine their final errors.

It is a command line tool writen in Python 3 that allows you to search the most common online databases for
photometric observations of your target, and fit theoretical atmosphere models to the obtained photometry. Speedyfit can
deal with both single and binary stars, and allows for the inclusion of constraints from other sources, as for example
the distance or reddening but equally so constraints from spectroscopic observations.

.. toctree::
   :maxdepth: 2
   :caption: User Guide:
   
   userguide/install
   userguide/quickstart
   userguide/command_line
   userguide/setup_file
   userguide/photometry_catalogs
   userguide/photometry_file
   userguide/model_grids
   userguide/making_figures


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
