Quickstart
==========

Searching for photometry
------------------------

Speedyfit can automatically download photometry from several of the large surveys, and many smaller studies available
on Vizier. By default only photometry from the following large and well calibrated surveys is obtained:

- Gaia DR2
- Sky Mapper DR1
- APASS DR9
- SDSS
- 2MASS
- WISE

You can download photometry for a given object as follows:

.. code-block:: bash

   speedyfit <object_name> -empty single --phot

This will create a setup file called object_name_single.yaml where you can setup the fitting parameters and a photometry
file called object_name.phot with the obtained photometry. The object name should be resolvable by simbad, or if that
is not possible, a J-type coordinate. For example, to obtain photometry of the star 'EO Ceti' the following object names
will work: 'EO Ceti', 'J012343.24-05 05 45.83' and 'J020.93019703-05.096064175' and any other alias that is recoginzed
by simbad.

Fitting the SED
---------------

The the -empty command, a default setup file is created where you can define the parameters of your fit. A default yaml
file for a single star fit looks as follows.

.. code-block:: yaml

   # photometry file with index to the columns containing the photbands, observations and errors
   objectname: <objectname>
   photometryfile: <photfilename>
   photband_index: band
   obs_index: flux
   err_index: eflux
   # parameters to fit and the limits on them in same order as parameters
   pnames: [teff, logg, rad, ebv]
   limits:
   - [20000, 60000]
   - [5.80, 5.80]
   - [0.01, 0.5]
   - [0, 0.10]
   # constraints on distance or other parameters if known
   constraints:
       parallax: [<plx>, <e_plx>]
   # constraints on derived properties as mass, luminosity, luminosity ratio  if known
   derived_limits: {}
   # path to the model grids with integrated photometry
   grids:
   - tmap
   # setup for the MCMC algorithm
   nwalkers: 100    # total number of walkers
   nsteps: 1000     # steps taken by each walker (not including burn-in)
   nrelax: 250      # burn-in steps taken by each walker
   a: 10            # relative size of the steps taken
   # set the percentiles for the error determination
   percentiles: [16, 50, 84] # 16 - 84 corresponds to 1 sigma
   # output options
   resultfile: <objectname>_results_single.csv   # filepath to write results
   plot1:
      type: sed_fit
      result: pc
      path: <objectname>_sed_single.png
   plot2:
      type: distribution
      show_best: true
      path: <objectname>_distribution_single.png
      parameters: ['teff', 'rad', 'L', 'ebv', 'd', 'mass']


In that file you can setup the photometry, the parameters and ranges to include in the fit, any constraints that are
known, the setup of the MCMC algorithm and what output you like to have.

to run a fit:

.. code-block::

   speedyfit <object_name_setupfile.yaml>