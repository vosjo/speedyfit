The setup file
==============

All details regarding the setup for the SED fit are provided to Speedyfit using a setup file in yaml format. Yaml was
chosen because it is easy to read by a human and by a computer. The default setup file with all its options are given
below:

.. code-block:: yaml

    # photometry file with index to the columns containing the photbands, observations and errors
    objectname: Vega
    photometryfile: Vega.phot
    photband_index: band
    obs_index: flux
    err_index: eflux
    photband_include: ['GAIA', 'APASS', '2MASS']
    photband_exclude: ['GALEX']
    # parameters to fit and the limits on them in same order as parameters
    pnames: [teff, logg, rad, ebv]
    limits:
    - [3500, 10000]
    - [4.31, 4.31]
    - [0.01, 2.5]
    - [0, 0.10]
    # constraints on distance, mass ratio or any of the fitted parameters when known
    constraints:
        parallax: [130.23, 0.36]
    # added constraints on derived properties as mass, luminosity, luminosity ratio
    derived_limits:
        mass: [1.5, 2.0]
    # path to the model grids with integrated photometry
    grids:
    - kurucz
    # setup for the MCMC algorithm
    nwalkers: 100    # total number of walkers
    nsteps: 1000     # steps taken by each walker (not including burn-in)
    nrelax: 250      # burn-in steps taken by each walker
    a: 10            # relative size of the steps taken
    # set the percentiles for the error determination
    percentiles: [16, 50, 84] # 16 - 84 corresponds to 1 sigma
    # output options
    resultfile: Vega_results.csv   # filepath to write results
    plot1:
     type: sed_fit
     result: pc
     path: Vega_sed.png
    plot2:
     type: distribution
     show_best: true
     path: Vega_distribution.png
     parameters: ['teff', 'rad', 'L', 'ebv', 'd', 'mass']

Photometry setup
----------------

.. option:: objectname (str)

    Name of the system being fitted. Can be whatever you want, and doesn't need to be resolvable by Simbad. It is used
    for the figures and can help in recognizing what model this file belongs too, but doesn't serve any other purpose.

.. option:: photometryfile (str)

    The path to the file where the photometry of your object is stored.

.. option:: photband_index (str or int)

    The name or index of the column in the photometry file that contains the names of the photometric bands. See
    :doc:`photometry_file` for more details.

.. option:: obs_index (str or int)

    The name or index of the column in the photometry file that contains the flux. See :doc:`photometry_file` for more
    details.

.. option:: err_index (str or int)

    The name or index of the column in the photometry file that contains the error on the flux. See
    :doc:`photometry_file` for more details.

.. option:: photband_include (list)

    A list of all the photometric bands to include in the fit. Any band not mentioned in this list is excluded. Bands
    are named as: SYSTEM.BAND (e.g. GAIA2.G for the G band from Gaia DR2). When the system is used (e.g. GAIA2) all
    bands from that system are included. If a specific band (e.g. GAIA2.G) is used, only that band is included.
    When both :option:`photband_exclude` and :option:`photband_include` are included in the setup file, photband_include
    is processed first, after which bands to be excluded are removed.

.. option:: photband_exclude (list)

    A list of which photometric bands to exclude from the fit. All bands mentioned here will be excluded. Bands are
    named as: SYSTEM.BAND (e.g. GAIA2.G for the G band from Gaia DR2). When the system is used (e.g. GAIA2) all
    bands from that system are excluded. If a specific band (e.g. GAIA2.G) is mentioned, only that band is excluded.
    When both :option:`photband_exclude` and :option:`photband_include` are included in the setup file, photband_include
    is processed first, after which bands to be excluded are removed.

Fit and constraints setup
-------------------------

.. option:: pnames (list)

    A list of the parameters to include in the fit. If you are performing a single fit, this will be:

    [teff, logg, rad, ebv]

    If you are performing a binary fit this will be:

    [teff, logg, rad, teff2, logg2, rad2, ebv]

    Note that for a binary fit, the reddening parameter (ebv) will only be included once. Both components in the fit
    need to have the same reddening.

.. option:: limits (list of lists)

    The limits on the parameters included in the fit. Limits have to be given in the same order as the parameters are
    given in :option:`pnames`.

    For example:

    | pnames = [teff, logg]
    | limits = [(3500, 10000), (3.5, 4.5)]

    means that effective temperature is limited between 3500 and 10000 K, and logg is limited between 3.5 and 4.5 dex.

.. option:: constraints (dict)

    A dictionary containing all the constraints that you want to apply in the fit. These are treated as priors in the
    Bayesian fit. You can apply a prior on the distance/parallax, any of the parameters you are fitting (teff, logg,
    rad, ebv), any of the derived parameters (mass and luminosity) and in case of a binary fit also the mass ratio.

    A constraint is not the same as a limit. A constraint consists of a value and an error. You can provide a double
    sided error as well. e.g.

    | parallax: [130.23, 0.36] -> parallax = 130.23 +- 0.36
    | teff: [5630, 150, 270] -> teff = 5630 -150 +270

    Note: there constraints can be valuable if for example the effective temperature and surface gravity are
    known from a spectroscopic fit, and you want to derive the radius from the SED and the parallax while correctly
    propagating the errors.

.. option:: derived_limits (dict)

    A dictionary containing all the limits that you want to apply to derived parameters. These are also treated as
    priors, but as limiting priors were the prior is uniform between the lower and upper limit. You can apply limits
    on all derived parameters (mass and luminosity), on the distance and in case of a binary fit also on the mass
    ratio (q). e.g:

    q: [0.5, 0.9]

    means that the mass ratio is limited between 0.5 and 0.9

.. option:: grids (list)

    List of the model atmosphere grids used in the fit. For a single star fit only one grid name should be provided,
    for a binary fit, two grid names should be provided. You can provide the name of a grid as recognized by speedyfit,
    or the direct path to the grid. Recognized names are:

    - kurucz
    - munari
    - tmap
    - koester
    - blackbody

    For more info on the different available grids see :doc:`model_grids`.

.. option:: percentiles (list)

    The percentiles to use when calculating the error on the fitted parameters. By default the percentiles are set to:

    [16, 50, 85]

    The best fitting value is the median of all values, and the lower and upper error given by speedyfit will then
    correspond to 1 sigma. This means that 68% of all steps taken by the walkers are within the error given by
    Speedyfit. See `the 68–95–99.7 rule <https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule>`_ for more
    information on this. Quick summary:

    | 1 sigma: [16, 50, 85]
    | 2 sigma: [2.5, 50, 97.5]
    | 3 sigma: [0.15, 50, 99.85]

    If you don't do any burn-in (via option :option:`nrelax`), setting the percentiles to correspond with for example
    3 sigma will result in very large unrealistic errors.

MCMC setup
----------

The MCMC chain is implemented using the EMCEE package. See `emcee <https://emcee.readthedocs.io/>`_ for more
information on the parameters governing the MCMC fitting process.

.. option:: nwalkers (int)

    Number of walkers used in the MCMC chain.

.. option:: nsteps (int)

    Number of steps that each walker takes in the MCMC chain.

.. option:: nrelax (int)

    Number of burn in steps that each walker takes. Walkers will start of at random locations in the parameter space
    and will slowly converge to an optimal solution. The nrelax steps are the first number of steps that will not be
    taken into account when determining the errors on the parameters from the walker chains.

.. option:: a (int)

    Relative step size taken by each walker.


Output and figures
------------------

.. option:: resultfile (str)

    The name of the file to store the results of the fit. This file will contain the best fitting value and error of
    each fitted and derived parameters in csv format.

.. option:: plot<n> (dict)

    In the result file you can also define what figures you want speedyfit to make. There are 3 types of figures that
    speedyfit can make:

    - sed_fit
    - constraint
    - distribution

    In the setup file you can define up to 10 different configurations of those 3 figures. All figures will follow the same
    structure:

    .. code-block:: yaml

        plot1: # number of the plot, up to 10
            type: <plot type> # which type of plot you want
            path: <figurename>.png # filename in which to store the figure
            parameter: value # any other parameters to configure this plot

    For more information on available plots and their configuration see :doc:`making_figures`
