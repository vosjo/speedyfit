Command line options
====================

Speedyfit can run in three different modes:

- setup: create a setup file for an SED fit.
- phot: dowload photometry for an object.
- fit: fit a photometric SED based on a setup file.
- checkgrids: check which model atmosphere grids are installed.

Each mode has its own options, which are described below.

speedyfit setup
---------------

The 'setup' mode allows you to create a default setup file for the object that you want to fit. This setup file will
contain the detailed setup for the SED fit including the parameters to fit, their limits, where to find the photometry,
what model grids to use, the settings for the mcmc sampler and the output options. Speedyfit will create a default file
based on a few arguments that you can then further edit to your liking. More details on the different possible options
available in the setup file are given in :doc:`setup_file`.

**speedyfit** setup <*object_name*>  [*options*]

.. program:: speedyfit setup

.. option:: object_name

    The name of the object of which you want to fit the SED. This name needs to be resolvable by simbad, or it needs to
    be the coordinates in J<nnnnnn.n+-nnnnnn.n> format. As an example, consider the star Vega.

    The following names for Vega (ra = 18:36:56.33635, dec = +38:47:01.2802) will be recognized:

    - Vega
    - any other of the 62 identifiers recognized by simbad:
      `* alf Lyr: Identifiers <http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=vega&submit=SIMBAD+search#lab_ident>`_
    - J183656.3+384701.2

.. option:: -grid <model_grid>

    The name of the model grid that you want to use in the SED fit. This name is used to fill in the 'grid' parameter,
    and the 'limits' parameter. By default the limits are set roughly at the limits of the grid. As some grids are not
    rectangular, this is not always possible.

    If the option :option:`-grid` binary is used, a binary setup file will be created with as default a hot component
    fitted with the TMAP grid, and a cold component fitted with the Kurucz grid.

.. option:: --phot

    If provided, speedyfit will also download photometry from Vizier and Tap archives. If you don't specify this option,
    you can download the photometry in a separate step with 'speedyfit phot'.

.. option:: --nopx

    If provided, speedyfit setup will NOT include the Gaia DR 2 parallax as a constraint in the fitting process. By
    default the parallax is always included.

speedyfit phot
--------------

The 'photometry' or 'phot' mode allows you do download photometry from several vizier and tap catalogs for the objects
you are interested in. The default catalogs that Speedyfit will query can be found in :doc:`photometry_catalogs`. More
information on the formatting of the photometry file can be found in :doc:`photometry_file`. This command only needs an
object name to work.

**speedyfit** phot <*object_name*>  [*options*]

.. program:: speedyfit phot

.. option:: object_name

    The name of the object of which you want to download photometry. This name needs to be resolvable by simbad, or it
    needs to be the coordinates in J<nnnnnn.n+-nnnnnn.n> format. As an example, consider the star Vega.

    The following names for Vega (ra = 18:36:56.33635, dec = +38:47:01.2802) will be recognized:

    - Vega
    - any other of the 62 identifiers recognized by simbad:
      `* alf Lyr: Identifiers <http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=vega&submit=SIMBAD+search#lab_ident>`_
    - J183656.3+384701.2

.. option:: -o <output_file>, -output <output_file>

    You can specify the name of the file in which speedyfit will store the downloaded photometry. If no name is given,
    the default name is: <object_name>.phot

speedyfit fit
-------------

Finally, the 'fit' mode will use all information provided in the setup file together with the downloaded photometry to
fit the SED and store and display the results. This command requires a setup yaml file, and requires you to have
collected photometry of your target object and describe where it can be found in the setup file.

**speedyfit** fit <*setup_file*>  [*options*]

.. program:: speedyfit fit

.. option:: setup_file

    The name of the setup yaml file in which all parameters necessary for the SED fit process are specified. A default
    setup file can be created with the 'speedyfit setup' command.

.. option:: --noplot

    When provided, speedyfit fit will NOT display any of the plots, and will only save them to file. By default all
    plots requested in the setup file are shown when the fit is finished.

speedyfit checkgrids
--------------------

The 'checkgrids' mode will check and print which model atmosphere grids are installed and can be used by speedyfit.

**speedyfit** checkgrids [*options*]

.. program:: speedyfit checkgrids

.. option:: --bands

    When provided, speedyfit will also list the photometric bands included in the integrated grids.