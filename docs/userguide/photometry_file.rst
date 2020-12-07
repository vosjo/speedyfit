The photometry file
===================

Speedyfit can download photometry from Vizier and use Tap queries to obtain photometry from Skymapper. It saves the
downloaded photometry to a text file using the
`astropy.io.ascii.write <https://docs.astropy.org/en/stable/api/astropy.io.ascii.write.html>`_ function using the
'fixed_width' format.

A photometry file will look like:

.. code-block:: bash

    |        band |  meas | emeas | unit | distance |             bibcode |      flux |     eflux |
    | SKYMAPPER.U | 11.37 | 0.01  |  mag |    0.209 | 2018PASA...35...10W | 2.502e-13 | 2.996e-15 |
    | SKYMAPPER.V | 11.34 | 0.00  |  mag |    0.209 | 2018PASA...35...10W | 2.110e-13 | 1.166e-15 |
    | SKYMAPPER.G | 10.64 | 0.00  |  mag |    0.209 | 2018PASA...35...10W | 2.427e-13 | 6.708e-16 |
    | SKYMAPPER.R | 10.32 | 0.01  |  mag |    0.209 | 2018PASA...35...10W | 2.203e-13 | 3.856e-15 |
    | SKYMAPPER.Z | 10.04 | 0.01  |  mag |    0.209 | 2018PASA...35...10W | 1.258e-13 | 2.203e-15 |
    |     GAIA2.G | 10.29 | 0.00  |  mag |    0.144 | 2018A&A...616A...1G | 1.896e-13 | 6.987e-17 |
    |    GAIA2.BP | 10.65 | 0.00  |  mag |    0.144 | 2018A&A...616A...1G | 2.219e-13 | 2.658e-16 |
    |    GAIA2.RP |  9.77 | 0.00  |  mag |    0.144 | 2018A&A...616A...1G | 1.619e-13 | 1.193e-16 |
    |     2MASS.J |  9.00 | 0.02  |  mag |    0.126 | 2003yCat.2246....0C | 7.621e-14 | 1.895e-15 |
    |     2MASS.H |  8.52 | 0.04  |  mag |    0.126 | 2003yCat.2246....0C | 4.474e-14 | 1.895e-15 |
    |    2MASS.KS |  8.40 | 0.01  |  mag |    0.126 | 2003yCat.2246....0C | 1.857e-14 | 3.250e-16 |

Where the columns have the following meaning:

- band: The photometric pass band of the measurement provided in SYSTEM.BAND format.
- meas: The measurement as how it is included in the catalog. Usually magnitudes, but can be fluxes as well.
- emeas: The error on the catalog measurement.
- unit: The unit of the catalog measurement.
- distance: The distance between the target and the catalog source in arcsec.
- bibcode: The bibcode of the article describing the source when available.
- flux: The flux of the measurement in erg/s/cm2/AA.
- eflux: The error on the flux in erg/s/cm2/AA.

The only collumns that are important for the fit, and the only ones that speedyfit will actually read are:

- band
- flux
- eflux

The names of these columns are the default values. You can change them, if you also update the corresponding names in
the setup file of the fit: photband_index, obs_index, err_index. See :doc:`setup_file`

You don't have to use photometry provided by speedyfit, you can use your own photometry as long as you have the
photometric bands in a format that speedyfit can read, the flux and the errors. You can have any other number of columns
in the photometry file, but they are not needed and are ignored by speedyfit.

.. note::

    The flux and its error have to be given in units of : erg/s/cm2/AA!
    The fitting part of Speedyfit will not preform any unit conversions.

Speedyfit allows two possible formats for the photometry file. The 'fixed_width' format described above, and a
headerless format that is read by
`astropy.io.ascii.write <https://docs.astropy.org/en/stable/api/astropy.io.ascii.write.html>`_ with the default format
value, and assuming no header (data_start=0, header_start=None). In the headerless format, the position of the
photometry band, flux and error needs to be specified as an integer.

An example of a headerless photometry file is:

.. code-block:: bash

    #    meas  e_meas flag unit   photband  source        cmeas       e_cmeas     cunit         bibcode            comments
    # float64 float64 |S20 |S30   |S30      |S50          float64     float64     |S50          |S19               |S74
      11.0835  0.0014 nan  mag    GAIA2.G   I/345/gaia2   9.37669e-14 9.37669e-16 erg/s/cm2/AA  2018yCat.1345....0G -
      11.0258  0.0055 nan  mag    GAIA2.BP  I/345/gaia2   1.61477e-13 1.61477e-15 erg/s/cm2/AA  2018yCat.1345....0G -
      11.0571   0.001 nan  mag    GAIA2.RP  I/345/gaia2   4.96983e-14 4.96983e-16 erg/s/cm2/AA  2018yCat.1345....0G -
       11.123   0.012 nan  mag    APASS.B   II/336/apass9  2.2233e-13 2.45728e-15 erg/s/cm2/AA  2015AAS...22533616H -
        11.17   0.009 nan  mag    APASS.V   II/336/apass9 1.23896e-13 1.23896e-15 erg/s/cm2/AA  2015AAS...22533616H -
       11.063   0.033 nan  ABmag  APASS.G   II/336/apass9 1.86108e-13 5.65658e-15 erg/s/cm2/AA  2015AAS...22533616H -
       11.251   0.013 nan  ABmag  APASS.R   II/336/apass9 8.73103e-14  1.0454e-15 erg/s/cm2/AA  2015AAS...22533616H -
       11.377   0.028 nan  ABmag  APASS.I   II/336/apass9 5.41327e-14 1.39603e-15 erg/s/cm2/AA  2015AAS...22533616H -

Where the columns of interest are the photband, cmeas and e_cmeas. To use this file you would specify:
photband_index = 4, obs_index = 6, err_index = 7 in the setup file.

.. note::

    Headerless does not necessarily mean without a header. As long as the header is commented it is ignored by
    speedyfit, and the file is considered headerless.

