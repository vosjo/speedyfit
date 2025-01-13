Adding photometric bands to a grid
==================================

Speedyfit uses pre integrated atmosphere grid to make the fitting process efficient. If you want to include photometric
bands in your fit that are not included by default, you will need to create a new integrated atmosphere grid that
includes those bands. Currently it is not possible to add bands to an existing integrated grid, you will need to make a
new one with all the bands you want, and overwrite the old one.

File naming convention
----------------------

Original atmosphere models are named as follows:

    kurucz93_z0.0_k2odfnew_sed.fits

The integrated atmosphere models start with an 'i', and also include the reddening law (e.g. lawfitzpatrick2004) and the
reddening parameters used for this law (e.g. Rv3.10) in the name. The integrated values in these files are also already
de-reddened.

    ikurucz93_z0.0_k2odfnew_sed_lawfitzpatrick2004_Rv3.10.fits

You should have downloaded these files when you installed speedyfit.

Creating a new integrated grid
------------------------------

If your photometry band is included in speedyfit, you have everything you need to create a new integrated atmosphere
grid. It is not difficult to do this, but it can take a long time depending on how many bands you want to include,
what the reddening range is and the hardware you are running it on.

To started, create a new python script in whatever folder you want to work in. Add the following code to that script:

.. code-block:: python

    import numpy as np
    from speedyfit import integrate_grid

    # responses is a list of all the passbands that you want to include in the grid
    responses = ['GALEX', 'IUE', 'STROMGREN', 'JOHNSON', 'GAIA3E', 'GAIA2', 'SKYMAPPER', 'APASS', 'SDSS', '2MASS', 'WISE']

    # ebvs is an array of all the distinct reddening values you want to include in the grid. The smaller the step size,
    # the more accurate the result, but the longer the calculation will take.
    evbs = np.r_[0:1.02:0.05]

    # the gridname of the atmosphere models you want to integrate
    gridname = 'kurucz'

    # Create the new integrated grid
    integrate_grid.calc_integrated_grid(
        threads=6, # number of parallel threads to use when calculating the grid
        law='fitzpatrick2004', # the reddening law to use (best not to change this)
        Rv=3.1, # additional reddening parameter (best not to change this)
        ebvs=evbs,
        responses=responses,
        grid=gridname
    )


The run the script, and when it completes it will have created a new file in the current directory with the integrated
atomosphere models. In the case of the code above, the file will be
'ikurucz93_z0.0_k2odfnew_sed_lawfitzpatrick2004_Rv3.10.fits'. Expected output will be similar to this

.. code-block:: bash

    python make_new_grid.py
    ... subselection: GALEX
    ... subselection: IUE
    ... subselection: STROMGREN
    ... subselection: JOHNSON
    ... subselection: GAIA3E
    ... subselection: GAIA2
    ... subselection: SKYMAPPER
    ... subselection: APASS
    ... subselection: SDSS
    ... subselection: 2MASS
    ... subselection: WISE
    Selected response curves: GALEX.FUV, GALEX.NUV, IUE.B1, IUE.B2, IUE.B3, STROMGREN.B, STROMGREN.HBN, STROMGREN.HBW, STROMGREN.U, STROMGREN.V, STROMGREN.Y, JOHNSON.B, JOHNSON.H, JOHNSON.I, JOHNSON.J, JOHNSON.K, JOHNSON.L, JOHNSON.M, JOHNSON.N, JOHNSON.R, JOHNSON.U, JOHNSON.V, GAIA3E.BP, GAIA3E.G, GAIA3E.RP, GAIA2.BP, GAIA2.G, GAIA2.RP, SKYMAPPER.G, SKYMAPPER.I, SKYMAPPER.R, SKYMAPPER.U, SKYMAPPER.V, SKYMAPPER.Z, APASS.B, APASS.G, APASS.I, APASS.R, APASS.V, APASS2.B, APASS2.G, APASS2.I, APASS2.R, APASS2.V, SDSS.G, SDSS.GP, SDSS.I, SDSS.IP, SDSS.R, SDSS.RP, SDSS.U, SDSS.UP, SDSS.Z, SDSS.ZP, 2MASS.H, 2MASS.J, 2MASS.KS, WISE.W1, WISE.W2, WISE.W3, WISE.W4
    Total number of tables: 475
    3500.0 0.5 1 475: ET 2198 seconds

    ...

    49000.0 5.0 474 475: ET 0 seconds
    Written output to ikurucz93_z0.0_k2odfnew_sed_lawfitzpatrick2004_Rv3.10.fits
    Encountered 0 exceptions!



Add this file to the same directory where your other
atmosphere models are stored. And now you can use your new passband in your fits.
