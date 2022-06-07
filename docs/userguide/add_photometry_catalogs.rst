Adding Photometry catalogs
==========================

Speedyfit is by default configured with a few standard catalogs that cover many different uses. If this is not
sufficient for you, you can activate other catalogs, or add your own. This is an easy process.

The catalog files are stored with the program data. Changing those is annoying, and an update would overwrite them.
Therefor, the first step is to copy them to your user space:

.. code-block:: bash

   speedyfit copy_catalogs

This will copy two photometry catalog files: tap_cats_phot.cfg and vizier_cats_phot.cfg to the SPEEDYFIT_MODELS
directory where you can edit them. Speedyfit will first check the SPEEDYFIT_MODELS directory, and if it can't find any
catalogs there, it will use the ones in its installation directory. So even when updating speedyfit, the catalogs in
your SPEEDYFIT_MODELS will not be affected.

You can now activate any of the other catalogs in those files by uncommenting them. Or you can add your own catalogs.

Vizier catalogs
---------------
The vizier catalogs are stored in the file: vizier_cats_phot.cfg. Most catalogs that you are going to use will be
available through `VIZIER <https://vizier.u-strasbg.fr/>`_. This is by far the easiest system to work with. Speedyfit
uses the `astroquery <https://astroquery.readthedocs.io/en/latest/>`_ library to query vizier.

The syntax to use is straightforward, see e.g. the APASS catalog:

.. code-block:: bash

    [II/336/apass9] # APASS all sky survey
    Bmag = APASS.B
    Vmag = APASS.V
    g_mag = APASS.G
    g_mag_unit = ABmag
    r_mag = APASS.R
    r_mag_unit = ABmag
    i_mag = APASS.I
    i_mag_unit = ABmag
    bibcode = 2015AAS...22533616H

in square brackets is the vizier name of the catalog: [II/336/apass9], and then for each included band there is a
mapping between the vizier name and the band name used by speedyfit e.g.: Bmag = APASS.B. You can also specify units
that are used to convert magnitudes to fluxes if necessary e.g.: g_mag_unit = ABmag. A bibcode can be added for
informational purposes: bibcode = 2015AAS...22533616H. For most catalogs you will only need to include the the band
mappings to have a functional system.

If the error keyword in the catalog is not the default, you can specify it with e.g.: <error_key_in_catalog> =
e_<band_name_in_catalog>. See e.g. the J/A+A/600/A50/sdcat catalog below for an example of that:

.. code-block:: bash

    [J/A+A/600/A50/sdcat] # subdwarf catalog Geier+2017: galex for sdB stars
    FUVmagc = GALEX.FUV
    FUVmagc_unit = ABmag
    FUVmagc_err = e_FUVmag
    NUVmagc = GALEX.NUV
    NUVmagc_unit = ABmag
    NUVmagc_err = e_NUVmag
    bibcode = 2017A&A...600A..50G

Here the error on the FUV band in the catalog is given by FUVmagc_err, instead of the default e_FUVmag.

TAP catalogs
------------
There is one other way of obtaining photometry, which is through the tap protocol. The Table Access Protocol specified
by the International Virtual Observatory Alliance: http://www.ivoa.net/documents/TAP/. It uses the Astronomical Data
Query Language which is based on SQL. This is harder to work with as there are many more possibilities for querying
data. Catalogs that use this format are stored in the tap_cats_phot.cfg file.

Currently there is only 1 catalog that uses the TAP protocol. You can try to add others, but there is no guarantee that
it will work.
