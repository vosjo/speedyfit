Photometry
==========

Speedyfit can automatically obtain photometry from several large online databases, using both Vizier and TAP queries.
It will use simbad to find the coordinates for your system, and will do a search based on those coordinates. The
included catalogs are given below:

Photometry Catalogs
-------------------

Gaia DR2 (Optical)
^^^^^^^^^^^^^^^^^^

The second data release of the Gaia Spacecraft:  https://gea.esac.esa.int/archive/
Included photometric bands are:

- GAIA2.G
- GAIA2.BP
- GAIA2.RP

This catalog is queried through Vizier: `I/345/gaia2 <https://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/345/gaia2>`_

Cite as: `Gaia collaboration et al. 2018, A&A <https://ui.adsabs.harvard.edu/abs/2018A%26A...616A...1G/abstract>`_

Skymapper DR1 (Optical)
^^^^^^^^^^^^^^^^^^^^^^^

The first public data release of the Skymapper survey: http://skymapper.anu.edu.au/
Included photometric bands are:

- SKYMAPPER.U
- SKYMAPPER.V
- SKYMAPPER.G
- SKYMAPPER.R
- SKYMAPPER.I
- SKYMAPPER.Z

This catalog is queried through TAP at: http://api.skymapper.nci.org.au/public/tap/

Cite as: `Wolf et al. 2018, PASA <https://ui.adsabs.harvard.edu/abs/2018PASA...35...10W>`_

APASS (Optical)
^^^^^^^^^^^^^^^

The 9th data release of The AAVSO Photometric All-Sky Survey: https://www.aavso.org/apass. The APASS survey is valid
from roughly 7th to 17th magnitude. Included photometric bands are:

- JOHNSON.B
- JOHNSON.V
- SLOAN.U (only available for very few stars)
- SLOAN.G
- SLOAN.R
- SLOAN.I
- SLOAN.Z

Even though these bands are supposed to be Johnson and Sloan, they differ slightly from the standard bands, and
speedyfit uses specific APASS passbands.

This catalog is queried through Vizier: `II/336/apass9 <https://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=II/336/apass9>`_

Cite as: `Henden et al. 2015, AAS <https://ui.adsabs.harvard.edu/abs/2015AAS...22533616H>`_

Note: The 10th datarelease of APASS is public, but can not be automatically queried. You can include this data by hand by obtaining it from:
https://www.aavso.org/apass-dr10-download

2MASS (IR)
^^^^^^^^^^

The Two Micron All Sky Survey: https://irsa.ipac.caltech.edu/Missions/2mass.html. Included photometric bands are:

- 2MASS.J (Magnitude limit = 15.8)
- 2MASS.H (Magnitude limit = 15.1)
- 2MASS.KS (Magnitude limit = 14.3)

This catalog is queried through Vizier: `II/246/out <https://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=II/246/out>`_

Cite as: `Skrutskie et al. 2006, AJ <https://ui.adsabs.harvard.edu/abs/2006AJ....131.1163S>`_

WISE (IR)
^^^^^^^^^

The Wide-field Infrared Survey Explorer: http://wise.ssl.berkeley.edu/. Included photometric bands are:

- WISE.W1
- WISE.W2
- WISE.W3
- WISE.W4

This catalog is queried through Vizier: `II/311/wise <https://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=II/311/wise>`_

Cite as: `Cutri et al. 2012, yCat <http://cdsads.u-strasbg.fr/abs/2012yCat.2311....0C>`_

Downloading photometry
----------------------

Downloading photometry for a target can be done through the command line by adding the '--phot' option when using the
'empty' function of speedyfit. Or by using the 'get_photometry' function.

Command line use
^^^^^^^^^^^^^^^^

.. code-block:: bash

   speedyfit <object_name> -empty single --phot

Will download all photometry that speedyfit can find from the above catalogs and store it in a file called
'object_name.phot'