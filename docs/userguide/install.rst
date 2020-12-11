 
Installation
============

The installation of speedyfit requires two steps, installing the python package, and downloading the required atmosphere
models. The speedyfit package can be installed with pip from the pypi repository as follows:

.. code-block:: bash

   pip install speedyfit

If you want install the current development version directly from github you can also use pip:

.. code-block:: bash

   pip install git+https://github.com/vosjo/speedyfit.git#egg=speedyfit

The atmosphere models that speedyfit uses to fit the SEDs can be downloaded from:

    http://www.astro.physik.uni-potsdam.de/~jorisvos/Speedyfit/modelgrids.tar.gz

Download them and unpack them in a directory of your choice. The last step is to store the path to the atmosphere models
in an environment variable so that Speedyfit will know where to get them. In a bash shell this is done as follows:

.. code-block:: bash

   export SPEEDYFIT_MODELS="<path to extracted atmosphere models>"

Where the path could be something like: '/home/user/speedyfit/modelgrids/'. To not always have to type this in the
terminal when running Speedyfit, you can add it to your '.bashrc' file.

To check that speedyfit can find all models run:

.. code-block:: bash

   speedyfit checkgrids

Which if everything went well should give you the following output:

.. code-block:: bash

   Checking which atmosphere models are available...
   Checking for models in <path to extracted atmosphere models>
   kurucz2
            raw: available
            integrated: available
   munari
            raw: available
            integrated: available
   tmap
            raw: available
            integrated: available
   blackbody
            raw: available
            integrated: available

If you get an error like:

.. code-block:: bash

   SPEEDYFIT_MODELS environmental variable not set. CAN NOT find models!
   Please point the SPEEDYFIT_MODELS variable to the directory where you stored models.
   On bash use:
   export SPEEDYFIT_MODELS='<path to extracted atmosphere models>'

Check that your "SPEEDYFIT_MODELS" variable is correctly set up. If you get "NOT FOUND" instead of "available" for any
of the models, check that the "SPEEDYFIT_MODELS" variable points to the correct directory, and that the directory
contains the extracted atmosphere models.

To uninstall Speedyfit, run:

.. code-block:: bash

   pip uninstall speedyfit