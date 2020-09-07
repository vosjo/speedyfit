 
# Speedyfit

A python package to fit the photometric spectral energy distribution of stars. Uses a Markov chain Monte Carlo approach 
to determine the errors on the derived parameters.

Speedyfit is a command line tool writen in Python 3 that allows you to search the most common online databases for 
photometric observations of your target, and fit theoretical atmosphere models to the obtained photometry. Speedyfit can
deal with both single and binary stars, and allows for the inclusion of constraints from other sources, as for example
the distance or reddening. 

## Installation

The installation of speedyfit requires two steps, installing the python package, and downloading the required atmosphere
models. The speedyfit package can be installed with pip from the pypi repository as follows:

    pip install speedyfit
    
The atmosphere models that speedyfit uses to fit the SEDs can be downloaded from:

    http://www.astro.physik.uni-potsdam.de/~jorisvos/Speedyfit/modelgrids.tar.gz
    
Download them and unpack them in a directory of your choice. The last step is to store the path to the atmosphere models
in an environment variable so that Speedyfit will know where to get them. In a bash shell this is done as follows:

    export SPEEDYFIT_MODELS="<path to extracted atmosphere models>"
    
Where the path could be something like: '/home/user/speedyfit/modelgrids/'. To check that speedyfit can find all models 
run:

    python -c "from speedyfit.model import check_grids; check_grids()"

Which if everything went well should give you the following output:

    Checking which atmosphere models are available...
    kurucz2
             raw: available
             integrated: available
    munari
             raw: available
             integrated: available
    tmap
             raw: available
             integrated: available
    koester
             raw: available
             integrated: available
    blackbody
             raw: available
             integrated: available

If you get "NOT FOUND" for any of the models, check that the "SPEEDYFIT_MODELS" variable is correctly set up.

To uninstall Speedyfit, run:

    pip uninstall speedyfit
    
## Searching for photometry

Speedyfit can automatically download photometry from several of the large surveys, and many smaller studies available 
on Vizier. By default only photometry from the following large and well calibrated surveys is obtained:

- Gaia DR2: https://gea.esac.esa.int/archive/, [Gaia collaboration et al. 2018, A&A](https://ui.adsabs.harvard.edu/abs/2018A%26A...616A...1G/abstract)
- Sky Mapper DR1: http://skymapper.anu.edu.au/, [Wolf et al. 2018, PASA](https://ui.adsabs.harvard.edu/abs/2018PASA...35...10W)
- APASS DR9: https://www.aavso.org/apass, [Henden et al. 2015, AAS](https://ui.adsabs.harvard.edu/abs/2015AAS...22533616H)
- 2MASS: https://irsa.ipac.caltech.edu/Missions/2mass.html, [Skrutskie et al. 2006, AJ](https://ui.adsabs.harvard.edu/abs/2006AJ....131.1163S/abstract)
- WISE: http://wise.ssl.berkeley.edu/, [Cutri et al. 2012, yCat](http://cdsads.u-strasbg.fr/abs/2012yCat.2311....0C)

You can download photometry for a given object as follows:

    speedyfit <object_name> -empty single --phot
    
This will create a setup file called object_name_single.yaml where you can setup the fitting parameters and a photometry
file called object_name.phot with the obtained photometry. The object name should be resolvable by simbad, or if that 
is not possible, a J-type coordinate. For example, to obtain photometry of the star 'EO Ceti' the following object names
will work: 'EO Ceti', 'J012343.24-05 05 45.83' and 'J020.93019703-05.096064175' and any other alias that is recoginzed
by simbad.

## Fitting the SED

The the -empty command, a default setup file is created where you can define the parameters of your fit. A default yaml
file for a single star fit looks as follows. 

```yaml
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
```

In that file you can setup the photometry, the parameters and ranges to include in the fit, any constraints that are 
known, the setup of the MCMC algorithm and what output you like to have.

to run a fit:

    speedyfit <object_name_setupfile.yaml>
    
## constraints

Speedyfit allows you to add constraints on fitted parameters as the distance or redenning from other sources. As well as
constraints on derived parameters as the mass or luminosity if those are known. The constraints with their errors will 
be propagated in the fit, and are included in the log likelyhood function in the MCMC algorithm. If photometry was 
obtained with speedyfit, the gaia distance will have been obtained automatically as well. 

A typical use of constraints on fitted parameters as effective temperature and surface gravity is to use the SED to
only determine the radius of the star is the atmospheric parameters are already obtained from for example as 
spectroscopic analysis.

## Example usage

Lets try to fit the sdB+G type binary PG1104+243. We start with obtaining the photometry:

    speedyfit PG1104+243 -empty binary --phot
    
The '-empty binary' option creates a default setup file for a binary fit, and the '--phot' option downloads photometry 
from the standard sources. This command will have created 2 files: 'PG1104+243_binary.yaml' and 'PG1104+243.phot'. Lets 
start with looking at the photometry file:

```
|      band |               meas |                 emeas |  unit |           distance |             bibcode |                   flux |                  eflux |
|   GAIA2.G |            11.2007 |                0.0012 |   mag |               1.05 | 2018A&A...616A...1G |  8.248395510256044e-14 |    9.1164636206566e-17 |
|  GAIA2.BP |            11.2266 |                0.0106 |   mag |               1.05 | 2018A&A...616A...1G | 1.3151989765830283e-13 |  1.284023604507669e-15 |
|  GAIA2.RP |            11.0464 |                0.0018 |   mag |               1.05 | 2018A&A...616A...1G | 5.0236794630799464e-14 |  8.328563599441116e-17 |
|   APASS.B | 11.359999656677246 |   0.04699999839067459 |   mag | 0.7320000000000001 | 2015AAS...22533616H |  1.785653711432151e-13 |   7.72984461564579e-15 |
|   APASS.V | 11.331999778747559 |  0.017999999225139618 |   mag | 0.7320000000000001 | 2015AAS...22533616H | 1.0751193682636523e-13 | 1.7823986812698607e-15 |
|   APASS.G | 11.210000038146973 |    0.0989999994635582 | ABmag | 0.7320000000000001 | 2015AAS...22533616H |  1.613585138824778e-13 | 1.4713051584537062e-14 |
|   APASS.R | 11.310999870300293 |  0.014000000432133675 | ABmag | 0.7320000000000001 | 2015AAS...22533616H |  8.454580460944508e-14 | 1.0901739261158757e-15 |
|   APASS.I | 11.350000381469727 |  0.027000000700354576 | ABmag | 0.7320000000000001 | 2015AAS...22533616H |  5.398668203455588e-14 | 1.3425364709722915e-15 |
|   2MASS.J | 10.767999649047852 |  0.026000000536441803 |   mag |              0.132 |                   - | 1.5039297610503095e-14 |  3.601443372959058e-16 |
|   2MASS.H | 10.520000457763672 |  0.027000000700354576 |   mag |              0.132 |                   - |  7.091360781455302e-15 | 1.7634739011803824e-16 |
|  2MASS.KS | 10.510000228881836 |  0.023000000044703484 |   mag |              0.132 |                   - |  2.674932088409217e-15 |  5.666518062432707e-17 |
|    SDSS.U |  11.87399959564209 | 0.0020000000949949026 |   mag |              0.354 |                   - | 1.6258783227745728e-13 | 2.9949786934881625e-16 |
|    SDSS.G | 11.284000396728516 | 0.0010000000474974513 |   mag |              0.354 |                   - | 1.5247601537219582e-13 | 1.4043560668439271e-16 |
|    SDSS.R | 11.440999984741211 | 0.0010000000474974513 |   mag |              0.354 |                   - |  7.587446336921016e-14 |  6.988296663640909e-17 |
|    SDSS.I | 11.449999809265137 | 0.0010000000474974513 |   mag |              0.354 |                   - |  5.095632528557113e-14 |  4.693251222769927e-17 |
|    SDSS.Z | 12.795999526977539 |  0.004999999888241291 |   mag |              0.354 |                   - | 1.0300768439060776e-14 |  4.743679064803486e-17 |
|   WISE.W1 | 10.456000328063965 |  0.023000000044703484 |   mag |              0.636 |                   - |  5.373832766226388e-16 | 1.1383810664301285e-17 |
|   WISE.W2 | 10.482999801635742 |  0.020999999716877937 |   mag |              0.636 |                   - | 1.5478084540820753e-16 | 2.9937269251093825e-18 |
|   WISE.W3 | 10.451000213623047 |   0.07800000160932541 |   mag |              0.636 |                   - | 4.3005202859946986e-18 | 3.0895220013709453e-19 |
|   WISE.W4 |  9.104000091552734 |                   nan |   mag |              0.636 |                   - | 1.1617863723639043e-18 | 2.1400895857990046e-20 |

```

We likely want to remove the SDSS photometry as it is not very reliable at the bright end, and the WISE.W4 band as it 
doesn't have an error. The WISE.W3 band does have an error stated but is very close to the detection limit, and likely
not reliable.

In the 'PG1104+243_binary.yaml' we don't have to change anything as the default settings are good for this object, and 
the parallax was filled automatically. The defaults assume a cool companion for which Speedyfit will use the Kurucz 
 model grid, and a hot component for which the TMAP grid is used. Also by default the logg is not fit, but fixed at 
 reasonable values as an SED fit seldom can constrain logg.
 
 Lets fit the SED:
 
    speedyfit PG1104+243_binary.yaml
    
 This will output:
 
```
Applied constraints: 
         distance = 274.02515550927575 - 6.17988937560995 + 6.17988937560995
100%|██████████████████████████████████████████████████████████████████████| 750/750 [01:12<00:00, 10.29it/s]
================================================================================

Resulting parameter values and errors:
   Par             Best        Pc       emin       emax
   teff       =    6201      6134   -     74   +     60
   rad        =    0.89      0.90   -   0.02   +   0.02
   teff2      =   47542     37532   -   6553   +   5984
   rad2       =    0.12      0.13   -   0.01   +   0.03
   ebv        =   0.090     0.068   -  0.033   +  0.023
   mass       =    0.59      0.60   -   0.03   +   0.03
   mass2      =    0.31      0.38   -   0.08   +   0.21
   q          =   1.883     1.594   -  0.572   +  0.431
   lr         =   0.040     0.050   -  0.007   +  0.010
   rr         =   7.628     7.018   -  1.400   +  0.891
   d          =     273       275   -      6   +      7
   L          =    1.05      1.03   -   0.06   +   0.06
   L2         =   26.49     20.62   -   4.00   +   3.86
   scale      =   0.000     0.000   -  0.000   +  0.000
   chi2       =   8.411    13.424   -  2.532   +  4.364

```
    
 And will produce a figure of the SED with the best fitting model, and a distribution of the parameters:
 
 ![example SED fit](https://raw.githubusercontent.com/vosjo/speedyfit/master/docs/PG1104%2B243_sed_binary.png)
 
 ![example parameter distribution]( https://raw.githubusercontent.com/vosjo/speedyfit/master/docs/PG1104%2B243_distribution_primary.png)
