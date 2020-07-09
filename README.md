 
# Speedyfit

A python package to fit the photometric spectral energy distribution of stars. Uses a Markov chain Monte Carlo approach 
to determine the errors on the derived parameters.

Speedyfit is a command line tool writen in Python 3 that allows you to search the most common online databases for 
photometric observations of your target, and fit theoretical atmosphere models to the obtained photometry. Speedyfit can
deal with both single and binary stars, and allows for the inclusion of constraints from other sources, as for example
the distance or reddening. 

## Instalation

The simplest way to install Speedyfit is using pip from the terminal. This will install Speedyfit in the local folder.

    pip install git+https://github.com/vosjo/speedyfit.git#egg=speedyfit

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