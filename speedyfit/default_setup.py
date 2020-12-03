
default_single = """
# photometry file with index to the columns containing the photbands, observations and errors
objectname: <objectname>
photometryfile: <photfilename>
photband_index: band
obs_index: flux
err_index: eflux
photband_exclude: <photband_exclude>
# parameters to fit and the limits on them in same order as parameters
pnames: [teff, logg, rad, ebv]
limits: <parameter_limits>
# constraints on distance and mass ratio is known
constraints: <constraints>
# added constraints on derived properties as mass, luminosity, luminosity ratio
derived_limits: {}
# path to the model grids with integrated photometry
grids: 
<model_grids>
# setup for the MCMC algorithm
nwalkers: 100    # total number of walkers
nsteps: 1000     # steps taken by each walker (not including burn-in)
nrelax: 250      # burn-in steps taken by each walker
a: 10            # relative size of the steps taken
# set the percentiles for the error determination 
percentiles: [0.2, 50, 99.8] # 16 - 84 corresponds to 1 sigma
# output options
resultfile: <objectname>_results_<postfix>.csv   # filepath to write results
plot1:
 type: sed_fit
 result: pc
 path: <objectname>_sed_<postfix>.png
plot2:
 type: distribution
 show_best: true
 path: <objectname>_distribution_<postfix>.png
 parameters: ['teff', 'rad', 'L', 'ebv', 'd', 'mass']
"""

default_binary = """
# photometry file with index to the columns containing the photbands, observations and errors
objectname: <objectname>
photometryfile: <photfilename>
photband_index: band
obs_index: flux
err_index: eflux
photband_exclude: ['GALEX', 'SDSS', 'WISE.W3', 'WISE.W4']
# parameters to fit and the limits on them in same order as parameters
pnames: [teff, logg, rad, teff2, logg2, rad2, ebv]
limits:
- [3500, 10000]
- [4.31, 4.31]
- [0.01, 2.5]
- [20000, 80000]
- [5.8, 5.8]
- [0.01, 0.5]
- [0, 0.10]
# constraints on distance and mass ratio if known
constraints: <constraints>
# added constraints on derived properties as mass, luminosity, luminosity ratio
derived_limits: {}
# path to the model grids with integrated photometry
grids: 
- kurucz
- tmap
# setup for the MCMC algorithm
nwalkers: 100    # total number of walkers
nsteps: 500     # steps taken by each walker (not including burn-in)
nrelax: 250      # burn-in steps taken by each walker
a: 10            # relative size of the steps taken
# set the percentiles for the error determination 
percentiles: [16, 50, 84] # 16 - 84 corresponds to 1 sigma
# output options
resultfile: <objectname>_results_<postfix>.csv   # filepath to write results
plot1:
 type: sed_fit
 path: <objectname>_sed_<postfix>.png
plot2:
 type: distribution
 path: <objectname>_distribution_<postfix>.png
 parameters: ['teff', 'rad', 'teff2', 'rad2', 'ebv', 'd']
"""


default_tmap = """
# photometry file with index to the columns containing the photbands, observations and errors
objectname: <objectname>
photometryfile: <photfilename>
photband_index: band
obs_index: flux
err_index: eflux
photband_exclude: ['GALEX', 'SDSS', 'WISE.W3', 'WISE.W4']
# parameters to fit and the limits on them in same order as parameters
pnames: [teff, logg, rad, ebv]
limits:
- [20000, 60000]
- [5.80, 5.80]
- [0.01, 0.5]
- [0, 0.10]
# constraints on distance and mass ratio is known
constraints: <constraints>
# added constraints on derived properties as mass, luminosity, luminosity ratio
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
percentiles: [0.2, 50, 99.8] # 16 - 84 corresponds to 1 sigma
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
"""