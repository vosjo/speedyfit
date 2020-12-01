default_single = """
# photometry file with index to the columns containing the photbands, observations and errors
objectname: <objectname>
photometryfile: <photfilename>
photband_index: band
obs_index: flux
err_index: eflux
photband_exclude: []
# parameters to fit and the limits on them in same order as parameters
pnames: [teff, logg, rad, ebv]
limits:
- [20000, 60000]
- [5.80, 5.80]
- [0.01, 0.5]
- [0, 0.10]
# constraints on distance and mass ratio is known
constraints: 
  parallax: [<plx>, <e_plx>]
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

default_double = """
# photometry file with index to the columns containing the photbands, observations and errors
objectname: <objectname>
photometryfile: <photfilename>
photband_index: band
obs_index: flux
err_index: eflux
photband_exclude: []
# parameters to fit and the limits on them in same order as parameters
pnames: [teff, logg, rad, teff2, logg2, rad2, ebv]
limits:
- [3500, 10000]
- [4.31, 4.31]
- [0.01, 2.5]
- [20000, 50000]
- [5.8, 5.8]
- [0.01, 0.5]
- [0, 0.10]
# constraints on distance and mass ratio if known
constraints:
  parallax: [<plx>, <e_plx>]
# added constraints on derived properties as mass, luminosity, luminosity ratio
derived_limits: {}
# path to the model grids with integrated photometry
grids: 
- kurucz2
- tmap
# setup for the MCMC algorithm
nwalkers: 100    # total number of walkers
nsteps: 500     # steps taken by each walker (not including burn-in)
nrelax: 250      # burn-in steps taken by each walker
a: 10            # relative size of the steps taken
# set the percentiles for the error determination 
percentiles: [16, 50, 84] # 16 - 84 corresponds to 1 sigma
# output options
resultfile: <objectname>_results_binary.csv   # filepath to write results
plot1:
 type: sed_fit
 path: <objectname>_sed_binary.png
plot2:
 type: distribution
 path: <objectname>_distribution_primary.png
 parameters: ['teff', 'rad', 'teff2', 'rad2', 'ebv', 'd']
#plot3:
# type: distribution
# #path: <objectname>_distribution_derived.png
# parameters: ['mass', 'L', 'mass2', 'L2', 'q']
"""