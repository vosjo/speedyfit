
import numpy as np
import pylab as pl

import matplotlib.cm as cm

from scipy.stats import gaussian_kde

from ivs.sed import filters, model

import statfunc
 
 
def plot_distribution(data, parameters=None, percentiles=[16, 50, 84]):
   """
   Plots a histogram of the requested parameters, or all parameters if none are given.
   """
   
   if parameters is None:
      parameters = data.dtype.names
   
   x = len(parameters)
   
   for i, par in enumerate(parameters):
      d = data[par]
      
      pc  = np.percentile(d, percentiles, axis=0)
      
      pl.subplot(1,x,i+1)
      
      pl.hist(d, 50, normed=True)
      
      for l in pc:
         pl.axvline(x=l, ls='--', color='k')
         
         emin, emax = pc[1]-pc[0], pc[2]-pc[1]
         
         if np.min([emin, emax]) > 10:
            temp = "{}={:0.0f} -{:0.0f} +{:0.0f}"
         elif np.min([emin, emax]) > 1:
            temp = "{}={:0.1f} -{:0.1f} +{:0.1f}"
         else:
            temp = "{}={:0.2f} -{:0.2f} +{:0.2f}"
         
      pl.title(temp.format(par, pc[1], emin, emax))

def plot_distribution_2d(data, xpar, ypar, percentiles=[16, 50, 84]):
   
   x, y = data[xpar], data[ypar]
   
   xpc  = np.percentile(x, percentiles, axis=0)
   ypc = np.percentile(y, percentiles, axis=0)
   
   pl.hist2d(x, y, bins=50, normed=True)
   
   for i in xpc:
      pl.axvline(x=i, ls='--', color='w')
   for i in ypc:
      pl.axhline(y=i, ls='--', color='w')
   
   pl.xlabel(xpar)
   pl.ylabel(ypar)
   

def plot_distribution_density(data, xpar, ypar, percentiles=[16, 50, 84]):
   
   x, y = data[xpar][::10], data[ypar][::10]
   
   xpc  = np.percentile(x, percentiles, axis=0)
   ypc = np.percentile(y, percentiles, axis=0)
   
   # Calculate the point density
   xy = np.vstack([x,y])
   z = gaussian_kde(xy)(xy)

   # Sort the points by density, so that the densest points are plotted last
   idx = z.argsort()
   x, y, z = x[idx], y[idx], z[idx]
   
   pl.scatter(x, y, c=z, s=50, edgecolor='')
   
   
 
def plot_fit(obs, obs_err, photbands, pars={}, constraints={}):
   
   grid1 = dict(grid='kuruczsdb', z=0, Rv=3.1)
   grid2 = dict(grid='tmapsdb', z=0, Rv=3.1)
   model.set_defaults_multiple(grid1,grid2)
   
   
   colors = np.array([filters.is_color(p) for p in photbands])
   psystems = np.array([p.split('.')[0] for p in photbands])
   waves = np.array([filters.eff_wave(p) for p in photbands[~colors]])
   
   #-- def system colors
   all_systems = set(psystems)
   colorcycle = cm.rainbow(np.linspace(0, 1, len(all_systems)))
   system_colors = {}
   for p, c in zip(all_systems, colorcycle):
      system_colors[p] = c
   
   abs_xlim = (0.9 * np.min(waves), 1.1 * np.max(waves) )
   
   if not len(pars.keys()) == 0:
      ipars = pars.copy()
      for key, value in pars.items():
         ipars[key] = [value]
         
      syn, Labs = model.get_itable_pix(photbands=photbands, **ipars)
      syn = syn[:,0]
      
      ratio = (obs/syn)[~colors]
      weights = (obs/obs_err)[~colors]
      scale = np.average(ratio,weights=weights)
      
      
   
   # plot the fit of the absolute data
   #====================================
   ax = pl.subplot2grid((3,3), (0, 0), colspan=2, rowspan=2)
   
   
   if not len(pars.keys()) == 0:
      
      #-- synthetic integrated fluxes
      pl.plot(waves, scale*syn[~colors], 'xr')
      
      #-- synthetic model
      wave, flux = model.get_table(**pars)
      pl.plot(wave, scale*flux, '-r')
   
   for system in all_systems:
      s = np.where((psystems == system) & (~colors))
      if len(obs[s]) == 0: continue
   
      w = np.array([filters.eff_wave(p) for p in photbands[s]])
      pl.errorbar(w, obs[s], yerr=obs_err[s], ls='', marker='o', 
                  color=system_colors[system], label=system)
   
      
   pl.xlim(abs_xlim)
   pl.ylim([0.9*np.min(obs[~colors]), 1.1*np.max(obs[~colors])])
   pl.legend(loc='best', prop={'size': 9}, numpoints=1)
   ax.set_xscale("log", nonposx='clip')
   ax.set_yscale("log", nonposy='clip')
   
   pl.ylabel('Absolute Flux')
   
   #-- add band names to top of axis
   ax2 = ax.twiny()
   ax2.set_xscale("log", nonposx='clip')
   ax2.set_xlim(abs_xlim)
   waves = np.array([filters.eff_wave(p) for p in photbands[~colors]])
   bandnames = np.array([b.split('.')[-1] for b in photbands[~colors]])
   ax2.set_xticks(waves)
   ax2.set_xticklabels(bandnames)
   ax2.xaxis.tick_top()
   pl.minorticks_off()
   
   
   
   # plot O-C of the absolute measurements
   #====================================
   ax = pl.subplot2grid((3,3), (2, 0), colspan=2, rowspan=1)
   
   for system in all_systems:
      s = np.where((psystems == system) & (~colors))
      if len(obs[s]) == 0: continue
      
      w = np.array([filters.eff_wave(p) for p in photbands[s]])
      y = ((obs - syn*scale) / obs)[s]
      yerr = (obs_err / obs)[s]
      pl.errorbar(w, y , yerr=yerr, ls='', marker='o', color=system_colors[system])
   
   pl.axhline(y=0, color='k', ls='--')
   
   pl.xlim(abs_xlim)
   pl.gca().set_xscale('log',nonposx='clip')
   
   pl.ylabel('(O-C) / O')
   pl.xlabel('Wavelength (AA)')
   
   # plot O-C of the color fits
   #====================================
   ax = pl.subplot2grid((3,3), (0, 2), colspan=1, rowspan=2)
   
   y = np.array(range(len(obs[colors])))
   x = obs[colors] - syn[colors]
   x_err = obs_err[colors]
   oc_sys = psystems[colors]
   
   for system in all_systems:
      s = np.where(oc_sys == system)
      if len(x[s]) == 0: continue
   
      pl.errorbar(x[s], y[s], xerr=x_err[s], ls='', marker='o',
                  color=system_colors[system])
   
   pl.axvline(x=0, color='k', ls='--')
   
   yticknames = [b.split('.')[1] for b in photbands[colors]]
   pl.yticks(y, yticknames)
   pl.ylim([-y[-1]*0.1, y[-1]*1.1])
   
   ax.yaxis.tick_right()
   ax.xaxis.tick_top()
   ax.xaxis.set_label_position('top') 
   
   pl.xlabel('O-C')
   
   
   # plot fractional O-C of the constraints
   #=======================================
   ax = pl.subplot2grid((3,3), (2, 2), colspan=1, rowspan=1)
   
   if len(constraints.keys()) > 0:
      
      #-- get derived properties
      theta, pnames = [], []
      for n, v in pars.items():
         pnames.append(n)
         theta.append(v)
      
      derived_properties = statfunc.get_derived_properties(theta, pnames)
      derived_properties['distance'] = np.sqrt(1/scale) 
      
      x = range(len(constraints.keys()))
      y, yerr, xticknames = [], [], []
      for n, v in constraints.items():
         y.append( (v[0] - derived_properties[n] ) / v[0] )
         yerr.append(v[1]/ v[0])
         xticknames.append(n)
      
      pl.errorbar(x, y, yerr=yerr, ls='', marker='o')
   
   pl.axhline(y=0, color='k', ls='--')
   
   ax.yaxis.tick_right()
   ax.xaxis.tick_bottom()
   ax.yaxis.set_label_position('right') 
   
   pl.ylabel('(O-C) / O')
   pl.xlim([-x[-1]*0.1, x[-1]*1.1])
   pl.xticks(x, xticknames)
   