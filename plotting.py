
import numpy as np
import pylab as pl

from ivs.sed import filters, model

import statfunc
 
def plot_fit(obs, obs_err, photbands, pars={}, constraints={}):
   
   grid1 = dict(grid='kuruczsdb', z=0, Rv=3.1)
   grid2 = dict(grid='tmapsdb', z=0, Rv=3.1)
   model.set_defaults_multiple(grid1,grid2)
   
   
   colors = np.array([filters.is_color(p) for p in photbands])
   waves = np.array([filters.eff_wave(p) for p in photbands[~colors]])
   
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
   
   pl.errorbar(waves, obs[~colors], yerr=obs_err[~colors], ls='', marker='o')
   
   if not len(pars.keys()) == 0:
      
      #-- synthetic integrated fluxes
      pl.plot(waves, scale*syn[~colors], 'xr')
      
      #-- synthetic model
      wave, flux = model.get_table(**pars)
      pl.plot(wave, scale*flux, '-r')
      
   pl.xlim(abs_xlim)
   pl.ylim([0.9*np.min(obs[~colors]), 1.1*np.max(obs[~colors])])
   ax.set_xscale("log", nonposx='clip')
   ax.set_yscale("log", nonposy='clip')
   
   
   # plot O-C of the absolute measurements
   #====================================
   ax = pl.subplot2grid((3,3), (2, 0), colspan=2, rowspan=1)
   
   y = ((obs - syn*scale) / obs)[~colors]
   yerr = (obs_err / obs)[~colors]
   pl.errorbar(waves, y , yerr=yerr, ls='', marker='o')
   
   pl.axhline(y=0, color='k', ls='--')
   
   pl.xlim(abs_xlim)
   pl.gca().set_xscale('log',nonposx='clip')
   
   pl.ylabel('(O-C) / O')
   
   # plot O-C of the color fits
   #====================================
   ax = pl.subplot2grid((3,3), (0, 2), colspan=1, rowspan=2)
   
   y = range(len(obs[colors]))
   
   pl.errorbar(obs[colors] - syn[colors], y, xerr=obs_err[colors], ls='', marker='o')
   
   #pl.plot(syn[colors], y, 'xr')
   
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
   