
import os 

import numpy as np
import pylab as pl

import matplotlib.cm as cm
import matplotlib.patches as patches
import matplotlib.ticker as ticker

from scipy.stats import gaussian_kde

from . import model, filters, reddening

from astropy.io import ascii

def format_parameter(name, value):
   
   temp = "{:7.3f}"
   if 'teff' in name or name == 'd':
      temp = "{:7.0f}"
   
   elif  'logg' in name or 'rad' in name or 'mass' in name or name=='L' or name == 'L2':
      temp = "{:7.2f}"
   
   elif 'ebv' in name or 'q' in name:
      temp = "{:7.3f}"
   
   if hasattr(value, '__iter__'):
      return [temp.format(v) for v in value]
   else:
      return temp.format(value)
   
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
   

def plot_constraints(constraints, samples, results):
   
   if 'distance' in constraints:
      constraints['d'] = np.array(constraints.pop('distance')) / 44365810.04823812
   
   pars = sorted(constraints.keys())
   
   for i, par in enumerate(pars):
   
      ax = pl.subplot(1, len(pars), i+1)
      
      pc = np.percentile(samples[par], [0.2, 16, 50, 84, 99.8])
      
      #-- plot 1 sigma range as box
      ax.add_patch(
         patches.Rectangle(
            (0.5, pc[1]),
            1.0,
            pc[3]-pc[1],
            fill=False 
         )
      )
         
      #-- plot best fit and 50 percentile fit
      pl.plot([0.5,1.5], [results[par][0], results[par][0]], '--r', lw=1.5)
      pl.plot([0.5,1.5], [results[par][1], results[par][1]], '-b', lw=1.5)
      
      #-- plot 3 sigma range as wiskers
      pl.plot([1.0, 1.0], [pc[0], pc[1]], '-k', lw=1.5, zorder=0)
      pl.plot([1.0, 1.0], [pc[3], pc[4]], '-k', lw=1.5, zorder=0)
      
      #pl.boxplot(samples[par], usermedians=usermedians)
      
      if par in constraints:
         pl.errorbar([1], constraints[par][0], 
                     yerr=[[constraints[par][1]],[constraints[par][2]]],
                     color='r', marker='x', mew=2, lw=2)
      
      ax.axes.get_xaxis().set_visible(False)
      
      if 'teff' in par:
         ticks = ticker.FuncFormatter(lambda x, pos: '{:0.1f}'.format(x/1000.))
         ax.yaxis.set_major_formatter(ticks)
      
      pl.title(par)

 
def plot_fit(obs, obs_err, photbands, pars={}, constraints={}, grids=[], gridnames=[], result='best', plot_components=True, plot_colors=False, plot_constraints=False):
   
   pars = pars.copy()
   
   # use model from 'best' results or 'pc' results
   resi = 0 if result == 'best' else 1

   print(gridnames)
   
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
   
   if not len(list(pars.keys())) == 0:
      ipars = pars.copy()
      for key, value in list(pars.items()):
         ipars[key] = [value[resi]]
         pars[key] = value[resi]
      _ = ipars.pop('d')
      _ = pars.pop('d')
      
      syn, Labs = model.get_itable(grid=grids, photbands=photbands, **ipars)
      syn = syn[:,0]
      
      scale = pars['scale']
   
   # plot the fit of the absolute data
   #====================================
   if not plot_colors and not plot_constraints:
      ax1 = pl.subplot2grid((3,3), (0, 0), colspan=3, rowspan=2)
   else:
      ax1 = pl.subplot2grid((3,3), (0, 0), colspan=2, rowspan=2)
   
   
   if not len(list(pars.keys())) == 0:
      
      #-- synthetic integrated fluxes
      pl.plot(waves, scale*syn[~colors], 'xr')
      
      #-- if possible, plot a non integrated model
      #   can only be done if provided gridnames are not paths to integrated files.
      if len(gridnames) > 0 and not os.path.isfile(gridnames[0]):
         #-- synthetic model
         ebv = pars['ebv']
         wave, flux = model.get_table(grid=gridnames, **pars)
         flux = reddening.redden(flux, wave=wave, ebv=ebv, rtype='flux', law='fitzpatrick2004')
         
         #ascii.write([wave, scale*flux], 'binary_model.txt', names=['wave', 'flux'], overwrite=True)
         
         pl.plot(wave, scale*flux, '-r')
         
         #-- plot components
         if plot_components and 'teff2' in pars:
            wave, flux = model.get_table(grid=gridnames[0], teff=pars['teff'], logg=pars['logg'], rad=pars['rad'], ebv=pars['ebv'])
            flux = reddening.redden(flux, wave=wave, ebv=ebv, rtype='flux', law='fitzpatrick2004')
            #ascii.write([wave, scale*flux], 'primary_model.txt', names=['wave', 'flux'], overwrite=True)
            pl.plot(wave, scale*flux, '--g')
            
            wave, flux = model.get_table(grid=gridnames[1], teff=pars['teff2'], logg=pars['logg2'], rad=pars['rad2'], ebv=pars['ebv'])
            flux = reddening.redden(flux, wave=wave, ebv=ebv, rtype='flux', law='fitzpatrick2004')
            #ascii.write([wave, scale*flux], 'secondary_model.txt', names=['wave', 'flux'], overwrite=True)
            pl.plot(wave, scale*flux, '--b')
         
         
      
   wave, flux, error, band = [], [], [], []
   for system in all_systems:
      s = np.where((psystems == system) & (~colors))
      if len(obs[s]) == 0: continue
   
      w = np.array([filters.eff_wave(p) for p in photbands[s]])
      pl.errorbar(w, obs[s], yerr=obs_err[s], ls='', marker='o', 
                  color=system_colors[system], label=system)
      
      wave.extend(w)
      flux.extend(obs[s])
      error.extend(obs_err[s])
      band.extend([system for i in w])
      
    
   ascii.write([wave, flux, error, band], 'observations.txt', names=['wave', 'flux', 'error', 'band'], overwrite=True)
   
      
   pl.xlim(abs_xlim)
   pl.ylim([0.9*np.min(obs[~colors]), 1.1*np.max(obs[~colors])])
   pl.legend(loc='best', prop={'size': 9}, numpoints=1)
   ax1.set_xscale("log", nonposx='clip')
   ax1.set_yscale("log", nonposy='clip')
   
   pl.ylabel('Absolute Flux')
   
   #-- add band names to top of axis
   ax2 = ax1.twiny()
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
   if not plot_colors and not plot_constraints:
      ax = pl.subplot2grid((3,3), (2, 0), colspan=3, rowspan=1)
   else:
      ax = pl.subplot2grid((3,3), (2, 0), colspan=2, rowspan=1)
   
   for system in all_systems:
      s = np.where((psystems == system) & (~colors))
      if len(obs[s]) == 0: continue
      
      w = np.array([filters.eff_wave(p) for p in photbands[s]])
      y = ((obs - syn*scale) / obs)[s]
      yerr = (obs_err / obs)[s]
      
      y = 2.5*np.log10(syn[s]*scale / obs[s])
      yerr = 2.5 / np.log(10.0) * obs_err[s] / obs[s]
      
      pl.errorbar(w, -y , yerr=yerr, ls='', marker='o', color=system_colors[system])
   
   pl.figtext(0.96, 0.96, '$\chi^2$ = {:0.2f}'.format(pars['chi2']), ha='right')
   
   pl.axhline(y=0, color='k', ls='--')
   
   pl.xlim(abs_xlim)
   pl.gca().set_xscale('log',nonposx='clip')
   
   pl.ylabel('(O-C) (mag)')
   pl.xlabel('Wavelength (AA)')
   

   # plot O-C of the color fits
   #====================================
   if len(obs[colors]) > 0 and plot_colors:
      
      ax = pl.subplot2grid((3,3), (0, 2), colspan=1, rowspan=2)
      
      y = np.array(list(range(len(obs[colors]))))
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
   if plot_constraints:
      ax = pl.subplot2grid((3,3), (2, 2), colspan=1, rowspan=1)
      
      x, xticknames = [], []
      
      if 'q' in constraints:
         
         q, ql, qu = constraints['q']
         
         ax.plot([1], pars['q'], 'k_', ms=15, mew=2)
         ax.errorbar([1], [q], yerr = [[ql], [qu]], marker='o')
         
         x.append(1)
         xticknames.append('q')
      
      
      if 'distance' in constraints:
         
         ax2 = ax.twinx()
         
         d  = constraints['distance'][0] / 44365810.04823812
         de = constraints['distance'][1] / 44365810.04823812
         
         
         ax2.plot([2], pars['d'], 'k_', ms=15, mew=2)
         ax2.errorbar([2], [d], yerr = de, marker='o')
         
         x.append(2)
         xticknames.append('distance')
      
      pl.xlim([0.5, 2.5])
      pl.xticks(x, xticknames)
   
   
   pl.figtext(0.03, 0.96, 'model = {}'.format(result))
