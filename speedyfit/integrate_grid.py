import os
import sys
import time
import shutil

import numpy as np

from astropy.io import fits
from astropy.table import Table

from multiprocessing import cpu_count, Manager, Process

from speedyfit import model, reddening, filters


def get_responses(responses=None, wave=(0, np.inf)):
    """
    Get a list of response functions and their information.

    You can specify bandpass systems (e.g. GENEVA) and then all Geveva filters
    will be collected. You can also specific specific bandpasses (e.g. GENEVA.V),
    in which case only that one will be returned. The default C{None} returns
    all systems except some Hubble ones.

    You can also set responses to C{None} and give a wavelength range for filter
    selection.

    Example input for C{responses} are:

    >>> responses = ['GENEVA.V','2MASS.J']
    >>> responses = ['GENEVA','2MASS.J']
    >>> responses = ['BOXCAR','STROMGREN']
    >>> responses = None

    @param responses: a list of filter systems of passbands
    @type responses: list of str
    @param add_spectrophotometry: add spectrophotometric filters to the filter list
    @type add_spectrophotometry: bool
    @param wave: wavelength range
    @type wave: tuple (float,float)
    """

    # -- if no responses are given, select using wavelength range
    if responses is None:
        responses = filters.list_response(wave_range=(wave[0], wave[-1]))
    else:
        responses_ = []
        for resp in responses:
            print('... subselection: {}'.format(resp))
            responses_ += filters.list_response(resp)
        responses = responses_
    # -- get information on the responses
    responses = [resp for resp in responses if not (
            ('ACS' in resp) or ('WFPC' in resp) or ('STIS' in resp) or ('ISOCAM' in resp) or ('NICMOS' in resp))]

    print('Selected response curves: {}'.format(', '.join(responses)))

    return responses


def get_threads(threads, max=np.inf):
    """
    Reads the threadcount, and returns an integer
    accepts: <integer>, 'max', 'half', 'safe'
       max: all cpus are used
       half: half of the cpus are used
       safe: all but one of the cpus are used
    """
    if threads == 'max':
        threads = cpu_count()
    elif threads == 'half':
        threads = cpu_count() / 2
    elif threads == 'safe':
        threads = cpu_count() - 1
    threads = int(threads)

    if threads > max:
        threads = max

    return threads


# def calc_integrated_grid(threads=1, ebvs=None, law='fitzpatrick2004', Rv=3.1, responses=None, update=False, grid=None):
#
#     # get the teff and logg combinations of the grid
#     teffs, loggs = model.get_grid_dimensions(grid=grid)
#
#     def process_ebvs(wave, flux, ebvs):
#
#         arr = []
#
#         for ebv in ebvs:
#             # redden the flux
#             flux_ = reddening.redden(flux, wave=wave, ebv=ebv, rtype='flux', law=law, Rv=Rv)
#
#             # calculate synthetic fluxes
#             synflux = model.synthetic_flux(wave, flux_, responses)
#
#             arr.append([np.concatenate(([ebv], synflux))])
#
#         return arr
#
#     for teff, logg in zip(teffs, loggs):
#
#         # get the spectrum
#         wave, flux = model.get_table(grid=grid, teff=teff, logg=logg)
#
#         results = process_ebvs(wave, flux, ebvs)

        

def calc_integrated_grid(threads=1, ebvs=np.r_[0:2.01:0.01], law='fitzpatrick2004', Rv=3.1, responses=None, grid=None):
    """
    Integrate an entire SED grid over all passbands and save to a FITS file.

    The output file can be used to fit SEDs more efficiently, since integration
    over the passbands has already been carried out.

    WARNING: this function can take a long time to compute!

    Extra keywords can be used to specify the grid.

    :param threads: number of threads
    :type threads; integer, 'max', 'half' or 'safe'
    :param ebvs: reddening parameters to include
    :type ebvs: numpy array
    :param law: interstellar reddening law to use
    :type law: string (valid law name, see C{reddening.py})
    :param Rv: Rv value for reddening law
    :type Rv: float
    :param responses: respons curves to add (if None, add all)
    :type responses: list of strings
    """

    # -- select number of threads
    threads = get_threads(threads, max=len(ebvs))

    # -- Get the grid, extract the different teff and logg combinations provided and load the first model to check
    #    which photbands can be calculated.
    gridfile = model.get_grid_file(integrated=False, grid=grid)
    ff = fits.open(gridfile)
    data = ff[1].data
    wave, flux = data['wavelength'], data['flux']

    teffs = []
    loggs = []
    hdus = []
    for hdu in ff[1:]:
        teffs.append(float(hdu.header['TEFF']))
        loggs.append(float(hdu.header['LOGG']))
        hdus.append(hdu)

    responses = get_responses(responses=responses, wave=wave)

    # -- definition of one process for multi processing:
    def do_ebv_process(ebvs, arr, responses):
        # Run over all reddening values, calculate the reddened model and integrate it over all photbands.
        for ebv in ebvs:
            # redden the model
            flux_ = reddening.redden(flux, wave=wave, ebv=ebv, rtype='flux', law=law, Rv=Rv)

            # calculate synthetic fluxes
            synflux = filters.synthetic_flux(wave, flux_, responses)

            # append to results
            arr.append([np.concatenate(([ebv], synflux))])

    # -- prepare the array containing the integrated fluxes
    #   (1 row per model, 1 column for each response curve and teff, logg, ebv and total luminosity)
    output = np.zeros((len(teffs) * len(ebvs), 4 + len(responses)))
    start = 0

    # also track possible errors
    exceptions = 0
    exceptions_logs = []

    print('Total number of tables: %i ' % (len(teffs)))
    c0 = time.time()
    for i, (teff, logg, hdu) in enumerate(zip(teffs, loggs, hdus)):
        if i > 0:
            print('%s %s %s %s: ET %d seconds' % (teff, logg, i, len(teffs), (time.time() - c0) / i * (len(teffs) - i)))

        # -- get model SED and absolute luminosity
        wave, flux = hdu.data['wavelength'], hdu.data['flux']
        Labs = model.luminosity(wave, flux)

        # -- threaded calculation over all E(B-V)s
        manager = Manager()
        arr = manager.list([])
        all_processes = []
        for j in range(threads):
            all_processes.append(Process(target=do_ebv_process, args=(ebvs[j::threads], arr, responses)))
            all_processes[-1].start()
        for p in all_processes:
            p.join()

        try:
            # -- collect the results and add them to 'output'
            arr = np.vstack([row for row in arr])
            sa = np.argsort(arr[:, 0]) # sort on ebv
            arr = arr[sa]
            output[start:start + arr.shape[0], :3] = teff, logg, Labs
            output[start:start + arr.shape[0], 3:] = arr
            start += arr.shape[0]
        except:
            print('Exception in calculating Teff=%f, logg=%f' % (teff, logg))
            print('Exception: %s' % (sys.exc_info()[1]))
            exceptions = exceptions + 1
            exceptions_logs.append(sys.exc_info()[1])

    # -- close the fits file with the model
    ff.close()

    # -- create Table object and add header info to the table
    output = Table(data=output, names=['teff', 'logg', 'Labs', 'ebv'] + responses)
    output.meta['gridfile'] = (os.path.basename(gridfile), 'original model file')
    output.meta['GRID'] = (grid, 'name of the model grid')
    output.meta['fluxtype'] = ('Flambda', 'units of the flux')
    output.meta['redlaw'] = (law, 'interstellar reddening law')
    output.meta['rv'] = (Rv, 'interstellar reddening parameter')

    # -- create the name of the output file and safe to disk
    outfile = 'i{0}'.format(os.path.basename(gridfile))
    outfile = os.path.splitext(outfile)
    outfile = outfile[0] + '_law{0}_Rv{1:.2f}'.format(law, Rv) + outfile[1]
    if os.path.isfile(outfile):
        print('Precaution: making original grid backup at {0}.backup'.format(outfile))
        shutil.copy(outfile, outfile + '.backup')
    output.write(outfile, overwrite=True)

    # # -- make FITS columns
    # output = output.T
    # cols = [fits.Column(name='teff', format='E', array=output[0]),
    #         fits.Column(name='logg', format='E', array=output[1]),
    #         fits.Column(name='ebv', format='E', array=output[3]),
    #         fits.Column(name='Labs', format='E', array=output[2])]
    # for i, photband in enumerate(responses):
    #     cols.append(fits.Column(name=photband, format='E', array=output[4 + i]))
    #
    # # -- make FITS extension and write grid/reddening specifications to header
    # table = fits.TableHDU.from_columns(fits.ColDefs(cols))
    # table.header.update(gridfile=os.path.basename(gridfile))
    # table.header.update(GRID=(grid, 'name of the model grid'),
    #                     FLUXTYPE=('Flambda', 'units of the flux'),
    #                     REDLAW=(law, 'interstellar reddening law'),
    #                     RV=(Rv, 'interstellar reddening parameter'))
    # # -- make/update complete FITS file
    # if os.path.isfile(outfile):
    #     os.remove(outfile)
    #     print('Removed existing file: %s' % (outfile))

    # hdulist = fits.HDUList([])
    # hdulist.append(fits.PrimaryHDU(np.array([[0, 0]])))
    # hdulist.append(table)
    # hdulist.writeto(outfile)
    print("Written output to %s" % outfile)

    print('Encountered %s exceptions!' % exceptions)
    for i in exceptions_logs:
        print('ERROR\n', i)

    return outfile


def check_grid(grid):
    """
    Check if the grid with integrated photometry calculated with calc_integrated_grid
    has all header information, and if all models succeeded.
    """
    print('Checking grid: {}'.format(grid))

    hdulist = fits.open(grid, mode='update')
    names = hdulist[1].columns.names
    for i, name in enumerate(names):
        if name.lower() in ['teff', 'logg', 'ebv', 'labs', 'vrad', 'rv', 'z']:
            names[i] = name.lower()
    cols = [fits.Column(name=name, format='E', array=hdulist[1].data.field(name)) for name in names]
    N = len(hdulist[1].data)

    keys = [key.lower() for key in hdulist[1].header.keys()]

    if 'z' not in hdulist[1].header:
        hdulist[1].header['Z'] = 0.0
        print('Adding metallicity (Z={}) to header!'.format(hdulist[1].header['Z']))

    if 'Rv' not in hdulist[1].header:
        hdulist[1].header['Rv'] = 3.1
        print('Adding Rv (Rv={}) to header!'.format(hdulist[1].header['Rv']))

    if 'z' not in names:
        z = hdulist[1].header.get('z', 0.0)
        print('Adding metallicity from header {}'.format(z))
        cols.append(fits.Column(name='z', format='E', array=np.ones(N) * z))
    else:
        print("Metallicity already in there")
    if 'vrad' not in names:
        vrad = 0.
        print('Adding radial velocity {}'.format(vrad))
        cols.append(fits.Column(name='vrad', format='E', array=np.ones(N) * vrad))
    else:
        print("Radial velocity already in there")

    fix_rv = False
    if 'rv' not in names:
        if 'rv' in keys:
            rv = hdulist[1].header['Rv']
            print("Adding interstellar Rv from header {}".format(rv))
        else:
            rv = 3.1
            print("Adding default interstellar Rv {}".format(rv))
        cols.append(fits.Column(name='rv', format='E', array=np.ones(N) * rv))
    elif not hdulist[1].header['Rv'] == hdulist[1].data.field('rv')[0]:
        rv = hdulist[1].header['Rv']
        fix_rv = rv
        print('Correcting interstellar Rv with {}'.format(rv))
    else:
        print("Interstellar Rv already in there")

    table = fits.TableHDU.from_columns(fits.ColDefs(cols))
    if fix_rv:
        table.data.field('rv')[:] = rv
    fake_keys = [key.lower() for key in table.header.keys()]

    # make sure all keywords set above are also in the new table header
    table.header.update(hdulist[1].header)

    # now update header keywords related to the created table
    for key in hdulist[1].header.keys():
        if not key.lower() in fake_keys:
            if len(key) > 8:
                key = 'HIERARCH ' + key
            table.header.update(key=hdulist[1].header[key])
    hdulist[1] = table

    logstring = "Axis:\n"
    for name in hdulist[1].columns.names:
        if name.islower() and not name == 'labs':
            ax = np.unique(hdulist[1].data.field(name))
            logstring += "\t\t{} {} {} {}\n".format(name, len(ax), min(ax), max(ax))
    print(logstring)

    keep = hdulist[1].data.field('teff') > 0
    print('Removing {}/{} false entries'.format(sum(~keep), len(keep)))
    hdulist[1].data = hdulist[1].data[keep]
    hdulist.flush()
    hdulist.close()


if __name__=="__main__":

    responses = ['GALEX', 'IUE', 'STROMGREN', 'JOHNSON', 'GAIA3E', 'GAIA2', 'SKYMAPPER', 'APASS', 'SDSS', '2MASS', 'WISE']
    evbs = np.r_[0:1.02:0.05]

    # Blackbody
    # calc_integrated_grid(threads=6, ebvs=evbs, law='fitzpatrick2004', Rv=3.1,
    #                      responses=responses, grid='blackbody')

    # # TMAP
    # calc_integrated_grid(threads=6, ebvs=evbs, law='fitzpatrick2004', Rv=3.1,
    #                      responses=responses, grid='tmap')

    # # Kurucz
    # calc_integrated_grid(threads=6, ebvs=evbs, law='fitzpatrick2004', Rv=3.1,
    #                      responses=responses, grid='kurucz')

    # # Munari
    # calc_integrated_grid(threads=6, ebvs=evbs, law='fitzpatrick2004', Rv=3.1,
    #                      responses=responses, grid='munari')

    # Koester
    calc_integrated_grid(threads=6, ebvs=evbs, law='fitzpatrick2004', Rv=3.1,
                         responses=responses, grid='koester')