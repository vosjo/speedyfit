
import os
import shutil
import logging
import configparser
import warnings

import numpy as np

from . import filters

from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astroquery.utils.tap.core import TapPlus
import pyvo

from astropy import units as u
from astropy.coordinates.angles import Angle
from astropy.coordinates import SkyCoord
from astropy.io import ascii

from numpy.lib.recfunctions import append_fields

from zero_point import zpt #gaiadr3-zeropoint
zpt.load_tables()

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
VIZIER_CAT = 'vizier_cats_phot.cfg'
TAP_CAT = 'tap_cats_phot.cfg'

# always request the distance of the observations to the search coordinates and sort on increasing distance
v = Vizier(columns=["*", '+_r'])

#-- read in catalog information
def load_photometry_catalogs():
    """
    Function to read the photometry catalogs: vizier_cats_phot.cfg and tap_cats_phot.cfg.
    The function will first look in the directory set with the SPEEDYFIT_MODELS environment variable. If no catalog
    files can be found in that directory, it will use the catalogs shipped with the speedyfit package.
    """
    user_dir = os.environ.get('SPEEDYFIT_MODELS', None)

    if user_dir is not None and os.path.exists(os.path.join(user_dir, VIZIER_CAT)):
        vizier_file = os.path.join(user_dir, VIZIER_CAT)
    else:
        vizier_file = os.path.join(FILE_DIR, VIZIER_CAT)

    if user_dir is not None and os.path.exists(os.path.join(user_dir, TAP_CAT)):
        tap_file = os.path.join(user_dir, TAP_CAT)
    else:
        tap_file = os.path.join(FILE_DIR, TAP_CAT)

    print(f"Going to read Vizier photometry catalog from {vizier_file}.")
    viz_info = configparser.ConfigParser()
    viz_info.optionxform = str # make sure the options are case sensitive
    with open(vizier_file) as ifile:
        viz_info.readfp(ifile)

    print(f"Going to read TAP photometry catalog from {tap_file}.")
    tap_info = configparser.ConfigParser()
    tap_info.optionxform = str  # make sure the options are case sensitive
    with open(tap_file) as ifile:
        tap_info.readfp(ifile)

    return viz_info, tap_info
VIZ_INFO, TAP_INFO = load_photometry_catalogs()

def copy_photometry_catalogs(overwrite=False):
    """
    Copy the two photometry catalogs: vizier_cats_phot.cfg and tap_cats_phot.cfg from the package source location
    to the SPEEDYFIT_MODELS directory, so that the user can edit them.
    """
    if os.environ.get('SPEEDYFIT_MODELS', None) is None:
        print("SPEEDYFIT_MODELS env variable is not set, can not copy the photometry catalogs.")
        return
    else:
        destination = os.environ.get('SPEEDYFIT_MODELS', None)

    vis_cat = os.path.join(FILE_DIR, VIZIER_CAT)
    if os.path.exists(vis_cat) and not overwrite:
        print(f"The file {vis_cat} already exists. To overwrite existing files use option --overwrite")
    elif (os.path.exists(vis_cat) and overwrite) or not os.path.exists(vis_cat):
        shutil.copy(vis_cat, destination)
        print(f"Copied {vis_cat} to {destination}")

    tap_cat = os.path.join(FILE_DIR, TAP_CAT)
    if os.path.exists(tap_cat) and not overwrite:
        print(f"The file {tap_cat} already exists. To overwrite existing files use option --overwrite")
    elif (os.path.exists(tap_cat) and overwrite) or not os.path.exists(tap_cat):
        shutil.copy(tap_cat, destination)
        print(f"Copied {tap_cat} to {destination}")


def get_coordinate(objectname):

    if objectname[0] == 'J' and objectname[1] != 'L':
        name = objectname[1:]
        if '+' in name:
            ra, dec = name.split('+')
            sign = '+'
        else:
            ra, dec = name.split('-')
            sign = '-'

        ra = "{}:{}:{}".format(ra[0:2], ra[2:4], ra[4:])
        dec = "{}{}:{}:{}".format(sign, dec[0:2], dec[2:4], dec[4:])

        ra = Angle(ra, unit='hour').degree
        dec = Angle(dec, unit='degree').degree

    else:
        data = Simbad.query_object(objectname)
        ra = Angle(data['RA'][0], unit='hour').degree
        dec = Angle(data['DEC'][0], unit='degree').degree

    return ra, dec


def get_vizier_photometry(objectname, radius=5):

    catalogs = VIZ_INFO.sections()

    data = v.query_object(objectname, catalog=catalogs, radius=radius*u.arcsec)

    photometry = []

    for catalog in list(data.keys()):

        distance = (data[catalog]['_r'][0] * u.arcmin).to(u.arcsec).value

        if VIZ_INFO.has_option(catalog, 'bibcode'):
            bibcode = VIZ_INFO.get(catalog, 'bibcode')
        else:
            bibcode = '-'

        print(catalog)

        for band in VIZ_INFO.options(catalog):

            if '_unit' in band: continue
            if '_err' in band: continue
            if 'bibcode' in band: continue
            if '-' in band: continue

            bandname = VIZ_INFO.get(catalog, band)
            value = data[catalog][band][0]

            errkw = VIZ_INFO.get(catalog, band + '_err') if VIZ_INFO.has_option(catalog, band + '_err') else 'e_' + band
            try:
                err = data[catalog][errkw][0]
            except:
                err = 0.05

            if VIZ_INFO.has_option(catalog, band + '_unit'):
                unit = VIZ_INFO.get(catalog, band + '_unit')
            else:
                unit = data[catalog][band].unit

            photometry.append(( bandname, value, err, unit, distance, bibcode ))

    dtypes = [('band', '<U20'), ('meas', 'f8'), ('emeas', 'f8'), ('unit', '<U10'), ('distance', 'f8'), ('bibcode', '<U20')]
    photometry = np.array(photometry, dtype=dtypes)

    return photometry


def tap_query_vo(ra, dec, catalog):
    """
    info see:
    http://docs.g-vo.org/pyvo/html/page011.html
    """

    service = pyvo.dal.TAPService(catalog)

    table = TAP_INFO.get(catalog, 'table')
    print(table)

    keywords = ""
    bands = []
    for band in TAP_INFO.options(catalog):
        if 'table' in band: continue
        if 'rakw' in band: continue
        if 'deckw' in band: continue
        if '_unit' in band: continue
        if '_err' in band: continue
        if 'bibcode' in band: continue

        errkw = TAP_INFO.get(catalog, band + '_err') if TAP_INFO.has_option(catalog, band + '_err') else 'e_' + band

        keywords += band + ", " + errkw + ", "
        bands.append(band)

    if keywords[-2:] == ", ":
        keywords = keywords[0:-2]

    rakw = TAP_INFO.get(catalog, 'rakw') if TAP_INFO.has_option(catalog, 'rakw') else 'raj2000'
    deckw = TAP_INFO.get(catalog, 'deckw') if TAP_INFO.has_option(catalog, 'deckw') else 'dej2000'

    # -- first try to query with distance
    query = """SELECT 
         DISTANCE(POINT('ICRS', {rakw:}, {deckw}),
                  POINT('ICRS', {ra:}, {dec:})) AS dist, {kws:}
         FROM {table:} AS m
         WHERE 
            1=CONTAINS(POINT('ICRS', {rakw:}, {deckw}),
                     CIRCLE('ICRS', {ra:}, {dec:}, 0.005 ))
         ORDER BY dist""".format(ra=ra, dec=dec, table=table, rakw=rakw, deckw=deckw, kws=keywords)

    results = service.run_sync(query).to_table()

    if len(results) == 0:
        dtypes = [('band', '<U20'), ('meas', 'f8'), ('emeas', 'f8'), ('unit', '<U10'), ('distance', 'f8'), ('bibcode', '<U20')]
        return np.array([], dtype=dtypes)

    # -- get the distance from the query result or calculate it.
    if 'dist' in results.dtype.names:
        distance = (results['dist'][0] * u.degree).to(u.arcsec).value
    else:
        c1 = SkyCoord(ra * u.degree, dec * u.degree)
        distance = SkyCoord(results['ra'][0] * u.degree, results['de'][0] * u.degree).separation(c1).to(u.arcsec).value

    bibcode = TAP_INFO.get(catalog, 'bibcode')

    # -- get the photometric measurements, errors and units
    photometry = []
    for band in bands:
        bandname = TAP_INFO.get(catalog, band)
        value = results[band][0]

        errkw = TAP_INFO.get(catalog, band + '_err') if TAP_INFO.has_option(catalog, band + '_err') else 'e_' + band
        err = results[errkw][0]

        if TAP_INFO.has_option(catalog, band + '_unit'):
            unit = TAP_INFO.get(catalog, band + '_unit')
        else:
            unit = results[band].unit

        photometry.append((bandname, value, err, unit, distance, bibcode))

    dtypes = [('band', '<U20'), ('meas', 'f8'), ('emeas', 'f8'), ('unit', '<U10'), ('distance', 'f8'), ('bibcode', '<U20')]
    photometry = np.array(photometry, dtype=dtypes)

    return photometry


def tap_query(ra, dec, catalog):

    tp = TapPlus(url=catalog)

    table = TAP_INFO.get(catalog, 'table')
    print(table)

    keywords = ""
    bands = []
    for band in TAP_INFO.options(catalog):
        if 'table' in band: continue
        if 'rakw' in band: continue
        if 'deckw' in band: continue
        if '_unit' in band: continue
        if '_err' in band: continue
        if 'bibcode' in band: continue

        errkw = TAP_INFO.get(catalog, band + '_err') if TAP_INFO.has_option(catalog, band + '_err') else 'e_' + band

        keywords += band + ", " + errkw + ", "
        bands.append(band)

    if keywords[-2:] == ", ":
        keywords = keywords[0:-2]

    rakw = TAP_INFO.get(catalog, 'rakw') if TAP_INFO.has_option(catalog, 'rakw') else 'raj2000'
    deckw = TAP_INFO.get(catalog, 'deckw') if TAP_INFO.has_option(catalog, 'deckw') else 'dej2000'

    try:
        #-- first try to query with distance
        query = """SELECT 
      DISTANCE(POINT('ICRS', {rakw:}, {deckw}),
               POINT('ICRS', {ra:}, {dec:})) AS dist, {kws:}
      FROM {table:} AS m
      WHERE 
         1=CONTAINS(POINT('ICRS', {rakw:}, {deckw}),
                  CIRCLE('ICRS', {ra:}, {dec:}, 0.005 ))
      ORDER BY dist""".format(ra=ra, dec=dec, table=table, rakw=rakw, deckw=deckw, kws=keywords)

        job = tp.launch_job(query)

        results = job.get_results()

    except Exception as e:

        print("problem", e)

        #-- not all catalogs accept a distance sorted query, we have to hope that the correct star is returned.
        query = """SELECT {rakw:} as RA, {deckw} as DE, {kws:}
      FROM {table:} AS m
      WHERE 
         1=CONTAINS(POINT('ICRS', {rakw:}, {deckw}),
                  CIRCLE('ICRS', {ra:}, {dec:}, 0.005 ))
      """.format(ra=ra, dec=dec, table=table, rakw=rakw, deckw=deckw, kws=keywords)

        job = tp.launch_job(query)

        results = job.get_results()

    if len(results) == 0:
        dtypes = [('band', 'a20'), ('meas', 'f8'), ('emeas', 'f8'), ('unit', 'a10'), ('distance', 'f8'), ('bibcode', 'a20')]
        return np.array([], dtype=dtypes)

    #-- get the distance from the query result or calculate it.
    if 'dist' in results.dtype.names:
        distance = (results['dist'][0] * u.degree).to(u.arcsec).value
    else:
        c1 = SkyCoord(ra * u.degree, dec * u.degree)
        distance = SkyCoord(results['ra'][0] * u.degree, results['de'][0] * u.degree).separation(c1).to(u.arcsec).value


    bibcode = TAP_INFO.get(catalog, 'bibcode')

    #-- get the photometric measurements, errors and units
    photometry = []
    for band in bands:
        bandname = TAP_INFO.get(catalog, band)
        value = results[band][0]

        errkw = TAP_INFO.get(catalog, band + '_err') if TAP_INFO.has_option(catalog, band + '_err') else 'e_' + band
        err = results[errkw][0]

        if TAP_INFO.has_option(catalog, band + '_unit'):
            unit = TAP_INFO.get(catalog, band + '_unit')
        else:
            unit = results[band].unit

        photometry.append(( bandname, value, err, unit, distance, bibcode))

    dtypes = [('band', '<U20'), ('meas', 'f8'), ('emeas', 'f8'), ('unit', '<U10'), ('distance', 'f8'), ('bibcode', '<U20')]
    photometry = np.array(photometry, dtype=dtypes)

    return photometry


def get_tap_photometry(ra, dec):
    catalogs = TAP_INFO.sections()

    dtypes = [('band', '<U20'), ('meas', 'f8'), ('emeas', 'f8'), ('unit', '<U10'), ('distance', 'f8'), ('bibcode', '<U20')]
    photometry = np.array([], dtype=dtypes)

    for catalog in catalogs:

        phot_ = tap_query_vo(ra, dec, catalog)
        photometry = np.hstack([photometry, phot_])

    return photometry


def get_photometry(objectname, filename=None):

    #-- get the object coodinates
    ra, dec = get_coordinate(objectname)

    #-- query tap catalogs
    photometry = get_tap_photometry(ra, dec)

    #-- query Vizier catalogs
    photometry_ = get_vizier_photometry("{} {}".format(ra, dec))

    print(photometry)
    print(photometry_)

    photometry = np.hstack([photometry, photometry_])

    #-- convert magnitudes to fluxes
    wave, flux, err = [], [], []
    for band, meas, emeas, unit in zip(photometry['band'], photometry['meas'], photometry['emeas'], photometry['unit']):
        if np.isnan(emeas) or emeas < 0: emeas = 0.02

        f_, e_ = filters.mag2flux(meas, emeas, band)

        wave.append(filters.eff_wave(band))
        flux.append(f_)
        err.append(e_)

    photometry = append_fields(photometry, ['flux', 'eflux'], data=[flux, err],
                               dtypes=['f8', 'f8'], usemask=False)

    if filename is not None:
        ascii.write(photometry, filename, format='fixed_width', overwrite=True)

    return photometry


def get_parallax(objectname, radius=5):

    v_gaia = Vizier(columns=["Plx", "e_Plx", '+_r', 'Gmag', 'nueff', 'pscol', 'ELAT', 'Solved'])

    data = v_gaia.query_object(objectname, catalog=['I/350/gaiaedr3'], radius=radius*u.arcsec)

    if len(data) == 0:
        return None, None

    data = data['I/350/gaiaedr3'][0]

    plx, plx_e = data['Plx'], data['e_Plx']

    if not data['pscol']:
        data['pscol'] = 0

    # apply parallax zero point correction of Lindgren. If not possible, use the average offset.
    # https://arxiv.org/pdf/2012.01742.pdf
    # https://gitlab.com/icc-ub/public/gaiadr3_zeropoint
    try:
        zp = zpt.get_zpt(data['Gmag'], data['nueff'], data['pscol'], data['ELAT'], data['Solved'])
    except Exception:
        warnings.warn("Could not calculate the parallax zero point offset based on Lindgren+2020, using average")
        zp = 0.02
    plx -= zp

    return plx, plx_e

