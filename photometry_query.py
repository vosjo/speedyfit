
import ConfigParser

import numpy as np

from astroquery.vizier import Vizier
from astroquery.utils.tap.core import TapPlus

from astropy import units as u

# always request the distance of the observations to the search coordinates and sort on increasing distance
v = Vizier(columns=["*", '+_r']) 

#-- read in catalog information
cat_info = ConfigParser.ConfigParser()
cat_info.optionxform = str # make sure the options are case sensitive
cat_info.readfp(open('vizier_cats_phot.cfg'))


tap_info = ConfigParser.ConfigParser()
tap_info.optionxform = str # make sure the options are case sensitive
tap_info.readfp(open('tap_cats_phot.cfg'))


def get_vizier_photometry(objectname, radius=5):
   
   catalogs = cat_info.sections()

   data = v.query_object(objectname, catalog=catalogs, radius=radius*u.arcsec)
   
   photometry = []
   
   for catalog in data.keys():
      
      distance = (data[catalog]['_r'][0] * u.arcmin).to(u.arcsec).value
      
      if cat_info.has_option(catalog, 'bibcode'):
         bibcode = cat_info.get(catalog, 'bibcode')
      else:
         bibcode = ''
      
      for band in cat_info.options(catalog):
         
         if '_unit' in band: continue
         if '_err' in band: continue
         if 'bibcode' in band: continue
         
         bandname = cat_info.get(catalog, band)
         value = data[catalog][band][0]
         
         errkw = cat_info.get(catalog, band+'_err') if cat_info.has_option(catalog, band+'_err') else 'e_'+band
         try:
            err = data[catalog][errkw][0]
         except:
            err = 0.05
         
         if cat_info.has_option(catalog, band+'_unit'):
            unit = cat_info.get(catalog, band+'_unit')
         else:
            unit = data[catalog][band].unit
         
         photometry.append(( bandname, value, err, unit, distance, bibcode ))
   
   dtypes = [('band', 'a20'), ('meas', 'f8'), ('emeas', 'f8'), ('unit', 'a10'), ('distance', 'f8'), ('bibcode', 'a20')]
   photometry = np.array(photometry, dtype=dtypes)
   
   return photometry


def tap_query(ra, dec, catalog):
   
   tp = TapPlus(url=catalog)
   
   table = tap_info.get(catalog, 'table')
   
   keywords = ""
   bands = []
   for band in tap_info.options(catalog):
      if 'table' in band: continue
      if '_unit' in band: continue
      if '_err' in band: continue
      if 'bibcode' in band: continue
      keywords += band + ", " + "e_" + band + ", "
      bands.append(band)
      
   if keywords[-2:] == ", ":
      keywords = keywords[0:-2]
   
   query = """SELECT 
   DISTANCE(POINT('ICRS', raj2000, dej2000),
            POINT('ICRS', {ra:}, {dec:})) AS dist, {kws:}
   FROM {table:} AS m
   WHERE 
      1=CONTAINS(POINT('ICRS', raj2000, dej2000),
                 CIRCLE('ICRS', {ra:}, {dec:}, 10/3600. ))
   ORDER BY dist""".format(ra=ra, dec=dec, table=table, kws=keywords)
   
   job = tp.launch_job(query) 
   
   results = job.get_results()
   
   distance = (results['dist'][0] * u.degree).to(u.arcsec).value
   bibcode = tap_info.get(catalog, 'bibcode')
   
   photometry = []
   for band in bands:
      bandname = tap_info.get(catalog, band)
      value = results[band][0]
      err = results["e_"+band][0]
      unit = results[band].unit
      
      
      photometry.append(( bandname, value, err, unit, distance, bibcode))
      
   dtypes = [('band', 'a20'), ('meas', 'f8'), ('emeas', 'f8'), ('unit', 'a10'), ('distance', 'f8'), ('bibcode', 'a20')]
   photometry = np.array(photometry, dtype=dtypes)

   return photometry

def get_tap_photometry(ra, dec):
   catalogs = tap_info.sections()
   
   dtypes = [('band', 'a20'), ('meas', 'f8'), ('emeas', 'f8'), ('unit', 'a10'), ('distance', 'f8'), ('bibcode', 'a20')]
   photometry = np.array([], dtype=dtypes)
   
   for catalog in catalogs:
      phot_ = tap_query(ra, dec, catalog)
      photometry = np.hstack([photometry, phot_])
      
   return photometry
      

photometry = get_tap_photometry(030.3931577216772, -53.7287603837835)

photometry_ = get_vizier_photometry("JL277")

photometry = np.hstack([photometry, photometry_])

print photometry
