import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import sys
from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astroquery.jplhorizons import Horizons
import warnings

warnings.filterwarnings("ignore")

which = "1-2"
path = f'/home/adina/.eleanor/ffis/s0092/{which}/raw'
files = np.sort([os.path.join(path, i) for i in os.listdir(path) if i.endswith('.fits')])

tab = Table(names=['file', 'mjd', 'jd', 'ra', 'dec', 'x', 'y'],
            dtype=[str, float, float, float, float, float, float])

for i in tqdm(range(len(files))):

    hdu = fits.open(files[i])
    w = WCS(hdu[1].header)

    tstart = Time(hdu[1].header['DATE-OBS'], scale='utc').mjd
    tend   = Time(hdu[1].header['DATE-END'], scale='utc').mjd
    
    obs_time = Time( (tstart+tend)/2.0, format='mjd', scale='utc')
    obs_time_jd = obs_time.jd

    obj = Horizons(id='3I', id_type='smallbody', location='@tess', epochs=obs_time_jd)
    eph = obj.ephemerides()

    sc = SkyCoord(eph['RA'][0], eph['DEC'][0], unit=(u.deg,u.deg), frame='icrs')
    x, y = w.world_to_pixel(sc)

    row = [files[i], obs_time.mjd, obs_time_jd,
           sc.ra.deg, sc.dec.deg, x, y]
    tab.add_row(row)

    hdu.close()

tab.write(f'3I_pixel_locations_{which}.csv', format='csv', overwrite=True)
