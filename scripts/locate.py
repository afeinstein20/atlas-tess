import os
import sys
import numpy as np
from astropy.io import fits
from astroquery.jplhorizons import Horizons
from astropy.wcs import WCS
from astropy.time import Time
import astropy.units as u
from tqdm import tqdm
from astropy.coordinates import SkyCoord
from astropy.table import Table
import warnings

warnings.filterwarnings('ignore')

which = '2-3'
path = f'{which}/raw'
files = np.sort([os.path.join(path, i) for i in os.listdir(path) if i.endswith('.fits')])
tab = Table(names=['file', 'x', 'y'], dtype=[str, float, float])

inframe = True

while inframe is True:

    for i in tqdm(range(len(files))):

        #if os.path.isfile(fn) == False:
        #    os.system(lines[i])

        hdu = fits.open(files[i])
        wcs = WCS(hdu[1].header)

        obs_date = Time(hdu[1].header['DATE-OBS'], scale='utc').mjd
        obs_end  = Time(hdu[1].header['DATE-END'], scale='utc').mjd

        #define the observing time for the middle point of the observations
        obs_obj = Time((obs_date+obs_end)/2,format='mjd',scale='utc')
        obstime_jd = obs_obj.jd

        #Query Horizons for the orbit
        obj = Horizons(id='A918 PE',id_type='smallbody',location='@tess',epochs=obstime_jd)
        eph = obj.ephemerides()

        #Create a SkyCoord for the object at the midpoint of the observations
        c = SkyCoord(eph['RA'][0], eph['DEC'][0], frame='icrs', unit='deg')

        cen_skycoord = SkyCoord(c.ra,c.dec,frame='icrs', unit=(u.deg, u.deg))

        try:
            # convert the objects SkyCoord to the frame pixels
            pix_cen = wcs.world_to_pixel(cen_skycoord)
        except:
            hdu.close()
            os.system(f'rm {fn}')

        hdu.close()

        if ( (pix_cen[0] < 0) | (pix_cen[0] > 2130) |
             (pix_cen[1] < 0) | (pix_cen[1] > 2130) ):
            os.system(f'rm {fn}')
            inframe = False

        else:
            tab.add_row([files[i], pix_cen[0], pix_cen[1]])
            tab.write(f'locs_{which}_asteroid.csv', format='csv', overwrite=True)

    inframe=False
