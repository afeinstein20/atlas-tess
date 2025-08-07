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

which = '2-3'
curl = f'{which}/tesscurl_sector_92_{which}.sh'

lines = np.array(open(curl, 'r').readlines())
args = np.arange(0, len(lines), 2, dtype=int)
lines = lines[args]

tab = Table(names=['file', 'x', 'y'], dtype=[str, float, float])

inframe = True

while inframe is True:

    for i in tqdm(range(len(lines))):
   
        fn = lines[i].split(' ')[5]

        if os.path.isfile(fn) == False:
            os.system(lines[i])

        hdu = fits.open(fn)
        wcs = WCS(hdu[1].header)

        obs_date = Time(hdu[1].header['DATE-OBS'], scale='utc').mjd
        obs_end  = Time(hdu[1].header['DATE-END'], scale='utc').mjd

        #define the observing time for the middle point of the observations
        obs_obj = Time((obs_date+obs_end)/2,format='mjd',scale='utc')
        obstime_jd = obs_obj.jd

        #Query Horizons for the orbit
        obj = Horizons(id='3I',id_type='smallbody',location='@tess',epochs=obstime_jd)
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
            tab.add_row([fn, pix_cen[0], pix_cen[1]])
            tab.write(f'locs_{which}.csv', format='csv', overwrite=True)

