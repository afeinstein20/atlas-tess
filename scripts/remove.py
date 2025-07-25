import os
import numpy as np
from tqdm import tqdm
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.time import Time
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astroquery.jplhorizons import Horizons

import warnings
warnings.filterwarnings("ignore")

path = '/home/adina/.eleanor/ffis/s0092/1-2/raw' #clean_64'
files = np.sort([os.path.join(path, i) for i in os.listdir(path) if i.endswith('.fits')])

bad = 0

for i in tqdm([0, len(files)-1]):#range(len(files))):

    hdu = fits.open(files[i])
    
    tmid = (Time(hdu[1].header['DATE-OBS']).jd + Time(hdu[1].header['DATE-END']).jd) / 2.0
    print(tmid)
    obj = Horizons(id='C/2025 N1',
                   id_type='smallbody',
                   location='@tess',
                   epochs=tmid)

    eph = obj.ephemerides()

    w = WCS(hdu[1].header)

    coord = SkyCoord(eph['RA'][0], eph['DEC'][0],
                     frame='icrs', unit='deg')

    try:
        x, y = w.world_to_pixel(coord)

        if (x > 0) & (x < 2048) & (y > 0) & (y < 2048):
            hdu.close()
        else:
            hdu.close()
            os.system('rm {}'.format(files[i]))
            bad += 1
    except ValueError:
        os.system('rm {}'.format(files[i]))

    #print(bad)

print(bad, len(files))

    #except ValueError:
    #    os.system('rm {}'.format(files[i]))
    #except OSError:
    #    os.system('rm {}'.format(files[i]))
    #except IndexError:
    #    os.system('rm {}'.format(files[i]))
