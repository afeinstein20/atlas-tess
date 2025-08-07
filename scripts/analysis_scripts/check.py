import os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.time import Time
from astropy.table import Table
from astropy.coordinates import SkyCoord

import warnings
warnings.filterwarnings("ignore")

path = '/home/adina/.eleanor/ffis/s0092/2-3'#/clean_64'
files = np.sort([os.path.join(path, i) for i in os.listdir(path) if i.endswith('.fits')])

locs = Table.read('iso_tess_locations_2-3.csv', format='csv')

for i in range(len(files)):

    try:
        hdu = fits.open(files[i])
    
        tmid = (Time(hdu[1].header['DATE-OBS']).mjd + Time(hdu[1].header['DATE-END']).mjd) / 2.0

        arg = np.where(locs['mjd'] >= tmid)[0][0]

        if arg > 0:
        
            coord = SkyCoord(locs['ra'][arg], locs['dec'][arg], unit=(u.deg, u.deg))
            w = WCS(hdu[1].header)
            x, y = w.world_to_pixel(coord)
            print(i, x, y)

            if (x > 0) & (x < 2048) & (y > 0) & (y < 2048):
                hdu.close()
            else:
                hdu.close()
#                os.system('rm {}'.format(files[i]))
        else:
            hdu.close()
#            os.system('rm {}'.format(files[i]))

    except:
        pass

#    except ValueError:
#        os.system('rm {}'.format(files[i]))
#    except OSError:
#        os.system('rm {}'.format(files[i]))
#    except IndexError:
#        os.system('rm {}'.format(files[i]))
