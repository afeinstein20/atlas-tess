import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import sys
from tqdm import tqdm
from astropy.time import Time

which = '2-3'
path = f'/home/adina/.eleanor/ffis/s0092/{which}/raw'

files = np.sort(os.listdir(path))

avg = np.zeros((len(files),2))

for i in tqdm(range(len(files))):

    hdu = fits.open(os.path.join(path, files[i]))
    
    obs_date = Time(hdu[1].header['DATE-OBS'], scale='utc').mjd
    obs_end  = Time(hdu[1].header['DATE-END'], scale='utc').mjd

    #define the observing time for the middle point of the observations
    obs_obj = Time((obs_date+obs_end)/2,format='mjd',scale='utc')
    obstime_jd = obs_obj.jd

    avg[i][0] = obstime_jd
    avg[i][1] = np.nanmedian(hdu[1].data)
    hdu.close()

np.save(f'avg_{which}.npy', avg)
