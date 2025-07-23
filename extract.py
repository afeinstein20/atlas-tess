from astropy.io import fits
import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation
import os
from tqdm import tqdm
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS
import astropy
from astropy.stats import sigma_clip
from tqdm import tqdm

import warnings
warnings.filterwarnings("ignore")

from astropy.coordinates import SkyCoord

##############
# Find files #
##############

which = '2-3'
ffi_dir = f'/home/adina/.eleanor/ffis/s0092/{which}'
ffi_fns = np.sort([os.path.join(ffi_dir, i) for i in os.listdir(ffi_dir)])[1:1500]

ffi_dir = f'/home/adina/.eleanor/ffis/s0092/{which}'
ffi_fns_raw = np.sort([os.path.join(ffi_dir, i) for i in os.listdir(ffi_dir)])[1:1500]

locs = Table.read(f'iso_tess_locations_{which}.csv', format='csv')

rc = Table.read('rcParams.txt', format='csv')
for name, val in zip(rc['name'], rc['value']):
    plt.rcParams[name] = val

#median = np.load(f'/home/adina/.eleanor/ffis/s0092/median_{which}_1000_4000.npy')
median = np.zeros((2078,2136))

######################
# Get Pixel Location #
######################

print('Loading the data . . .')

xsize, ysize = 50,50

postcard_sub = np.zeros((len(ffi_fns)-2, xsize, ysize))
postcard = np.zeros((len(ffi_fns)-2, xsize, ysize))
not_iso = np.zeros((len(ffi_fns)-2, xsize, ysize))
not_iso_sub = np.zeros((len(ffi_fns)-2, xsize, ysize))
new_locs = np.zeros((len(ffi_fns)-2, 6))

for i in tqdm(range(len(ffi_fns)-2)):
    try:
        hdr = fits.getheader(ffi_fns[i], ext=1)
    
        time = ( Time(hdr['DATE-OBS']).mjd + Time(hdr['DATE-END']).mjd) / 2.0
        arg = np.where(locs['mjd'] >= time)[0][0]
        coord = SkyCoord(locs['ra'][arg], locs['dec'][arg], unit=(u.deg, u.deg))
        coord = SkyCoord(295.016, -19.075, unit=(u.deg, u.deg))
        x, y = WCS(hdr).world_to_pixel(coord)
        #x, y = 1848, 1457
        #print(x, y)
        
        new_locs[i] = [x, y, coord.ra.deg, coord.dec.deg, locs['apmag'][arg], time]

        dat = fits.getdata(ffi_fns[i])
        
        postcard[i] = dat[int(y-xsize/2):int(y+xsize/2), int(x-ysize/2):int(x+ysize/2)]
        postcard_sub[i] = (dat-median)[int(y-xsize/2):int(y+xsize/2), int(x-ysize/2):int(x+ysize/2)]

        x += 100
        y += 100

        not_iso[i] = dat[int(y-xsize/2):int(y+xsize/2), int(x-ysize/2):int(x+ysize/2)]
        not_iso_sub[i] = (dat-median)[int(y-xsize/2):int(y+xsize/2), int(x-ysize/2):int(x+ysize/2)]
        del dat
        
    except TypeError:
        print(i)

    del hdr

#########################
# Plot the light curves #
#########################

print('Making the plot . . .')

q = new_locs[:,-1] > 60800.1

sc = sigma_clip(postcard, sigma=1.0)
bkg = np.nanmean(postcard*~sc.mask, axis=(1,2))

fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(14,7), sharex=True)
fig.set_facecolor('w')

ax1.plot(new_locs[:,-1][q], np.nansum(postcard[:,5:7,5:7], axis=(1,2))[q], color='k',
         label='3I/ATLAS', lw=2, zorder=10)
ax1.plot(new_locs[:,-1][q], np.nansum(not_iso[:,5:7,5:7], axis=(1,2))[q], color='deepskyblue',
         label='not 3I/ATLAS')

ax2.plot(new_locs[:,-1][q], bkg[q], color='darkorange', label='Average Background')

ax1.set_ylabel('Counts')
ax2.set_ylabel('Counts')

ax1.legend()
ax2.legend()

#plt.xlim(60816.1,)

ax1.set_ylim(1000,10000)
ax2.set_ylim(50,)

ax2.set_xlabel('Time [MJD]')
plt.savefig('3I-ATLAS_TESS_comparison_{}_v3.png'.format(which), bbox_inches='tight', dpi=300)

##############
# Make a GIF #
##############

print('making the GIF . . .')

fig = plt.figure(figsize=(10,10))
frames = []

for i in range(len(postcard)):
    frames.append([plt.imshow(postcard[i], origin='lower',
                              aspect='auto', 
                              vmin=-100, vmax=1000)])

print('Writing animation')
ani = animation.ArtistAnimation(fig, frames, interval=50, blit=True,
                                repeat_delay=1000)
ani.save("frame-ffi-{}_v3.gif".format(which))
