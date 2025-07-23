import numpy as np
from astropy.table import Table, Column
import matplotlib.pyplot as plt
from astropy.io import fits
import os, sys
from tqdm import tqdm
import pandas as pd
from scipy.signal import savgol_filter
import matplotlib.animation as animation


which='2-3'
target='A918PE'
inpath = f'/home/adina/.eleanor/ffis/s0092/{which}/{target}'
outpath = os.path.join(inpath, target)

files = np.sort([os.path.join(inpath, i) for i in os.listdir(inpath) if i.endswith('.fits')])
files = files

tab = Table.read(f'JPL_Horizons_{which}_{target}.csv', format='csv')
y1, y2 = int(np.nanmin(tab['y'])), int(np.nanmax(tab['y']))
x1, x2 = int(np.nanmin(tab['x'])), int(np.nanmax(tab['x']))

xpad = 48 # extra pixel padding
ypad = 30

xshape = (x2+xpad) - (x1-xpad)
yshape = (y2+ypad) - (y1-ypad)

pixel_filename = f'pixels_{which}_{target}_x1-{x1}_x2-{x2}_y1-{y1}_y2-{y2}.npy'
model_filename = f'models_{which}_{target}_x1-{x1}_x2-{x2}_y1-{y1}_y2-{y2}.npy'

######################
# Load in the pixels #
######################

if os.path.isfile(pixel_filename) == True:
    stacked_cutouts = np.load(pixel_filename) # Loads in file if it exists

else:
    for i in tqdm(range(len(files))):
        hdu = fits.open(files[i])

        cutout = hdu[1].data[x1-xpad:x2+xpad, y1-ypad:y2+ypad].flatten()

        if i == 0:
            stacked_cutouts = np.zeros((len(files), len(cutout)))

        stacked_cutouts[i] = cutout * 1.0

        hdu.close()

    np.save(pixel_filename, stacked_cutouts)

###########################################
# Create a table of models for each pixel #
###########################################

if os.path.isfile(model_filename) == True:
    models = np.load(model_filename) # Loads in file if it exists
else:
    models = np.zeros(stacked_cutouts.shape)

    for i in tqdm(range(stacked_cutouts.shape[1])):
        m = savgol_filter(stacked_cutouts[:,i], window_length=307, polyorder=2)
        models[:,i] = m

    np.save(model_filename, models)

##################################
# Remove the savgol filter trend #
##################################
removed = stacked_cutouts - models

removed = removed.reshape( (len(files), xshape, yshape))

print('Making the GIF ...')
fig = plt.figure(figsize=(12,10))
frames = []

for i in range(500,1000):#len(files)):
    frames.append([plt.imshow(removed[i], origin='lower',
                              aspect='auto', vmin=-10, vmax=10),
                   plt.title(i)])

ani = animation.ArtistAnimation(fig, frames, interval=50, blit=True,
                                repeat_delay=1000)
ani.save(f'cutout_{which}_{target}.gif')

#################################################
# Create the FITS files with the removed trends #
#################################################
"""
for i in tqdm(range(len(files))):
    hdu = fits.open(files[i])

    new_fits = files[i].split('/')[-1].split('.fits')[0] + '_rmv.fits'
    new_fits = os.path.join(outpath, new_fits)

    reshaped = removed[i].reshape((yshape, xshape))

    #fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(14,4), sharex=True, sharey=True)
    #ax1.imshow(hdu[1].data[y1-pad:y2+pad, x1-pad:x2+pad], vmin=0, vmax=1000)
    #ax2.imshow(stacked_cutouts[i].reshape((yshape, xshape)), vmin=0, vmax=1000)
    #ax3.imshow(reshaped, vmin=0, vmax=100)
    #plt.show()

    hdu[1].data[y1-pad:y2+pad, x1-pad:x2+pad] = reshaped
    hdu.writeto(new_fits, overwrite=True)
    hdu.close()

    os.system('rm {}'.format(files[i]))
"""
