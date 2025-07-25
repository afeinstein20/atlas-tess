import numpy as np
from astropy.table import Table, Column
import matplotlib.pyplot as plt
from astropy.io import fits
import os, sys
from tqdm import tqdm
import pandas as pd
from scipy.signal import savgol_filter
import matplotlib.animation as animation
from astropy.time import Time

which='1-2'
#which = '2-3'
target='3I'
flipped = True
inpath = f'/home/adina/.eleanor/ffis/s0092/{which}/raw'
outpath = os.path.join(inpath[:-4], 'pixel')

if which == '2-3':
    start, end = 0, 4399
elif which == '1-2':
    start, end = 280, 5718

files = np.sort([os.path.join(inpath, i) for i in os.listdir(inpath) if i.endswith('.fits')])
files = files[start:end]

coords = Table.read(f'/home/adina/.eleanor/ffis/s0092/locs_{which}.csv', format='csv')
coords.sort('file')
_, u = np.unique(coords['file'], return_index=True)
coords = coords[u]

search_fn_start = '/'.join(e for e in files[0].split('/')[-3:])
search_fn_end   = '/'.join(e for e in files[-1].split('/')[-3:])
start_arg = np.where(coords['file'] == search_fn_start)[0][0]
end_arg   = np.where(coords['file'] == search_fn_end)[0][0]
coords = coords[start_arg : end_arg]

pad = 20
xlow, xupp = np.nanmin(coords['x']) - pad, np.nanmax(coords['x']) + pad
ylow, yupp = np.nanmin(coords['y']) - pad, np.nanmax(coords['y']) + pad

print(start_arg, end_arg)
print(xlow, xupp, ylow, yupp)

if xlow < 0:
    xlow = 0
if ylow < 0:
    ylow = 0

if flipped:
    ext = '_flipped.npy'
else:
    ext = '.npy'

output_name = f'pixels_{which}_{target}_xlow-{int(xlow)}_ylow-{int(ylow)}{ext}'

output = {}

times = np.zeros(len(files))

######################
# Load in the pixels #
######################
for i in tqdm(range(len(files))):
    hdu = fits.open(files[i])

    if flipped:
        cutout = hdu[1].data[int(ylow):int(yupp), int(xlow):int(xupp)]
        output['shape'] = (len(files), cutout.shape[0], cutout.shape[1])
        output['extent'] = (ylow, yupp, xlow, xupp)
    else:
        cutout = hdu[1].data[int(xlow):int(xupp), int(ylow):int(yupp)]
        output['shape'] = (len(files), cutout.shape[0], cutout.shape[1])
        output['extent'] = (xlow, xupp, ylow, yupp)

    if i == 0:
        stacked_cutouts = np.zeros((len(files), len(cutout.flatten())))

    #stacked_cutouts[i] = cutout * 1.0
    stacked_cutouts[i] = cutout.flatten() * 1.0

    date_start = Time(hdu[1].header['DATE-OBS'], scale='utc').mjd
    date_end   = Time(hdu[1].header['DATE-END'], scale='utc').mjd
    times[i] = (date_start + date_end) / 2.0

    hdu.close()

output['pixels'] = stacked_cutouts
output['time'] = times

###########################################
# Create a table of models for each pixel #
###########################################
models = np.zeros(stacked_cutouts.shape)

for i in tqdm(range(stacked_cutouts.shape[1])):
    m = savgol_filter(stacked_cutouts[:,i], window_length=307, polyorder=2)
    models[:,i] = m

output['models'] = models
np.save(output_name, output)

##################################
# Remove the savgol filter trend #
##################################
"""
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
"""
#################################################
# Create the FITS files with the removed trends #
#################################################
"""
removed = stacked_cutouts.reshape(output['shape']) - models.reshape(output['shape'])

for i in tqdm(range(len(files))):
    hdu = fits.open(files[i])

    if flipped:
        ext = '_flipped_rmv.fits'
    else:
        ext = '_rmv.fits'

    new_fits = files[i].split('/')[-1].split('.fits')[0] + ext
    new_fits = os.path.join(outpath, new_fits)

    #fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(14,4), sharex=True, sharey=True)
    #ax1.imshow(hdu[1].data[y1-pad:y2+pad, x1-pad:x2+pad], vmin=0, vmax=1000)
    #ax2.imshow(stacked_cutouts[i].reshape((yshape, xshape)), vmin=0, vmax=1000)
    #ax3.imshow(reshaped, vmin=0, vmax=100)
    #plt.show()
    if flipped:
        hdu[1].data[int(ylow):int(yupp), int(xlow):int(xupp)] = removed[i]
    else:
        hdu[1].data[int(xlow):int(xupp), int(ylow):int(yupp)] = removed[i]

    hdu.writeto(new_fits, overwrite=True)
    hdu.close()
"""
