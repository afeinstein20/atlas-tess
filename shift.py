import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import os, sys
import matplotlib.animation as animation
from tqdm import tqdm
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from scipy.ndimage import shift

import warnings
warnings.filterwarnings("ignore")

def load_data(files):
    """
    Loads in the FITS file data and headers
    """
    images = np.zeros((len(files), 2048, 2048))
    headers= []
    tmid = np.zeros(len(files))

    for i in tqdm(range(len(files))):
        with fits.open(files[i]) as hdu:
            images[i] = hdu[0].data
            headers.append(hdu[0].header)

            st = Time(hdu[0].header['DATE-OBS']).mjd
            sp = Time(hdu[0].header['DATE-END']).mjd
            tmid[i] = (st + sp) / 2.0

        hdu.close()

    return images, headers, tmid

def center_object(files, table_file, size=51, method='sum',
                  start=0, end=1381):
    """
    A routine to shift the object into the center of the frame.

    files : np.array, list
       Array or list of FITS file names, including the path.
    table_file : str
       Name of the file where the (RA, Dec) positions are stored 
       for the object as a function of time.
    size : int, optional
       The size of the cutout. Default is 51 pixels.
    method : str, optional
       The routine for stacking the image. Default is 'sum'.
    start : int, optional
       The first FITS file to include in the analysis. Default 
       is 0.
    end : int, optional 
       The last FITS file to include in the analysis. Default is 
       1381.
    """

    images, wcs_list, img_time = load_data(files[start:end+1])
    tab = Table.read(table_file, format='csv', comment='#')

    # Calculate median frame
    med = np.nanmedian(images, axis=0)

    # Define the center of the image
    ny, nx = images[0].shape
    center = (nx // 2, ny // 2)

    # Turn the coordinates of the object into a SkyCoords object
    coords = SkyCoord(ra=tab['ra']*u.hourangle, dec=tab['dec']*u.deg)

    # Define where the output images will be stored
    shifted_images = np.zeros((len(images), size, size))

    # Loop through each FITS file
    i = 0
    for img, wcs, tmid in zip(images, wcs_list, img_time):

        # Find the closest (RA, Dec) to when the image was taken
        idx = np.where( tab['mjd'] >= tmid )[0][0]

        # Get pixel positions
        xpix, ypix = tab['xpixel'][idx], tab['ypixel'][idx]
        xpix = round(xpix)
        ypix = round(ypix)
        
        # Calculate offset from center of the image
        dx = xpix - center[0]# - xpix
        dy = ypix - center[1]# - ypix

        # Shift the image
        shifted = shift(img-med, shift=(dy, dx), order=1, mode='constant',
                        cval=np.nan)

        """
        if i % 20 == 0:
            kwargs = {'vmin':0, 'vmax':300,
                      'aspect':'auto', 'origin':'lower'}
            fig, (ax1, ax2) = plt.subplots(ncols=2, 
                                           figsize=(14,5),
                                           sharex=True, sharey=True)
            ax1.imshow(img, **kwargs) 
            ax1.scatter(ypix, xpix, facecolor=None, 
                     edgecolor='r', marker='o',
                     s = 50)
            ax2.imshow(shifted, **kwargs)
            plt.savefig('./figures/shifted_{0:04d}.png'.format(i),
                        dpi=300, bbox_inches='tight')
            plt.close()
        #"""

        shifted_images[i] = shifted[int(center[1]-size/2) : int(center[1]+size/2),
                                    int(center[0]-size/2) : int(center[0]+size/2)]
    
        i += 1

    return shifted_images


which = '2-3'
path = f'/home/adina/.eleanor/ffis/s0092/{which}/clean'
files = np.sort([os.path.join(path, i) for i in os.listdir(path) if i.endswith('.fits')])

shifty = center_object(files, 
                       'iso_tess_locations_{}.csv'.format(which),
                       start=700, end=1380)

summed_img = np.nansum(shifty, axis=0)
print(np.unique(shifty))

# Make a figure out of the shifted images
fig = plt.figure(figsize=(10,10))
plt.imshow(np.nansum(shifty, axis=0), vmin=-300, vmax=100,
           aspect='auto', origin='lower')
plt.colorbar()
plt.savefig('shifted-{}_v1.png'.format(which), dpi=300,
            bbox_inches='tight')
plt.close()

"""
# Make a movie out of the shifted images
fig = plt.figure(figsize=(10,10))
frames = []

extent = (y-size, y+size+1, x-size, x+size+1)

for i in range(len(cutout)):
    frames.append([plt.imshow(cutout[i]-med, origin='lower',
                              aspect='auto', extent=extent,
                              vmin=-300, vmax=50)])

print('Writing animation')


"""
