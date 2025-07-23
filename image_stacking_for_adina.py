#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 12 04:19:22 2025

@author: jwn0027
"""
from astropy.io import fits
from astroquery.jplhorizons import Horizons
from astropy.wcs import WCS
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import Angle,SkyCoord
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from glob import glob
from matplotlib import colors
from astropy.table import Table
#from wfc3tools import calwf3
import datetime
from matplotlib.patches import Ellipse
from time import sleep
import os
from reproject import reproject_interp, reproject_exact
from reproject.mosaicking import reproject_and_coadd
from reproject.mosaicking import find_optimal_celestial_wcs
import scipy.ndimage as ndimage
import scipy.signal as signal

from DIA.clean_function import clean_image_DIA

which = '2-3'
location = f'/home/adina/.eleanor/ffis/s0092/{which}/aligned'
record_id = '3I'
MEDIAN_FFI = np.load(f'/home/adina/.eleanor/ffis/s0092/median_{which}_aligned_all.npy')

def clean_image_noonan(data_slice):
    '''
    Parameters
    ----------
    data_slice : 2D image array or IFU slice

    Returns
    -------
    clean_slice : 2D image array or IFU slice with cosmic rays remove

    '''
    #Takes the IFU slice or 2D image array and the data quality slice to produce a masked array
    #Identify outliers

    clean_mask = np.ones_like(data_slice)
    rows,cols = data_slice.shape[0],data_slice.shape[1]
    clean_slice = data_slice
    #plt.imshow(clean_slice)
    for row in range(0,rows):
        for col in range(0,cols):
            perimeter_median = np.median(data_slice[row-1:row+1,col-1:col+1])
            if abs(data_slice[row,col]) >= abs(15*perimeter_median):
                clean_slice[row,col] = perimeter_median
                if abs(data_slice[row,col]) >= abs(25*perimeter_median):
                    #print('Flagged value: '+str(data_slice[row,col]))
                    #print('Replaced with: '+str(perimeter_median))
                    clean_mask[row,col]= 0
                    clean_slice[row,col] = perimeter_median
    return clean_slice

def find_image_shift(image1, image2):
    """
    Compute the best-fit alignment shift between two images using cross-correlation.
    :param image1: First image
    :param image2: Second image
    :return: (dx, dy) shift required to align image2 to image1
    """
    correlation = signal.correlate2d(image1, image2, boundary='symm', mode='same')
    y_max, x_max = np.unravel_index(np.argmax(correlation), correlation.shape)
    center_y, center_x = np.array(image1.shape) // 2
    dy = y_max - center_y
    dx = x_max - center_x
    return dx, dy

def get_ephemeris():
    """
    Query the JPL Horizons Small Body catalog to get the ephemeris for the object in
    question.
    """
    output_name = f"JPL_horizons_{which}_{'-'.join(e for e in record_id)}.csv"

    all_files = np.sort([os.path.join(location, fn) for fn in os.listdir(location) if fn.endswith('.fits')])
    all_files = all_files[500:]

    tab = Table(names=['filename', 'date', 'ra', 'dec', 'x', 'y'], dtype=[str, float, float, float, float, float])

    for i, obs in enumerate(all_files):
        print(i)
        hdu = fits.open(obs)
        header = hdu[1].header
        wcs = WCS(header)

        obs_date = Time(header['DATE-OBS'], scale='utc').mjd
        obs_end  = Time(header['DATE-END'], scale='utc').mjd

        #define the observing time for the middle point of the observations
        obs_obj = Time((obs_date+obs_end)/2,format='mjd',scale='utc')
        obstime_jd = obs_obj.jd

        #Query Horizons for the orbit
        obj = Horizons(id=record_id,id_type='smallbody',location='@tess',epochs=obstime_jd)
        eph = obj.ephemerides()

        #Create a SkyCoord for the object at the midpoint of the observations
        c = SkyCoord(eph['RA'][0], eph['DEC'][0], frame='icrs', unit='deg')

        cen_skycoord = SkyCoord(c.ra,c.dec,frame='icrs', unit=(u.deg, u.deg))

        # convert the objects SkyCoord to the frame pixels
        pix_cen = wcs.world_to_pixel(cen_skycoord)

        tab.add_row([obs, obstime_jd, eph['RA'][0], eph['DEC'][0], pix_cen[0], pix_cen[1]])

        hdu.close()

        tab.write(output_name, format='csv', overwrite=True)



def crop_and_coadd(method = 'mean',best_fit=False,use_full=False,verbose = False,save_int = True):
    c_skycoord_list = []
    err_list = []
    images = []

    if method == 'median':
        med_bool = True
    else:
        med_bool = False

    all_files = np.sort([os.path.join(location, fn) for fn in os.listdir(location) if fn.endswith('.fits')])

    if best_fit == True:
        #Pick a file with a WCS solution you trust here.
        comp_data= fits.getdata(all_files[2000])

    for i,obs in enumerate(all_files):

        print(i, obs)

        if 'DS_Store' in obs:
            continue

        ext = 1
        hdu = fits.open(obs)
        header = hdu[ext].header
        
        data = hdu[ext].data - MEDIAN_FFI
    
        mask = data == 0
        data[mask] = np.nan


        #Access the wcs_header
        #wcs_header = fits.getheader(obs, 1)
        #data = fits.getdata(obs)
        drz_wcs1 = WCS(header)

        #These may be different depending on the datasets. Double check
        #"""
        obs_date = Time(header['DATE-OBS'], scale='utc').mjd
        obs_end  = Time(header['DATE-END'], scale='utc').mjd

        #define the observing time for the middle point of the observations
        obs_obj = Time((obs_date+obs_end)/2,format='mjd',scale='utc')
        obstime_jd = obs_obj.jd

        #Query Horizons for the orbit
        obj = Horizons(id=record_id,id_type='smallbody',location='@tess',epochs=obstime_jd)
        eph = obj.ephemerides()

        #Create a SkyCoord for the object at the midpoint of the observations
        c = SkyCoord(eph['RA'][0], eph['DEC'][0], frame='icrs', unit='deg')
        c_skycoord_list.append(c)
        err_list.append([eph['RA_3sigma'][0],eph['DEC_3sigma'][0]])
        cen_skycoord = SkyCoord(c.ra,c.dec,frame='icrs', unit=(u.deg, u.deg))

        # convert the objects SkyCoord to the frame pixels
        pix_cen = drz_wcs1.world_to_pixel(cen_skycoord)
        pix_cen = [pix_cen[0], pix_cen[1]]
        #"""

        #data = clean_image_DIA(obs, xcen=int(pix_cen[0]), ycen=int(pix_cen[1]), bxs=256)

        if ((pix_cen[0] > 0) & (pix_cen[0] < 2136) &
            (pix_cen[1] > 0) & (pix_cen[1] < 2136)):

            #Clean the data if it needs it (you can also use LACosmic or somethin')
            data_clean = clean_image_noonan(data)

            if best_fit == True and i >= 1:
                #Perform iterative fitting to find best small shift to data
                if verbose == True:
                    print('Starting alignment process for: '+str(obs))

                dx,dy = find_image_shift(data_clean, images)
                pix_cen = [pix_cen[0]+dx,pix_cen[1]+dy]

            else:
                pix_cen = drz_wcs1.world_to_pixel(cen_skycoord)

            shape = 50
            x_low,x_high = int(pix_cen[0])-shape, int(pix_cen[0])+shape
            y_low,y_high = int(pix_cen[1])-shape, int(pix_cen[1])+shape

            data_trim = data_clean[x_low:x_high,y_low:y_high]

            images.append(data_trim)

            shape_small = 50

            if i == 0:
                data_og = data
                data_sum = data_trim
                #save a larger image for the fitting steps, centered on the best fit location for V6
                x_low,x_high = int(pix_cen[0])-shape_small, int(pix_cen[0])+shape_small
                y_low,y_high = int(pix_cen[1])-shape_small, int(pix_cen[1])+shape_small
                data_trim = data_clean[x_low:x_high,y_low:y_high]
                comp_data = data_trim
            else:
                data_og += data
                data_sum += data_trim
                if save_int == True:
                    output = np.median(np.array(images), axis = 0)
                    hdu_coadd = fits.PrimaryHDU(output)
                    hdu_list = fits.HDUList([hdu_coadd])
                    hdu_list.writeto(f'3I/obj_img_stacked_{which}_{i:04d}.fits',overwrite=True)
                    hdu_list.close()

                    os.system('rm 3I/obj_img_stacked_{0}_{1:04d}.fits'.format(which, i-1))
                if use_full == True:
                    x_low,x_high = int(pix_cen[0])-shape_small, int(pix_cen[0])+shape_small
                    y_low,y_high = int(pix_cen[1])-shape_small, int(pix_cen[1])+shape_small
                    data_trim = data_clean[x_low:x_high,y_low:y_high]

                    comp_data = data_trim

        hdu.close()
        del data
        del header

    #Save new stack to FITS file
    hdu_coadd = fits.PrimaryHDU(data_sum/i)
    hdu_list = fits.HDUList([hdu_coadd])
    hdu_list.writeto(f'3I/obj_img_stacked_{which}.fits',overwrite=True)

    data_med_array = np.array(images)

    if method == 'median':
        median = np.median(data_med_array, axis = 0)
        hdu_median = fits.PrimaryHDU(median)
        hdu_list = fits.HDUList([hdu_median])
        hdu_list.writeto('obj_img_median_asteroid.fits',overwrite=True)
        plt.figure(4)
        plt.imshow(median, vmin=0, vmax=1,cmap='viridis')
        plt.savefig('obj_img_median_asteroid.png', dpi=300, bbox_inches='tight')
        return median

    elif method == 'mean':
        mean = np.mean(data_med_array, axis = 0)
        hdu_mean = fits.PrimaryHDU(mean)
        hdu_list = fits.HDUList([hdu_mean])
        hdu_list.writeto('obj_img_mean_asteroid.fits',overwrite=True)
        plt.figure(4)
        plt.imshow(mean, vmin=0, vmax=1,cmap='viridis')
        plt.savefig('obj_img_mean_asteroid.png', dpi=300, bbox_inches='tight')
        return mean


crop_and_coadd(save_int=True)
#get_ephemeris()
