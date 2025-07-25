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

which = '2-3'
#which = '1-2'
path = f'/home/adina/.eleanor/ffis/s0092/{which}/clean_64'
tab = Table.read('iso_tess_locations_{}.csv'.format(which), format='csv')

tab = tab[(tab['xpixel'] > 0) & (tab['ypixel'] > 0)]

files = np.sort([os.path.join(path, i) for i in os.listdir(path) if i.endswith('.fits')])

start = 0
end = len(files)

size = 150
x = int(np.nanmedian(tab['xpixel']))
y = int(np.nanmedian(tab['ypixel']))

#x, y = 1024, 1024

cutout = np.zeros((end-start, size*2+1, size*2+1))
#error  = np.zeros((end-start, size*2+1, size*2+1))
#times  = np.zeros((end-start, 3))

bad = np.array([], dtype='U300')

for i in tqdm(range(start, end)):

#    try:
    hdu = fits.open(files[i])

        #hdu[1].data = hdu[1].data[y-size:y+size+1, x-size:x+size+1] + 0.0
        #hdu[2].data = hdu[2].data[y-size:y+size+1, x-size:x+size+1] + 0.0

        #newpath = path+'/cutout/'+files[i].split('/')[-1]
        #print(newpath)
        #hdu.writeto(newpath, overwrite=True)

    cutout[i] = hdu[0].data[y-size:y+size+1, x-size:x+size+1]

    hdu.close()


    #except:
    #    bad = np.append(bad, files[i])

    #os.system('rm {}'.format(files[i]))


med = np.nanmedian(cutout[:250,:,:], axis=0)

fig = plt.figure(figsize=(10,10))
frames = []

extent = (y-size, y+size+1, x-size, x+size+1)

for i in range(len(cutout)):
    frames.append([plt.imshow(cutout[i]-med, origin='lower',
                              aspect='auto', extent=extent, 
                              vmin=-300, vmax=50)])

plt.plot(180, 1600, 'rx')

print('Writing animation')
ani = animation.ArtistAnimation(fig, frames, interval=50, blit=True,
                                repeat_delay=1000)
ani.save("frame-ffi-{}_v2.gif".format(which))

