import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec

rcParams = Table.read('../rcParams.txt', format='csv')
for key, val in zip(rcParams['name'], rcParams['value']):
    plt.rcParams[key] = val

# Get the locations of 3I/ATLAS
first_locs = Table.read('../data/3I_pixel_locations_2-3.csv', format='csv')
first_locs.sort('file')
_, u = np.unique(first_locs['file'], return_index=True)
first_locs = first_locs[u][:4399]

second_locs = Table.read('../data/3I_pixel_locations_1-2.csv', format='csv')
second_locs.sort('file')
_, u = np.unique(second_locs['file'], return_index=True)
second_locs = second_locs[u][280:]

locs = [first_locs, second_locs]

# Get the observations
raw_first = fits.getdata('../data/tess2025134074247-s0092-2-3-0289-s_ffic.fits')
raw_second = fits.getdata('../data/tess2025148041926-s0092-1-2-0289-s_ffic.fits')

####################################################
# This figure is actually two separate PDF figures #
####################################################

########################################
# Plot the full FFIs oriented properly #
########################################

fig = plt.figure(figsize=(9,6))

gs = GridSpec(2, 2, figure=fig, height_ratios=[1,0.05])

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
cax = fig.add_subplot(gs[1, :])

axes = [ax1, ax2]

pdict = {'vmin':0, 'vmax':2000, 'aspect':'auto', 'cmap':'Greys_r'}

axes[0].imshow(raw_first, origin='upper', **pdict)

axes[0].set_xlim(2136,0)

im1 = axes[1].imshow(raw_second, origin='lower', **pdict)

# Add cutout patches
for i in range(2):
    xmin, xmax = np.nanmin(locs[i]['x']), np.nanmax(locs[i]['x'])
    ymin, ymax = np.nanmin(locs[i]['y']), np.nanmax(locs[i]['y'])

    rect = Rectangle((ymin, xmin), ymax-ymin, xmax-xmin, facecolor=None,
                     edgecolor='#de4f0d', lw=2, fill=False)
    axes[i].add_patch(rect)


axes[0].set_title('Camera 2 CCD 3')
axes[1].set_title('Camera 1 CCD 2')

plt.colorbar(im1, cax=cax, label='Counts', orientation='horizontal')

axes[0].text(x=2103, y=170, s='(a)', color='w', fontweight='bold', backgroundcolor='k')

axes[1].set_yticks([])

axes[0].set_ylabel('X Pixel Row')
axes[0].set_xlabel('Y Pixel Column')
axes[1].set_xlabel('Y Pixel Column')

for ax in axes:
    ax.set_rasterized(True)

plt.subplots_adjust(hspace=0.5, wspace=0.05)

plt.savefig('../figures/ffi_only.pdf', bbox_inches='tight', dpi=300)
plt.close()

###########################
# Plot the full postcards #
###########################
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(9,6))
axes = axes.reshape(-1)

pdict = {'vmin':0, 'vmax':2000, 'aspect':'auto', 'cmap':'Greys_r'}

axes[0].imshow(raw_first, origin='upper', **pdict)

axes[0].set_xlim(2136,0)

im1 = axes[1].imshow(raw_second, origin='lower', **pdict)

# Add cutout patches
for i in range(2):
    xmin, xmax = np.nanmin(locs[i]['x']), np.nanmax(locs[i]['x'])
    ymin, ymax = np.nanmin(locs[i]['y']), np.nanmax(locs[i]['y'])

    if i == 0:
        axes[i].set_xlim(ymax, ymin)
    else:
        axes[i].set_xlim(ymin, ymax)
    axes[i].set_ylim(xmin, xmax)


axes[0].set_title('Camera 2 CCD 3')
axes[1].set_title('Camera 1 CCD 2')

axes[0].text(x=300, y=1630, s='(b)', color='w', fontweight='bold', backgroundcolor='k')
axes[1].text(x=66, y=595, s='(c)', color='w', fontweight='bold', backgroundcolor='k')


axes[0].set_ylabel('X Pixel Row')
axes[1].set_xlabel('Y Pixel Column')
axes[1].set_ylabel('X Pixel Column')

for ax in axes:
    ax.set_rasterized(True)

plt.subplots_adjust(hspace=0.35, wspace=0.05)

plt.savefig('../figures/postcards.pdf', bbox_inches='tight', dpi=300)
