import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

rcParams = Table.read('../rcParams.txt', format='csv')
for key, val in zip(rcParams['name'], rcParams['value']):
    plt.rcParams[key] = val

# Load in the relevant data
data1 = np.load('../data/stacked_3I_2-3_v4.npy', allow_pickle=True).item()
tpf23 = data1['raw'] * 1.0
good23 = data1['good_frames'] + 0

data2 = np.load('../data/stacked_3I_1-2_v4.npy', allow_pickle=True).item()
tpf12 = data2['raw'] * 1.0
good12 = data2['good_frames'] + 0

# Plotting mask
mask = np.zeros(tpf12[0].shape)
mask[4:15, 4:15] = 1

fig, axes = plt.subplots(ncols=3, nrows=2, figsize=(8,8),
                         gridspec_kw={'width_ratios':[1,1,0.1]})
axes = axes.reshape(-1)

## Plot for Camera 2 CCD 3 ##
im = axes[0].imshow(tpf23[good23==0][300], aspect='auto', vmin=100, vmax=1000, origin='lower')
axes[0].set_ylabel('X Pixel Row')

axes[1].imshow(tpf23[good23==1][100], aspect='auto', vmin=100, vmax=1000, origin='lower')

axes[0].text(x=0., y=17, s='(a)', color='w', fontweight='bold')
axes[1].text(x=0., y=17, s='(b)', color='w', fontweight='bold')

plt.colorbar(im, cax=axes[2], label='Counts s$^{-1}$')

## Plot for Camera 1 CCD 2 ##
im = axes[3].imshow(tpf12[good12==0][650], aspect='auto', vmin=100, vmax=5000, origin='lower')
axes[3].set_ylabel('X Pixel Row')

axes[4].imshow(tpf12[good12==1][100], aspect='auto', vmin=100, vmax=5000, origin='lower')

axes[3].text(x=0., y=17, s='(c)', color='w', fontweight='bold')
axes[4].text(x=0., y=17, s='(d)', color='w', fontweight='bold')

plt.colorbar(im, cax=axes[5], label='Counts s$^{-1}$')

# Add the patches
for i in range(len(axes)):
    axes[i].imshow(mask, aspect='auto', cmap='Greys_r', alpha=0.2)
    patch = Rectangle((3.5,3.5), 11, 11, facecolor='none', linewidth=3, edgecolor='r')
    axes[i].add_patch(patch)

    if i != 2 and i != 5:
        axes[i].set_yticks([])
        axes[i].set_xticks([])

    if i == 3 or i == 4:
        axes[i].set_xlabel('Y Pixel Column')

    axes[i].set_rasterized(True)



plt.savefig('../figures/good_bad_together.pdf', bbox_inches='tight', dpi=300)
