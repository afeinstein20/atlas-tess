import numpy as np
from astropy import units
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.timeseries import LombScargle

rc = Table.read('../rcParams.txt', format='csv')
for name, val in zip(rc['name'], rc['value']):
    plt.rcParams[name] = val

# Load in the relevant data
data1 = np.load('../data/stacked_3I_2-3_v3.npy', allow_pickle=True).item()
tpf1 = data1['subtracted'] * 1.0
err1 = data1['err_sub'] * 1.0

data2 = np.load('../data/stacked_3I_1-2_v3.npy', allow_pickle=True).item()
tpf2 = data2['subtracted'] * 1.0
err2 = data2['err_sub'] * 1.0

#################
# Plot the data #
#################

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2, figsize=(8, 10))

im1 = ax1.imshow(np.nansum(tpf1, axis=0), aspect='auto',
                 origin='lower', vmin=4398.5)
im3 = ax3.imshow(np.sqrt(np.nansum(err1**2, axis=0)), aspect='auto', origin='lower',
                 cmap='Greys')

q = data2['good_frames'] == 0
im2 = ax2.imshow(np.nansum(tpf2[q], axis=0), aspect='auto',
                 origin='lower')

im4 = ax4.imshow(np.sqrt(np.nansum(err2[q]**2, axis=0)), aspect='auto',
                 origin='lower', cmap='Greys')

ims = [im1, im2, im3, im4]
letter = ['(a)', '(b)', '(c)', '(d)']

for i, ax in enumerate([ax1, ax2, ax3, ax4]):
    ax.plot(9,9,'o', ms=70, color="none", markeredgecolor='#de4f0d', markeredgewidth=3)
    ax.set_xticks([])
    ax.set_yticks([])
    cbar = plt.colorbar(ims[i], orientation='horizontal')
    if i < 2:
        cbar.set_label('Summed Counts')
        fc = 'w'
    else:
        cbar.set_label(r'$\sigma$ Counts')
        fc = 'k'

    ax.axhline(0.5, 0.06, 0.15, color=fc, zorder=10, lw=3)
    ax.text(0.5, 0.8, "41''", color=fc)
    ax.text(0., 16.5, letter[i], color=fc, fontweight='bold')#backgroundcolor='k')

ax1.set_title('Camera 2 CCD 3', fontsize=16)
ax2.set_title('Camera 1 CCD 2', fontsize=16)

plt.savefig('../figures/3I_deepstack_v2.pdf', bbox_inches='tight', dpi=300)
