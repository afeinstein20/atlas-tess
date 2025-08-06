import numpy as np
from astropy import units
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.timeseries import LombScargle

rc = Table.read('../rcParams.txt', format='csv')
for name, val in zip(rc['name'], rc['value']):
    plt.rcParams[name] = val

# Load in the relevant data
data1 = np.load('../data/stacked_3I_2-3_v4.npy', allow_pickle=True).item()
tpf1 = data1['subtracted'] * 1.0
err1 = data1['err_sub'] * 1.0

data2 = np.load('../data/stacked_3I_1-2_v4.npy', allow_pickle=True).item()
tpf2 = data2['subtracted'] * 1.0
err2 = data2['err_sub'] * 1.0

#################
# Plot the data #
#################

fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(ncols=2, nrows=3, figsize=(8, 16))

q = data1['good_frames'] == 0
im1 = ax1.imshow(np.nanmedian(tpf1[q], axis=0), aspect='auto',
                 origin='lower', vmin=0)
im3 = ax3.imshow(np.nanmean(tpf1[q], axis=0), aspect='auto', origin='lower',
                 cmap='magma', vmin=0, vmax=0.3)
im5 = ax5.imshow(np.nansum(tpf1[q], axis=0), aspect='auto', origin='lower',
                 cmap='Greys_r', vmin=0, vmax=1000)

q = data2['good_frames'] == 0
im2 = ax2.imshow(np.nanmedian(tpf2[q], axis=0), aspect='auto',
                 origin='lower', vmin=0, vmax=0.5)

im4 = ax4.imshow(np.nanmean(tpf2[q], axis=0), aspect='auto',
                 origin='lower', cmap='magma', vmin=1, vmax=3)

im6 = ax6.imshow(np.nansum(tpf2[q], axis=0), aspect='auto', origin='lower',
                 cmap='Greys_r', vmax=13000, vmin=8000)

ims = [im1, im2, im3, im4, im5, im6]
letter = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']

for i, ax in enumerate([ax1, ax2, ax3, ax4, ax5, ax6]):
    ax.plot(9,9,'o', ms=70, color="none", markeredgecolor='w', markeredgewidth=3)
    ax.set_xticks([])
    ax.set_yticks([])
    cbar = plt.colorbar(ims[i], orientation='horizontal')
    if i < 2:
        cbar.set_label('Median Counts s$^{-1}$')
        fc = 'w'
    elif i >= 2 and i < 4:
        cbar.set_label(r'Mean Counts s$^{-1}$')
        fc = 'w'
    else:
        cbar.set_label(r'Summed Counts s$^{-1}$')
        fc = 'w'

    ax.axhline(0.5, 0.06, 0.15, color=fc, zorder=10, lw=3)
    ax.text(0.5, 0.8, "41''", color=fc)
    ax.text(0., 16.5, letter[i], color=fc, fontweight='bold')#backgroundcolor='k')

ax1.set_title('Camera 2 CCD 3', fontsize=16)
ax2.set_title('Camera 1 CCD 2', fontsize=16)

plt.savefig('../figures/3I_deepstack_v2.pdf', bbox_inches='tight', dpi=300)
