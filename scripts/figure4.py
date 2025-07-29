import numpy as np
from astropy import units
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.timeseries import LombScargle

rc = Table.read('../rcParams.txt', format='csv')
for name, val in zip(rc['name'], rc['value']):
    plt.rcParams[name] = val

# Load in the relevant data
data1 = np.load('../data/stacked_3I_2-3.npy', allow_pickle=True).item()
tpf1 = data1['subtracted'] * 1.0

data2 = np.load('../data/stacked_3I_1-2.npy', allow_pickle=True).item()
tpf2 = data2['subtracted'] * 1.0

#################
# Plot the data #
#################

fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(8, 3.5))

ax1.imshow(np.nansum(tpf1, axis=0)[2:,2:], aspect='auto',
                     origin='lower', vmin=4398.5)

q = data2['good_frames'] == 1
ax2.imshow(np.nansum(tpf2[q], axis=0)[2:,2:], aspect='auto',
                     origin='lower')

for ax in [ax1, ax2]:
    ax.plot(8,8,'o', ms=70, color="none", markeredgecolor='#de4f0d', markeredgewidth=3)
    ax.set_xticks([])
    ax.set_yticks([])

    ax.axhline(0.5, 0.07, 0.172, color='w', zorder=10, lw=3)
    ax.text(0.8, 0.8, "41''", color='w')

    ax.set_xlim(-0.5, 16.5)
    ax.set_ylim(-0.5, 16.5)

ax1.set_title('Camera 2 CCD 3', fontsize=16)
ax2.set_title('Camera 1 CCD 2', fontsize=16)

plt.savefig('../figures/3I_deepstack.pdf', bbox_inches='tight', dpi=300)
