import os
import sys
from astropy import units
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

tab = Table.read('../rcParams.txt', format='csv')

for key, val in zip(tab['name'], tab['value']):
    plt.rcParams[key] = val

dat1 = np.load('../data/avg_2-3.npy')
time1 = dat1[:,0] - 2457000

dat2 = np.load('../data/avg_1-2.npy')
time2 = dat2[:,0] - 2457000

fig, ax = plt.subplots(figsize=(10,4))
fig.set_facecolor('w')

dia_ignored = '#fd5901'
ignored = '#f6cb52'
lw = 0
alpha = 0.5

plt.semilogy(time1, dat1[:,1], '.', label='Camera 2 CCD 3', color='#005f60')
plt.axvspan(time1[0]+10.4, time1[-1]+0.1, alpha=alpha, lw=lw, color=ignored)

plt.semilogy(time2, dat2[:,1], '.', label='Camera 1 CCD 2', color='#249ea0')
plt.axvspan(time2[0]+14.5, time2[-1]+0.1, alpha=alpha, lw=lw, color=ignored)
plt.axvspan(time2[0]-0.1, time2[0]+1, alpha=alpha, lw=lw, color=ignored, label='Ignored')

plt.ylabel('Median Counts per FFI')
plt.xlabel('Time [BJD - 2457000]')
plt.legend(markerscale=3, ncol=3, bbox_to_anchor=(0,1.02,1,0.102), loc=3, mode="expand", borderaxespad=0)

plt.xlim(3802.5, 3830)
ax.set_rasterized(True)

plt.savefig('../figures/ffi_cutoffs.png', bbox_inches='tight', dpi=300)
