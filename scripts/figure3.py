import numpy as np
from astropy import units
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.timeseries import LombScargle

rc = Table.read('../rcParams.txt', format='csv')
for name, val in zip(rc['name'], rc['value']):
    plt.rcParams[name] = val

# Load in the relevant data
data = np.load('../data/stacked_A918PE_2-3.npy', allow_pickle=True).item()
tpf = data['subtracted'] * 1.0
time = (data['time']+2400000.5)-2457000

#################
# Plot the data #
#################

fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(14,3),
                               gridspec_kw={'width_ratios':[0.3,1]})

im = ax1.imshow(np.nansum(tpf, axis=0)[2:,2:], aspect='auto',
                origin='lower', vmin=4400, vmax=5000)
ax1.plot(8,8,'ro', ms=70, color="none", markeredgecolor='r', markeredgewidth=3)
ax1.set_xticks([])
ax1.set_yticks([])
ax2.set_title('896 Sphinx (A918 PE)', fontsize=16)
ax1.axhline(0.5, 0.07, 0.172, color='w', zorder=10, lw=3)
ax1.text(0.8, 0.8, "41''", color='w')

ax1.set_xlim(-0.5, 16.5)
ax1.set_ylim(-0.5, 16.5)

mask = np.zeros(tpf[0].shape)
mask[9:12, 9:12] = 1.0
lc = np.nansum(tpf*mask, axis=(1,2))
ax2.plot(time, lc / np.nanmedian(lc), color='k')

q = time > 3804
freq, power = LombScargle(time*units.day, tpf[:,10,10]).autopower(minimum_frequency=1.0/(50*units.hour),
                                                                        maximum_frequency=1.0/(5.0*units.hour))


ax2.set_xlabel('Time [BJD - 2457000]', fontsize=16)
ax2.set_ylabel('Normalized Flux', fontsize=16)
plt.subplots_adjust(wspace=0.2)

plt.savefig('../figures/asteroid_recovery.pdf', bbox_inches='tight', dpi=300)
