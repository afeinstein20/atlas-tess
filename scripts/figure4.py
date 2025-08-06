import numpy as np
from astropy import units
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.timeseries import LombScargle

rc = Table.read('../rcParams.txt', format='csv')
for name, val in zip(rc['name'], rc['value']):
    plt.rcParams[name] = val

# Load in the relevant data
data = np.load('../data/stacked_A918PE_2-3_v3.npy', allow_pickle=True).item()
tpf = data['subtracted'] * 1.0
error = data['error'] * 1.0
time = (data['time']+2400000.5)-2457000

#################
# Plot the data #
#################

fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=3, figsize=(14,14))

im = ax1.imshow(np.nanmedian(tpf, axis=0), aspect='auto', origin='lower', vmin=0, vmax=60)

ax1.plot(9, 9,'o', ms=70, color="none", markeredgecolor='#de4f0d', markeredgewidth=3)
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_title('896 Sphinx (A918 PE)', fontsize=16)
ax1.axhline(0.5, 0.07, 0.172, color='w', zorder=10, lw=3)
ax1.text(0.8, 0.8, "41''", color='w')
ax1.set_xlabel('Y Pixel Column')
ax1.set_ylabel('X Pixel Column')

q = (time > 3804) & (time < 3813)

mask = np.zeros(tpf[0].shape)
lc = np.nansum(tpf[:,8:11, 8:11], axis=(1,2))[q]
lc_err = np.sqrt(np.nansum(error**2, axis=(1,2)))[q]

ax2.errorbar(time[q], lc[q] / np.nanmedian(lc[q]),
             yerr=lc_err[q]/np.nanmedian(lc[q]), marker='.',
             linestyle='', color='k', alpha=0.1)

freq, power = LombScargle(time*units.day, lc).autopower(minimum_frequency=1.0/(50*units.hour),
                                                        maximum_frequency=1.0/(5.0*units.hour))


ax3.plot(1.0/freq, power, color='k')
ax3.axvline(21.038, color=parula[60], lw=3, label=r'$P_{rot} = 21.038 \pm 0.008$ hours;'+
                                                   '\nPolakis (2018)')

ax3.axvline(21.223, color=parula[160], lw=3, label=r'$P_{rot} = 21.223 \pm 0.646$ hours;'+
                                                   '\nThis work')

ax2.set_xlabel('Time [BJD - 2457000]', fontsize=16)
ax2.set_ylabel('Normalized Flux', fontsize=16)
ax3.set_ylabel('Power', fontsize=16)
ax3.set_xlabel('Period [hours]', fontsize=16)

ax3.set_xlim(1,50)
ax3.set_ylim(0,0.05)
ax3.legend()

plt.subplots_adjust(wspace=0.2, hspace=0.4)

for ax in [ax1, ax2, ax3]:
    ax.set_rasterized(True)

plt.savefig('../figures/asteroid_recovery_v2.pdf', bbox_inches='tight', dpi=300)
