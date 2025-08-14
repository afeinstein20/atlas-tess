import numpy as np
from astropy import units
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.timeseries import LombScargle

rc = Table.read('../rcParams.txt', format='csv')
for name, val in zip(rc['name'], rc['value']):
    plt.rcParams[name] = val

plt.rcParams['font.size'] = 18

parula = np.load('../data/parula_colors.npy')

# Load in the relevant data
data = np.load('../data/stacked_A918PE_2-3_v3.npy', allow_pickle=True).item()
tpf = data['subtracted'] * 1.0
error = data['error'] * 1.0
time = (data['time']+2400000.5)-2457000

#################
# Plot the data #
#################

fig, (ax1, ax2, ax3) = plt.subplots(ncols=1, nrows=3, figsize=(12,14))

im = ax1.imshow(np.nanmedian(tpf, axis=0), origin='lower', vmin=0, vmax=60)

ax1.plot(9, 9,'o', ms=70, color="none", markeredgecolor='w', markeredgewidth=3)
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_title('896 Sphinx (A918 PE)\nMedian Deepstack Image', fontsize=16)
ax1.axhline(0.5, 0.07, 0.172, color='w', zorder=10, lw=3)
ax1.text(0.8, 0.8, "41''", color='w')
ax1.set_xlabel('Y Pixel Column')
ax1.set_ylabel('X Pixel Row')
ax1.text(x=0, y=16.5, s='(a)', fontweight='bold', color='w')

q = (time > 3804) & (time < 3813)

mask = np.zeros(tpf[0].shape)
lc = np.nansum(tpf[:,8:11, 8:11], axis=(1,2))[q]
lc_err = np.sqrt(np.nansum(error**2, axis=(1,2)))[q]

ax2.errorbar(time[q], lc / np.nanmedian(lc),
             yerr=lc_err/np.nanmedian(lc)*100, marker='.',
             linestyle='', color='k', alpha=0.1)

ax2.errorbar(time[q], lc / np.nanmedian(lc), marker='.',
             linestyle='', color='k')

ax2.text(x=3803.6, y=1.5, s='(b)', fontweight='bold')
ax2.set_title('Light Curve of 896 Sphinx')

freq, power = LombScargle(time[q]*units.day, lc).autopower(minimum_frequency=1.0/(50*units.hour),
                                                        maximum_frequency=1.0/(1.0*units.hour))


ax3.plot(1.0/freq, power, color='k')
ax3.axvline(21.038, color=parula[60], lw=3, label=r'$P_{rot} = 21.038 \pm 0.008$ hours;'+
                                                   '\nPolakis (2018)')

ax3.axvline(21.223, color=parula[160], lw=3, label=r'$P_{rot} = 21.223 \pm 0.646$ hours;'+
                                                   '\nThis work')

ax2.set_xlabel('Time [BJD - 2457000]')
ax2.set_ylabel('Normalized Flux')
ax3.set_ylabel('Power')
ax3.set_xlabel('Period [hours]')

ax3.set_xlim(1,50)
ax3.set_ylim(0,0.05)
leg = ax3.legend()
for legobj in leg.get_lines():
    legobj.set_linewidth(6.0)

ax3.text(x=1.5, y=0.045, s='(c)', fontweight='bold')

plt.subplots_adjust(wspace=0.2, hspace=0.4)

for ax in [ax1, ax2, ax3]:
    ax.set_rasterized(True)


plt.savefig('../figures/asteroid_recovery_v2.pdf', bbox_inches='tight', dpi=300)
