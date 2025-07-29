import numpy as np
from astropy import units
import matplotlib.pyplot as plt
from astropy.table import Table
from scipy.signal import medfilt
from matplotlib.gridspec import GridSpec
from astropy.timeseries import LombScargle


rc = Table.read('../rcParams.txt', format='csv')
for name, val in zip(rc['name'], rc['value']):
    plt.rcParams[name] = val

def convert_mag(counts):
    # Converts counts/second to magnitude
    # From the TESS handbook
    dur = (len(counts) * 200 * units.s)
    dur = dur.value
    return -2.5 * np.log10(counts/200) + 20.44

first = np.load('../data/stacked_3I_2-3.npy', allow_pickle=True).item()
second= np.load('../data/stacked_3I_1-2.npy', allow_pickle=True).item()

fig = plt.figure(figsize=(12,7))

gs = GridSpec(2, 2, figure=fig, height_ratios=[1,0.6])
ax1 = fig.add_subplot(gs[0,:])
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[1, 1])

axes = [ax2, ax3]
colors = ['#005f60', '#249ea0']

for i, ccd in enumerate([first, second]):
    lc = ccd['subtracted'][:,10,10] * np.nanmedian(convert_mag(ccd['raw'][:,10,10]))
    time = ccd['time'] + 2400000.5 - 2457000

    if i == 1:
        q = (time > 3817.3) & (time < 3827.5)
    else:
        q = time > 3800

    ax1.plot(time[q], lc[q], color=colors[i])

    m = medfilt(lc[q], 11)
    ax1.plot(time[q], m, 'k')

    frequency, power = LombScargle(time[q]*units.day, m).autopower(minimum_frequency=1.0/(120.0*units.hour),
                                                                       maximum_frequency=1.0/(1.0*units.hour))
    axes[i].plot(1.0/frequency, power, color='k')
    axes[i].set_xscale('log')
    axes[i].set_xlabel('Period [hours]', fontsize=16)


ax1.set_xlabel('Time [BJD - 2457000]', fontsize=16)

ax1.set_ylabel('TESS Magnitude', fontsize=16)
ax1.set_ylim(20.4, 19.4)
ax1.set_xlim(3802.5, 3828.5)

ax2.set_ylabel('Power', fontsize=16)

ax2.set_ylim(0,0.035)
ax3.set_ylim(0,0.035)

plt.subplots_adjust(hspace=0.3)

plt.savefig('../figures/tess_lightcurve.pdf', dpi=300, bbox_inches='tight')
