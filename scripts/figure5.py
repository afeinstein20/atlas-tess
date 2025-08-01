import numpy as np
from astropy import units
import matplotlib.pyplot as plt
from astropy.table import Table
from scipy.signal import medfilt
from matplotlib.gridspec import GridSpec
from astropy.timeseries import LombScargle
from lightkurve.lightcurve import LightCurve as LC


rc = Table.read('../rcParams.txt', format='csv')
for name, val in zip(rc['name'], rc['value']):
    plt.rcParams[name] = val

first = np.load('../data/stacked_3I_2-3_v3.npy', allow_pickle=True).item()
second= np.load('../data/stacked_3I_1-2_v3.npy', allow_pickle=True).item()

fig = plt.figure(figsize=(12,7))

gs = GridSpec(2, 2, figure=fig, height_ratios=[1,0.6])
ax1 = fig.add_subplot(gs[0,:])
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[1, 1])

axes = [ax2, ax3]
colors = ['#005f60', '#249ea0']

for i, ccd in enumerate([first, second]):

    good = ccd['good_frames'] == 0

    lc = ccd['subtracted'][:,9,9]
    lc_err = ccd['err_sub'][:,9,9]/100
    time = ccd['time'] + 2400000.5 - 2457000

    if i == 1:
        q = (time >= 3818) & (time < 3827.8)
    else:
        q = time > 3804

    print(time[q])
    ax1.errorbar(time[q], lc[q], yerr=lc_err[q],
                 marker='.', linestyle='', alpha=0.4, color=colors[i],
                 zorder=0)

    # binned to 36 minutes
    lk = LC(time=time[q]*units.day, flux=lc[q],
            flux_err=lc_err[q]).bin(time_bin_size=36*units.min)
    ax1.errorbar(lk.time.value, lk.flux.value, yerr=lk.flux_err.value,
                 marker='o', color='k', linestyle='', alpha=0.7)

    min_freq = 1.0/(70.0*units.hour)
    max_freq = 1.0/(1.0*units.hour)

    ls = lk.to_periodogram(minimum_frequency=min_freq,
                           maximum_frequency=max_freq)

    axes[i].plot(ls.period.to(units.hour), ls.power, color='k')
    axes[i].set_xlabel('Period [hours]', fontsize=16)
    axes[i].set_xlim(1,70)



ax1.set_xlabel('Time [BJD - 2457000]', fontsize=16)

ax1.set_ylabel('TESS Magnitude', fontsize=16)
#ax1.set_ylim(20.0, 19.2)
ax1.set_xlim(3803.5, 3828.5)

ax2.set_ylabel('Power', fontsize=16)

ax2.set_ylim(0,0.0006)
ax3.set_ylim(0,0.0006)
ax3.set_yticklabels([])

plt.subplots_adjust(hspace=0.3)

plt.savefig('../figures/tess_lightcurve_v2.pdf', dpi=300, bbox_inches='tight')
