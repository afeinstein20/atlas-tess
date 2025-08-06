import numpy as np
from astropy import units
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.timeseries import LombScargle
from lightkurve.lightcurve import LightCurve

parula = np.load('../data/parula_data.npy')

# Periodogram function
def periodogram(t, f, maxp=1*units.hour, minp=70*units.hour):
    lk = LightCurve(time=t*units.day, flux=f+1.0).bin(time_bin_size=36*units.minute)
    ls = lk.to_periodogram(minimum_frequency=1.0/minp,
                           maximum_frequency=1.0/maxp)

    return ls, lk

# Set the rcParams
rc = Table.read('../rcParams.txt', format='csv')
for name, val in zip(rc['name'], rc['value']):
    plt.rcParams[name] = val

# Load in the data
first = np.load('../data/stacked_3I_2-3_v4.npy', allow_pickle=True).item()
second= np.load('../data/stacked_3I_1-2_v4.npy', allow_pickle=True).item()

time1 = first['time'] + 2400000.5 - 2457000
time2 = second['time'] + 2400000.5 - 2457000

fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(14,8), gridspec_kw={'width_ratios':[0.4,1]})
axes = axes.reshape(-1)

# Plot median deepstack image for Camera 2 CCD 3
q = first['good_frames'] == 0
im1 = axes[0].imshow(np.nanmedian(first['subtracted'][q], axis=0), aspect='auto',
                 origin='lower', vmin=0, cmap='Greys_r')

q = (time1 > 3804) & (first['good_frames'] == 0)

x = 1
for i in range(7, 12):
    lc = first['subtracted'][:,i, i]
    p = periodogram(time1[q], lc[q])
    if i != 9:
        color, lw, zorder = parula[int(x*14)], 2, 0
        x += 1
    else:
        color, lw, zorder = 'k', 3, 100

    axes[1].plot(p[0].period.to(units.hour), p[0].power, color=color, lw=lw,
                 zorder=zorder)
    axes[0].plot(i, i, 'x', color=color, zorder=100, markeredgewidth=2)


# Plot median deepstack image for Camera 1 CCD 2
q = second['good_frames'] == 0
im2 = axes[2].imshow(np.nanmedian(second['subtracted'][q], axis=0), aspect='auto',
                     origin='lower', vmin=0, vmax=0.5, cmap='Greys_r')

q = (time2 > 3818) & (time2 < 3828) & (second['good_frames'] == 0)

x = 1
for i in range(7, 12):
    lc = second['subtracted'][:,i, i]
    p = periodogram(time2[q], lc[q])
    if i != 9:
        color, lw, zorder = parula[int(x*14)], 2, 0
        x += 1
    else:
        color, lw, zorder = 'k', 3, 100

    axes[3].plot(p[0].period.to(units.hour), p[0].power, color=color, lw=lw,
                 zorder=zorder)
    axes[2].plot(i, i, 'x', color=color, zorder=100, markeredgewidth=2)

axes[3].set_xlabel('Period [hours]')
for i in [1,3]:
    axes[i].set_ylabel('Power')
    axes[i].set_xlim(1,70)
    axes[i].axvline(16.16, color='k', linestyle=':', zorder=0, label='16.16 hours')
    axes[1].legend()

for i in [0,2]:
    axes[i].set_yticks([])
    axes[i].set_xticks([])
    axes[i].set_ylabel('X Pixel Row')
axes[2].set_xlabel('Y Pixel Row')

plt.savefig('../figures/trail_lead.pdf', bbox_inches='tight', dpi=300)
