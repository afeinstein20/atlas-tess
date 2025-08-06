import numpy as np
from astropy import units
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from astropy.timeseries import LombScargle
from lightkurve.lightcurve import LightCurve


# Set the rcParams
rc = Table.read('../rcParams.txt', format='csv')
for name, val in zip(rc['name'], rc['value']):
    plt.rcParams[name] = val

# Load in the data
first = np.load('../data/stacked_3I_2-3_v4.npy', allow_pickle=True).item()
second= np.load('../data/stacked_3I_1-2_v4.npy', allow_pickle=True).item()

# Camera 2 CCD 3 data
time1 = first['time'] + 2400000.5 - 2457000
q = (time1 > 3804) & (first['good_frames'] == 0)
time1 = time1[q]
lc1 = first['subtracted'][:,8:11,8:11][q]
bkg1 = first['subtracted'][:,3:6,3:6][q]

# Camera 1 CCD 2 data
time2 = second['time'] + 2400000.5 - 2457000
q = (time2 > 3818) & (time2 < 3828) & (second['good_frames'] == 0)

time2 = time2[q]
lc2 = second['subtracted'][:,8:11,8:11][q]
bkg2 = second['subtracted'][:,3:6,3:6][q]

# Periodogram function
def periodogram(t, f, maxp=1*units.hour, minp=70*units.hour):
    lk = LightCurve(time=t*units.day, flux=f+1.0).bin(time_bin_size=36*units.minute)
    ls = lk.to_periodogram(minimum_frequency=1.0/minp,
                           maximum_frequency=1.0/maxp)

    return ls, lk

#########################
#### Make the Figure ####
#########################
fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(14,5), gridspec_kw={'width_ratios':[0.23,1]})
fig.set_facecolor('w')
axes = axes.reshape(-1)

dc = '#de4f0d'
bc = '#ffbe35'

q = (first['good_frames'] == 0)

# Imshow the TPF
axes[0].imshow(np.nanmedian(first['subtracted'][q], axis=0), aspect='auto', cmap='Greys_r')

# Plot the apertures
rect = Rectangle((7.25, 7.25), 3.45, 3.45, edgecolor=dc, facecolor='none', lw=3)
axes[0].add_patch(rect)

rect = Rectangle((2.25, 2.25), 3.45, 3.45, edgecolor=bc, facecolor='none', lw=3)
axes[0].add_patch(rect)


# Make the Lomb-Scargle Periodograms
bkg_23, bk23 = periodogram(time1, bkg1)
axes[1].plot((bkg_23.period).to(units.hour), bkg_23.power, color=bc, lw=2)

ls_23, m23 = periodogram(time1, lc1)
axes[1].plot((ls_23.period).to(units.hour), ls_23.power, color=dc, lw=2)

##################
# Camera 1 CCD 2 #
##################

q = second['good_frames'] == 0 #& (lc2 < 20.) & (bkg2 > 19.5)

# Plot the TPF
axes[2].imshow(np.nanmedian(second['subtracted'][q], axis=0),
               aspect='auto',
               cmap='Greys_r')

# Plot the apertures
rect = Rectangle((7.25, 7.25), 3.45, 3.45, edgecolor=dc, facecolor='none', lw=3)
axes[2].add_patch(rect)

rect = Rectangle((2.25, 2.25), 3.45, 3.45, edgecolor=bc, facecolor='none', lw=3)
axes[2].add_patch(rect)

# Make the Lomb-Scargle Periorograms
bkg_12, mb12 = periodogram(time2, bkg2)
axes[3].plot((bkg_12.period).to(units.hour), bkg_12.power, color=bc, lw=2)

ls_12, m12 = periodogram(time2, lc2)
axes[3].plot((ls_12.period).to(units.hour), ls_12.power, color=dc, lw=2)


axes[3].set_xlabel('Period [hours]')

for i in [1,3]:
    axes[i].set_xlim(1,70)
    axes[i].set_ylabel('Power')
    if i == 1:
        label = '16.16 hours'
    else:
        label=''
    axes[i].axvline(16.16, color='k', linestyle=':', label=label)

axes[1].legend()

for i in [0,2]:
    axes[i].set_xticks([])
    axes[i].set_yticks([])
#axes[3].set_ylim(0,0.06)
#axes[1].set_ylim(0,0.04)
axes[1].set_xticklabels([])

axes[0].set_ylabel('Camera 2 CCD 3')
axes[2].set_ylabel('Camera 1 CCD 2')

for i in range(4):
    axes[i].set_rasterized(True)

plt.subplots_adjust(hspace=0.3, wspace=0.25)
plt.savefig('../figures/ls_comparison_v2.pdf', dpi=300, bbox_inches='tight')
