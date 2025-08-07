import numpy as np
from astropy import units
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.gridspec import GridSpec
from astropy.timeseries import LombScargle
from lightkurve.lightcurve import LightCurve


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

fig = plt.figure(figsize=(14,4))
fig.set_facecolor('w')

colors = ['#005f60', '#249ea0']

for i, ccd in enumerate([first, second]):
    lc = ccd['subtracted'][:,10,10] * np.nanmedian(convert_mag(ccd['raw'][:,10,10]))
    time = ccd['time'] + 2400000.5 - 2457000

    if i == 1:
        q = (time > 3816.9) & (time < 3828.2)
    else:
        q = time > 3800

    plt.plot(time[q], lc[q], '.', color=colors[i])

    binned = LightCurve(time=time[q] * units.day, flux=lc[q]).bin(1*units.hour)
    plt.errorbar(binned.time.value, binned.flux.value,
                 yerr=binned.flux_err.value, linestyle='',
                 marker='o', markeredgecolor='k',
                 color='none')


plt.xlabel('Time [BJD - 2457000]', fontsize=16)

plt.ylabel('TESS Magnitude', fontsize=16)
plt.ylim(20.4, 19.4)
plt.xlim(3802.5, 3828.5)

plt.savefig('../figures/lightcurve_only.png', dpi=300, bbox_inches='tight')
