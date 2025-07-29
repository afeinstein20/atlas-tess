"""
Figure originally created by Henry Hsieh
"""

import numpy as np
from astropy import units
from astropy.time import Time
from sbpy.photometry import HG
import matplotlib.pyplot as plt
from astropy.table import Table
from astroquery.jplhorizons import Horizons
from matplotlib.ticker import NullFormatter, MaxNLocator, MultipleLocator, FormatStrFormatter

rcParams = Table.read('../rcParams.txt', format='csv')
for key, val in zip(rcParams['name'], rcParams['value']):
    plt.rcParams[key] = val

data_filepath = '../data/ztf_data.txt'
data2i_filepath = '../data/ztf2i_data.txt'

majorFormatter0 = FormatStrFormatter('%d')
majorFormatter1 = FormatStrFormatter('%.1f')
majorFormatter2 = FormatStrFormatter('%.2f')
nullfmt         = NullFormatter()         # no labels

# 3I data
i3_data = np.loadtxt(data_filepath)
jd      = i3_data[:,0]
jdoff   = i3_data[:,0] - 2460800.0
absmag  = i3_data[:,1] * 1.0
absmag_err = i3_data[:,2] * 1.0

# load 3I ephemeris (for JD <-> r conversion)
obj3i = Horizons(id='3I', location='I41',
                 epochs={'start' : Time(jd[0]-10,format='jd').iso,
                         'stop'  : Time(jd[-1]+10,format='jd').iso,
                         'step'  : '1d'})
eph3i = obj3i.ephemerides()

# 2I data
data2i = Table.read(data2i_filepath, format='ascii')

# load 2I ephemeris (for absmag calculation + JD <-> r conversion)
obj2i = Horizons(id='2I', location='I41', epochs=data2i['t_jd'])
eph2i = obj2i.ephemerides()

jdequiv2i = np.interp(eph2i['r'], eph3i['r'][::-1], eph3i['datetime_jd'][::-1])
absmag2i = data2i['rmag'] - 5*np.log10(eph2i['r']) - 5*np.log10(eph2i['delta']) \
         - HG(H=0, G=0.15)(eph2i['alpha'].quantity) \
         + (4.81 - 4.65) # r' mag -> V mag

dt_factor = eph2i['r_rate'] / np.interp(eph2i['r'], eph3i['r'][::-1], eph3i['r_rate'][::-1])

plt.figure(1,figsize=(14,14))
rect_plot = [0.1,0.1,0.5,0.25]
axPlot    = plt.axes(rect_plot)

xmin, xmax = 2, 50
ymin, ymax = 10, 15
xticks_major = np.arange(xmin, xmax, 5)
xticks_minor = np.arange(xmin, xmax, 1)
yticks_major = np.arange(ymin, ymax, 0.2)
yticks_minor = np.arange(ymin, ymax, 0.05)

c1 = '#de4f0d'
c2 = '#555a5a'

# Plot 3I/ATLAS data (non-TESS)
axPlot.errorbar(jdoff[2:],absmag[2:],
                yerr=absmag_err[2:], ecolor='k', lw=1.5,
                fmt='o', color='k', mfc='k', ms=8, zorder=9, markeredgecolor='k', markeredgewidth=1, label='3I/ATLAS')


# Plot 3I/ATLAS data (TESS)
axPlot.errorbar(jdoff[:2],absmag[:2],
                yerr=absmag_err[:2], xerr=[5.24, 6.55], ecolor=c1,
                fmt='o', mfc='k', ms=10, zorder=9, markeredgecolor=c1, markeredgewidth=1,
                lw=2)

# Plot 2I/Borisov data
axPlot.errorbar(jdequiv2i - 2460800, absmag2i,
                xerr=np.array([data2i['dt1'], data2i['dt2']])*dt_factor, yerr=data2i['rerr'],
                fmt='o', color=c2, mfc='w', ms=8, zorder=9, label='2I/Borisov', markeredgecolor=c2)

for hline in [12.38, 12.28, 12.48]:
    if hline != 12.38:
        style='dotted'
    else:
        style='dashed'

    axPlot.hlines(hline,0,50,colors='#000000',linestyles=style, alpha=0.5)

axPlot.tick_params(axis='both',which='major',labelsize=18,length=7)
axPlot.tick_params(axis='both',which='minor',labelsize=0,length=3)
axPlot.tick_params(axis='x',which='major',labeltop=False) #length=0
axPlot.tick_params(axis='x',which='minor',labeltop=False)

axPlot.xaxis.set_ticks_position('both')
axPlot.yaxis.set_ticks_position('both')
axPlot.set_xticks(xticks_major)
axPlot.set_xticks(xticks_minor,minor=True)
axPlot.set_yticks(yticks_major)
axPlot.set_yticks(yticks_minor,minor=True)
axPlot.tick_params(which='both',direction='in')
axPlot.xaxis.set_major_formatter(majorFormatter0)
axPlot.yaxis.set_major_formatter(majorFormatter1)
axPlot.set_xlim([2,47])
axPlot.set_ylim([11.9,13.0])
axPlot.invert_yaxis()
axPlot.set_xlabel('Time [JD - 2460800]',fontsize=18)
axPlot.set_ylabel(r'$H_V$ [mag]',fontsize=18)
axPlot.legend(fontsize=18)

xticks2_major = np.arange(4,8,0.2)
xticks2_minor = np.arange(4,8,0.05)

jdoff_to_r = lambda jdoff: np.interp(jdoff + 2460800, eph3i['datetime_jd'], eph3i['r'])
r_to_jdoff = lambda r: np.interp(r, eph3i['r'][::-1], eph3i['datetime_jd'][::-1]) - 2460800
ax2 = axPlot.secondary_xaxis('top', functions=(jdoff_to_r, r_to_jdoff))
ax2.set_xticks(xticks2_major)
ax2.set_xticks(xticks2_minor,minor=True)
ax2.set_xlabel('$r$ [au]', fontsize=18)
#ax2.tick_params(axis='both',which='major',labelsize=18)

axPlot.tick_params(axis='both',which='major',direction='in',labelsize=18,top=False,bottom=True,right=True,left=True)
axPlot.tick_params(axis='both',which='minor',direction='in',labelsize=0,top=False,bottom=True,right=True,left=True)
ax2.tick_params(axis='both',which='major',direction='in',labelsize=18,top=True,bottom=False,right=True,left=True,length=7)
ax2.tick_params(axis='both',which='minor',direction='in',labelsize=0,top=True,bottom=False,right=True,left=True,length=3)

#Show the plot
plt.draw()

# Save to a File
plt.savefig('../figures/secular_lc_atlas.pdf',format = 'pdf',dpi=400,bbox_inches='tight')
