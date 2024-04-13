"""
Simple script to download topofile
"""

import os
from clawpack.geoclaw import topotools

plot_topo = True

url = 'https://www.ngdc.noaa.gov/thredds/dodsC/regional/hawaii_6_mllw_2005.nc'

extent = [-157, -156, 20, 21]
coarsen = 1
topo = topotools.read_netcdf(url, extent=extent, 
                             coarsen=coarsen, verbose=True)

name = 'hawaii6s_nuu'
print('name = ',name)
fname = name + '.asc'
topo.write(fname, topo_type=3, header_style='asc', 
                     grid_registration='llcorner', Z_format='%.0f')



if plot_topo:
    # plot the topo and save as a png file...
    import matplotlib.pyplot as plt
    topo.plot()
    plt.title('Topo file %s' % name)
    fname = name + '.png'
    plt.savefig(fname)
    print('Created %s' % fname)
