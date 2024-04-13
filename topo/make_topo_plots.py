"""
Plot the topography and contours of the slump.
1/9 arcsecond topo is used to create nice hillshade plots,
and NuuHillshadeOffshore.npy is saved for use in plotting GeoClaw results.
"""

from pylab import *
from clawpack.geoclaw import topotools, dtopotools
from clawpack.visclaw import plottools
import convert_tif
import os

fname = 'ncei19_n20x75_w156x25_2021v1.tif'
topo = convert_tif.read_tif(fname)

extent = [-156.182, -156.172, 20.62,20.63]
topoNuu = topo.crop(extent)

# Landslide slump:

slump = dtopotools.DTopography('../dtopo/slump_10m.dtt3',3)

offshore_extent = [-156.25, -156.12, 20.54, 20.64]
topoNuu2 = topo.crop(offshore_extent)

# Figure with contours:

figure(figsize=(12,12))
ax = axes()
#cb_kwargs = {'shrink': 0.6, 'extend':'both'}
#topoNuu2.plot(axes=ax, limits=[-100,15],cb_kwargs=cb_kwargs)
tN = topoNuu2
#contour(tN.X, tN.Y, tN.Z, linspace(0,400,9), colors='g')
#contour(tN.X, tN.Y, tN.Z, linspace(-1600,0,33), colors='b', linestyles='-')
contour(tN.X, tN.Y, tN.Z, arange(0,101,10), colors='g',linewidths=0.8)
contour(tN.X, tN.Y, tN.Z, arange(-2000,-1,100), colors='b', linestyles='-',linewidths=0.8)
title("Nu'u study area and offshore slump", fontsize=16);

#NuuBox = [-156.182, -156.172, 20.62,20.63]
NuuBox= [-156.1825,-156.175,20.625,20.629]  # to agree with fgout plots
print('NuuBox = ',NuuBox)
plottools.plotbox(NuuBox,{'linewidth':1.5, 'color':'r'})

contour(slump.X, slump.Y, slump.dZ[0,:,:], arange(-9,10,2), colors='k', linewidths=2)

ax.set_aspect(1/cos(20.6*pi/180))
xticks(fontsize=12)
yticks(fontsize=12)
xlabel('longitude', fontsize=14)
ylabel('latitude', fontsize=14)

if 0:
    fname = 'slump_offshore_nuu.jpg'
    savefig(fname, bbox_inches='tight')
    print('Created ',fname)
    
# Figure with hillshade:
# Based on https://matplotlib.org/stable/gallery/specialty_plots/topographic_hillshading.html

import matplotlib.pyplot as plt
from matplotlib.colors import LightSource

X,Y,z = tN.X, tN.Y, flipud(tN.Z)
dx = tN.delta[0] * 111e3 * cos(21*pi/180)
dy = tN.delta[1] * 111e3

# Shade from the northwest, with the sun 45 degrees from horizontal
ls = LightSource(azdeg=315, altdeg=45)
cmap = plt.cm.gist_earth

# compute hillshade with blended colors:
rgb = ls.shade(z, cmap=cmap, blend_mode='hsv',
                       vert_exag=15, dx=dx, dy=dy)
                       
figure(figsize=(12,12))
ax = axes()

imshow(rgb, extent=offshore_extent)

contour(tN.X, tN.Y, tN.Z, [0], colors='b', linestyles='-',linewidths=0.8)
contour(slump.X, slump.Y, slump.dZ[0,:,:], arange(-9,10,2), colors='k', linewidths=2)

title("Nu'u study area and offshore slump", fontsize=16);

#NuuBox = [-156.182, -156.172, 20.62,20.63]
plottools.plotbox(NuuBox,{'linewidth':2.5, 'color':'r'})

contour(slump.X, slump.Y, slump.dZ[0,:,:], arange(-9,10,2), colors='yellow', linewidths=2)

ax.set_aspect(1/cos(20.6*pi/180))
xticks(fontsize=12)
yticks(fontsize=12)
xlabel('longitude', fontsize=14)
ylabel('latitude', fontsize=14)

if 1:
    fname = 'slump_offshore_nuu_hillshade.jpg'
    savefig(fname, bbox_inches='tight')
    print('Created ',fname)

# Zoom on study area as background:                       
figure(figsize=(12,12))
ax = axes()

imshow(rgb, extent=offshore_extent)
contour(tN.X, tN.Y, tN.Z, range(4), colors='b', linestyles='-',linewidths=0.8)
axis(NuuBox)
title("Nu'u Study Area, contours at 0, 1, 2, 3 meters", fontsize=15);
ax.set_aspect(1/cos(20.6*pi/180))
ticklabel_format(useOffset=False)
xticks(fontsize=12)
yticks(fontsize=12)
xlabel('longitude', fontsize=14)
ylabel('latitude', fontsize=14);

if 1:
    fname = 'nuu_topo_hillshade.jpg'
    savefig(fname, bbox_inches='tight')
    print('Created ',fname)
    
# Save rgb for use as a background image:
fname = 'NuuHillshadeOffshore.npy'
save(fname, rgb)
print('Created ',fname)

fname_extent = fname.split('.')[0] + '_extent.txt'
savetxt(fname_extent, offshore_extent)
print('Created ',fname_extent)

# Test reloading:
extent = loadtxt(fname_extent)
rgb2 = load(fname)
figure(figsize=(8,8))
imshow(rgb2, extent=extent)
title('After reloading %s' % fname)
        
