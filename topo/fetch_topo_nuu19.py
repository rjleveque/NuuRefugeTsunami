"""
Create topofile nuu19s.asc, a plot of this topo.
Optionally, also create a kmz file that can be viewed in Google Earth.
"""

from pylab import *
import os,sys
from clawpack.visclaw import colormaps, plottools
from clawpack.geoclaw import topotools, kmltools
from clawpack.clawutil.data import get_remote_file
import glob
import zipfile
import convert_tif

topo_name = 'nuu19s'

# Read and crop DEM:
fname = 'ncei19_n20x75_w156x25_2021v1.tif'
if os.path.isfile(fname):
    print('Found %s, no need to download' % fname)
else:
    url = 'https://chs.coast.noaa.gov/htdata/raster2/elevation/' + \
          'NCEI_ninth_Topobathy_Hawaii_9428/tiles/' + fname
    get_remote_file(url, output_dir='.', file_name=fname, verbose=True)
    
topo = convert_tif.read_tif(fname)

#extent = [-156.186, -156.172, 20.618,20.63]
extent = [-156.22, -156.16, 20.61,20.63]
topo = topo.crop(extent)

# Create topofile:
fname = topo_name + '.asc'
topo.write(fname, topo_type=3, Z_format=' %.3f', header_style='asc')
print('Created ',fname)

# Plot topo:
fname = topo_name + '.png'
fig,ax = subplots(figsize=(9,5))
topo.plot(axes=ax,limits=(-50,20),cb_kwargs={'extend':'both', 'shrink':0.5})
title('%s' % topo_name)
savefig(fname)
print('Created ',fname)


# If desirec, also create a kmz file that can be opened in Google Earth:

if 0:
    # Create kmz file:
    X,Y,Z = topo.X, topo.Y, topo.Z

    zmin = -15.
    zmax = 15.

    cmap_land = colormaps.make_colormap({ 0.0:[0.1,0.4,0.0],
                                         0.25:[0.0,1.0,0.0],
                                          0.5:[0.8,1.0,0.5],
                                          1.0:[0.8,0.5,0.2]})

    cmap_sea = colormaps.make_colormap({ 0.0:[0,0,1], 1.:[.8,.8,1]})

    cmap_topo, norm_topo = colormaps.add_colormaps((cmap_land, cmap_sea),
                                         data_limits=(zmin,zmax),
                                         data_break=0.)

    cmap_sea_dry = colormaps.make_colormap({ 0.0:[1.0,0.7,0.7], 1.:[1.0,0.7,0.7]})
    cmap_dry, norm_dry = colormaps.add_colormaps((cmap_land, cmap_sea_dry),
                                         data_limits=(zmin,zmax),
                                         data_break=0.) 
                                         
                                         

    kml_dir = 'kmlfiles_%s' % topo_name
    os.system('mkdir -p %s' % kml_dir)
    print('Will put png and kml files in %s' % kml_dir)

    Z_land = ma.masked_where(Z<0., Z)
    png_filename = '%s/%s_land.png' % (kml_dir, topo_name)
    fig,ax,png_extent,kml_dpi = kmltools.pcolorcells_for_kml(X,Y,Z_land,
                     png_filename=png_filename, dpc=2, 
                      cmap=cmap_topo, norm=norm_topo) 

    Z_water = ma.masked_where(Z>=0, Z)
    png_filename = '%s/%s_water.png' % (kml_dir, topo_name)
    fig,ax,png_extent,kml_dpi = kmltools.pcolorcells_for_kml(X,Y,Z_water,
                 png_filename=png_filename, dpc=2, cmap=cmap_topo, norm=norm_topo)

    kmltools.kml_build_colorbar('%s/colorbar.png' % kml_dir, cmap_topo, 
               norm=norm_topo, label='meters', title='topo', extend='both')

    png_files=['%s_water.png' % topo_name, 
               '%s_land.png' % topo_name]
    png_names=['%s_water' % topo_name,
               '%s_land' % topo_name]
    cb_files = ['colorbar.png']
    cb_names = ['colorbar_topo']

    fname = os.path.join(kml_dir, topo_name+'.kml')
    kmltools.png2kml(png_extent, png_files=png_files, png_names=png_names, 
                     name=topo_name, fname=fname,
                     radio_style=False,
                     cb_files=cb_files, cb_names=cb_names)

    savedir = os.getcwd()
    os.chdir(kml_dir)
    files = glob.glob('*.kml') + glob.glob('*.png')
    print('kmz file will include:')
    for file in files:
        print('    %s' % os.path.split(file)[-1])

    fname_kmz = '%s.kmz' % topo_name
    with zipfile.ZipFile(fname_kmz, 'w') as zip:
        for file in files:
            zip.write(file) 
        print('Created %s' % os.path.abspath(fname_kmz))
    os.chdir(savedir)
              
