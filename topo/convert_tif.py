from pylab import *

def read_tif(path):
    """
    Read a geotiff topography DEM and return a topotools.Topography object.
    """
    
    try:
        import gdal
    except:
        from osgeo import gdal
    gdal.UseExceptions()

    from clawpack.geoclaw import topotools
    
    ds = gdal.Open(path)

    srcband = ds.GetRasterBand(1)
    z = srcband.ReadAsArray()

    transform = ds.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    dx = transform[1]
    dy = transform[5]
    
    dy = -dy
    Z = flipud(z)
    
    x = linspace(xOrigin, xOrigin + (z.shape[0]-1)*dx, z.shape[0])
    y = linspace(yOrigin - (z.shape[1]-1)*dy, yOrigin, z.shape[1])
    print('DEM has shape ', z.shape)
    print('x range:  %.6f, %.6f' % (x[0], x[-1]))
    print('y range:  %.6f, %.6f' % (y[0], y[-1]))
    
    topo = topotools.Topography()
    topo.set_xyZ(x,y,Z)
    return topo
    
def convert_tif_to_asc(fname_tif, fname_asc=None, extent=None, coarsen=1):
    """
    Read a tif and write a .asc raster file for use in GeoClaw 
    with topo_type==3.  Rounds off DEM z values to millimeter precision to
    minimize file size.
    """
    
    import os

    print('Reading ',fname_tif)
    topo = read_tif(fname_tif)
    
    if (extent is not None) or (coarsen > 1):
        topo = topo.crop(extent, coarsen=coarsen)
        
    if fname_asc is None:
        fname_asc = os.path.splitext(fname_tif)[0]
        if extent is not None:
            fname_asc = fname_asc + '_cropped'
        if coarsen > 1:
            fname_asc = fname_asc + '_coarsen%s' % coarsen
        fname_asc = fname_asc + '.asc'
        
    # write .asc file with millimeter precision:
    topo.write(fname_asc, topo_type=3, header_style='asc', Z_format='%.3f')
    print('Created ',fname_asc)
    
if __name__ == '__main__':
    # if executed at the Unix command line....
    import sys
    args = sys.argv[1:]   # any command line arguments
    fname_tif = args[0]   # always assume fname_tif is given
    if len(args) > 1:
        coarsen = int(args[1])
        print('Will coarsen by %i' % coarsen)
    else:
        coarsen = 1
    convert_tif_to_asc(fname_tif, coarsen=coarsen)

        
