"""
Use PyVista to plot the tsunami in the Nu'u pond region from fgout grid #5.
"""

from pylab import *
import os
import pyvista as pv
from clawpack.geoclaw import fgout_tools
import debris_tools
from clawpack.visclaw import animation_tools
from nuu_debris import compute_debris_paths

global etamesh, debris_spheres

fgno = 5  # which fgout grid

mstr = '10'
outdir = '_output_%sm' % mstr
format = 'binary32'  # format of fgout grid output
qmap = 'geoclaw-bouss'  # defines mapping for fgout variables


# where to save a png file for each frame, for constructing an animation:
framedir = '_frames_%sm' % mstr
os.system('mkdir -p %s' % framedir)
os.system('rm %s/*' % framedir)  # remove frames from previous version

tfinal = 15*60.  # create animation up to this time (or end of computation)

# compute debris paths based on fgout2:
debris_paths, dbnos_water, dbnos_land = compute_debris_paths(mstr,tfinal,
                                                             dxd=0.0002)

debris_spheres = len(debris_paths)*[' ']  # for plot actors

# Instantiate object for reading fgout frames:
fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, format, qmap) 
#fgout_grid.read_fgout_grids_data_pre511(fgno)
try:
    fgout_grid.read_fgout_grids_data(fgno)
except:
    fgout_grid.read_fgout_grids_data_pre511(fgno)
    fgout_grid.times = linspace(fgout_grid.tstart, fgout_grid.tend,
                                fgout_grid.nout)
fgout = fgout_grid.read_frame(1)

# convert x,y to meters, roughly:

x = (fgout.x - fgout.x[0]) * 111e3 * cos(fgout.y.mean()*pi/180)
y = -(fgout.y - fgout.y[0]) * 111e3
z = array([0.])
X,Y,Z = meshgrid(x, y, z, indexing='ij')
topoxyz = pv.StructuredGrid(X,Y,Z)

B = fgout.B
B = flipud(B)
Bmax = 25.
B = minimum(B, Bmax)
topoxyz.point_data['B'] = B.flatten(order='F')

warpfactor = 8  # amplification of elevations

topowarp = topoxyz.warp_by_scalar('B', factor=warpfactor)

texture = pv.read_texture('./GE_Nuu.jpg')
GE_extent = fgout.extent_edges
x2,x1 = (asarray(GE_extent[:2]) - fgout.x[0]) * 111e3 * cos(fgout.y.mean()*pi/180)
y2,y1 = (asarray(GE_extent[2:]) - fgout.y[0]) * 111e3

#x1,x2,y1,y2 = GE_extent

o = (x1, y1, 0.)
u = (x2, y1, 0.)
v = (x1, y2, 0.)

mapped_surf = topowarp.texture_map_to_plane(o, u, v)

if 0:
    # plot topo flat:
    p = pv.Plotter()
    p.add_mesh(topoxyz,cmap='gist_earth',clim=(-20,10))
    p.view_xy()
    p.show(window_size=(1500,1500))

if 0:
    # plot topo alone:
    p = pv.Plotter()
    p.add_mesh(topowarp,cmap='gist_earth',clim=(-20,20))
    p.camera_position = 'xz'
    p.camera.azimuth = 180
    p.camera.elevation = 20
    p.show(window_size=(1500,1500))

# replace land surface by nan in eta so it only shows water:
eta = where(fgout.h>0.1, fgout.eta, nan)
eta = flipud(eta)

make_animation = True

# select frames for fgout grid 5 every 10 seconds up to tfinal:
fgframes = [n+1 for n in range(len(fgout_grid.times)) \
               if ((fgout_grid.times[n] <= tfinal) and \
                  (mod(fgout_grid.times[n], 10) == 0))]

p = pv.Plotter(off_screen=make_animation, lighting='three lights')

#p.add_mesh(topowarp,cmap='gist_earth',clim=(-20,20))
p.add_mesh(mapped_surf,texture=texture)

topoxyz.point_data['eta'] = eta.flatten(order='F')
etawarp = topoxyz.warp_by_scalar('eta', factor=warpfactor)
etamesh = p.add_mesh(etawarp,color='c')

if 0:
    p.camera_position = 'xz'
    p.camera.azimuth = 120
    p.camera.elevation = 40
    
camera_position1 =  [(2351.473476407906, 223.43232607508872, 1290.8736572490081),
 (387.9651184082031, -442.28570556640625, -49.22932434082031),
 (-0.4995685365943228, -0.21553276761215126, 0.8390333147917629)]
 
p.camera_position = camera_position1

p.window_size=(1500,1500)


def set_frameno(fgframeno):
    global etamesh, debris_spheres
    #fgframeno = 6*int(floor(tmin)) + 1
    fgframeno = int(round(fgframeno))
    fgout = fgout_grid.read_frame(fgframeno)
    thour, remainder = divmod(fgout.t, 3600)
    tmin, tsec = divmod(remainder, 60)
    tstr = '%s:%s:%s' % (str(int(thour)).zfill(2),
                      str(int(tmin)).zfill(2),
                      str(int(tsec)).zfill(2))
    print('Frame %i, t = %s' % (fgframeno, tstr))

    eta = where(fgout.h>0.1, fgout.eta, nan)
    eta = flipud(eta)
    topoxyz.point_data['eta'] = eta.flatten(order='F')
    etawarp = topoxyz.warp_by_scalar('eta', factor=warpfactor)
    p.remove_actor(etamesh)
    etamesh = p.add_mesh(etawarp,color='c')
    
    p.add_title('Time %s\nSlump amplitude %s m' % (tstr,mstr))

    # add debris particles at this time:
    xd_water,yd_water = debris_tools.get_debris_xy(fgout.t, debris_paths,
                                                   dbnos_water)
    xd_land,yd_land = debris_tools.get_debris_xy(fgout.t, debris_paths, 
                                                   dbnos_land)
    # make function that returns eta(x,y) at this time:
    eta_fcn = fgout_tools.make_fgout_fcn_xy(fgout, 'eta')
    B_fcn = fgout_tools.make_fgout_fcn_xy(fgout, 'B')
    
    
    for k in range(len(xd_water)):
        xd = -(xd_water[k]- fgout.x[-1]) * 111e3 * cos(fgout.y.mean()*pi/180)
        yd = -(yd_water[k]- fgout.y[0]) * 111e3
        sphere_radius = 5    
        zd = warpfactor * max(eta_fcn(xd_water[k],yd_water[k]),
                 B_fcn(xd_water[k],yd_water[k])+sphere_radius/warpfactor)
        debris_sphere = pv.Sphere(radius=sphere_radius, center=(xd,yd,zd))
        p.remove_actor(debris_spheres[k])
        debris_spheres[k] = p.add_mesh(debris_sphere, color='r')
        #print('Added sphere at (%.3f, %.3f, %.3f)' % (xd,yd,zd))
        
    if not make_animation:
        # print camera position so that this can be copied and pasted
        # into this script after adjusting (and then sliding frameno)
        print('p.camera_position = ', p.camera_position)                   
    
        
    # to capture each frame after moving slider, uncomment this:
    #p.screenshot('PyVista_Nuu_Frame%s.png' % str(fgframeno).zfill(4))
    
#p.add_slider_widget(set_frameno, [0,20], value=0, title='Time (minutes)')
#p.add_slider_widget(set_frameno, [0,4], value=0, title='Time (minutes)')




if not make_animation:
    frames_range = [fgframes[0], fgframes[-1]]
    p.add_slider_widget(set_frameno, frames_range, value=1, title='Frameno',
                    pointa=(0.4,0.85), pointb=(0.9,0.85), color='blue',
                    slider_width=0.02, tube_width=0.005)

    p.show()
    
else:

    #for fgframeno in range(1,121,1):
    #for fgframeno in range(1,74,1):
    for fgframeno in fgframes:
        set_frameno(fgframeno)
        fname_png = '%s/PyVistaFrame%s.png' % (framedir, str(fgframeno).zfill(4))
        p.screenshot(fname_png)
        print('Created ',fname_png)
        
    p.close()

    anim = animation_tools.make_anim(framedir, fname_pattern='PyVistaFrame*.png')
    fname_mp4 = 'Nuu_debris_animation%sm_10sec.mp4' % mstr
    animation_tools.make_mp4(anim, fname_mp4)
    print('Created ',fname_mp4)
