#!/usr/bin/env python
# coding: utf-8

# # Make Slump Deformation
# 
# Supplementary material for the paper *Hawaiian legends of coastal devastation and paleotsunami reconstruction, Nuʻu, Kaupō, Maui, Hawaiʻi.,* by
# Scott Fisher, James Goff, Andrew Cundy, David Sear, James Terry Randall J LeVeque, and Loyce M Adams.  (Submitted for publication, April 2024).
# 
# This notebook explains how the landslide slump deformation used in the GeoClaw tsunami model was defined. It uses a number of Python tools from [GeoClaw](http://www.geoclaw.org) that are not described in detail here.  For information on these tools and on installing Clawpack, the open source package that includes GeoClaw, see www.clawpack.org.
# 
# This Jupyter notebook and and the GeoClaw code for running the tsunami simulations shown in the paper can be found at https://github.com/rjleveque/NuuRefugeTsunami.
# 

# ### Set up the notebook and import various packages:

# In[1]:


#get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


from pylab import *
from clawpack.geoclaw import util,topotools,dtopotools
from clawpack.visclaw import colormaps, animation_tools
from scipy.interpolate import RegularGridInterpolator, griddata
from IPython.display import HTML
import warnings
warnings.filterwarnings("ignore")


# The GeoClaw software requires one or more topography DEMs (topo files) to specify the underlying topography.  In addition, a "dtopo file" can be provided that gives changes to the topography at a series of times, as described in [the documentation](https://www.clawpack.org/topo.html#topography-displacement-files).
# 
# The approach used here is to specify the slump as a dipolar source that has a positive hump on the up-slope side and a negative hole on the down-slope side, which relaxes to zero over some time interval so that in effect the mass is transfered from the hump to the hole.  This dtopo function is added to the underlying topography, so that at the initial time t0 the topography is deformed and at the end of the slump duration it agrees with the topography specified by the present-day DEM.
# 
# This same approach has been used in several other papers.  See for example the citations from our paper to Okal and Hebert, 2007, Tappin et al. 2007, Watts et al. 2005, Fryer et al. 2004, Synolakis et al. 2002. We use a similar dipolar slump shape.
# 
# We first create this source with the desired dimensions specified in meters, and with the slump aligned with the x-axis.  Then we can convert it to longitude-latitude coordinates and place it where desired on the topography.  The resulting deformation at a sequence of times is then evaluated on a specified uniform rectantular grid in [UTM coordinates](https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system), and saved as a "dtopo file" for GeoClaw.

# Scroll to the bottom of this notebook to see the slump contours plots along with the topography.

# ### Specify the variation in the x direction along the centerline:
# 
# This shows the dipolar source along the x-axis. This will then be tapered in the y-direction to give a positive hump and negative hole.

# In[3]:


L = 3e3  # length in meters
A = 10.  # amplitude (approximate) in meters

slump_profile = lambda x: 2*A*where(abs(x)<L, sin(x*pi/L)/cosh(2.5*x/L)**2, 0.)

figure(figsize=(8,4))
xm = linspace(-L, L, 500)
plot(xm, slump_profile(xm),'b')
grid(True)
xlabel('x (meters)')
ylabel('vertical deformation (m)')
title('Slump profile')


# ### Taper in y
# The function shown above is now multiplied by a function $sech^2(3y/W)$ that tapers this down in the y-direction. 

# In[4]:


W = L    # width
ytaper = lambda y:  (1./cosh(3*(y)/W))**2
ym = linspace(-W,W,501)

figure(figsize=(8,4))
plot(ym, ytaper(ym), 'b')
grid(True)
xlabel('y (meters)')
ylabel('nondimensional')
title('Taper in y');


# # Plot how the initial slump deformation looks
# 
# The initial dtopo is the product `slump_profile(x) * ytaper(y)`, shown below on a contour plot after evaluating it on a 2D grid of points:
# 
# In this plot, `(x,y)` are in meters and the slump is aligned with the x-axis:

# In[5]:


Xm,Ym = meshgrid(xm,ym)
Zm0 = slump_profile(Xm) * ytaper(Ym)
pcolormesh(Xm,Ym,Zm0,cmap=colormaps.blue_white_red)
clim(-15,15)
colorbar();
contour(Xm,Ym,Zm0,linspace(-20,20,20),colors='k');
axis('equal')
xlabel('x (meters)')
ylabel('y (meters)')
title('Initial slump deformation in x-y plane');


# ## Slump placement on UTM grid
# 

# We chose a slump location based on the local offshore topography near Nu'u, to be at a spot where the topography is somewhat steep.  The point (0,0) in the plots shown above will be mapped to the location 
# `(x0,y0)` specified below. 

# In[6]:


x0 = -156.2  # longitude
y0 = 20.59   # latitude


# We assume the slump dimensions are small enough that linearizing the transformation from meters to degrees about the center point `(x0,y0)` is sufficient.  Define factors to map `dx` in degrees to `dx` in meters and similarly in y:

# In[7]:


Rearth = 6367500. # radius of earth
LON2M = Rearth*pi/180 * cos(y0*pi/180)
LAT2M = Rearth*pi/180


# Compute a bearing vector `(bx,by)` that points in the desired direction of the slump axis.  Here we take this to be the direction from `(x0,y0)`, the slump center to `(x_nuu, y_nuu)`, a point near the Nu'u Refuge on Maui. This vector is converted to meters and normalized to have length 1.  

# In[8]:


x_nuu = -156.187
y_nuu = 20.6275

theta = util.bearing(x0,y0,x_nuu,y_nuu,bearing_units='degrees')
theta = theta - 90 # rotation angle (not used below)
print('Rotation angle: %.1f degrees clockwise' % theta)

bx = (x_nuu - x0) * LON2M
by = (y_nuu - y0) * LAT2M
bb = sqrt(bx**2 + by**2)
bx = bx/bb
by = by/bb
print('Bearing vector b in meters: (%.3f, %.3f)' % (bx,by))
angle = -arctan(by/bx) * 180/pi
print('Angle computed from b vector: %.1f' % angle)


# ### Create uniform UTM grid on which to evaluate the slump, as needed in GeoClaw:

# In[9]:


# Choose a rectangle large enough to contain the slump:
slump_radius = sqrt(L**2 + W**2)
x1 = x0 - slump_radius / LON2M
x2 = x0 + slump_radius / LON2M
y1 = y0 - slump_radius / LAT2M
y2 = y0 + slump_radius / LAT2M
dx = dy = 1/3600.

x = arange(x1,x2+dx/2,dx)
y = arange(y1,y2+dy/2,dy)
X,Y = meshgrid(x,y)

print('Will evaluate slump on grid of dimensions ', X.shape)
print('For %.4f <= x <= %.4f and %.4f <= x <= %.4f' % (x1,x2,y1,y2))
print('With grid resolution %.3f arcseconds' % (dx*3600))


# Now convert the points `(X,Y)` on the UTM grid to distances (meters) from `(x0,y0)` and project onto the bearing vector and the perpendicular vector in order to compute `(Xm,Ym)`, distances along these two directions:

# In[10]:


dX0 = (X-x0) * LON2M
dY0 = (Y-y0) * LAT2M
Xm = bx*dX0 + by*dY0
Ym = by*dX0 - bx*dY0


# The initial slump shape can now be evaluated based on these distances:

# In[11]:


dZ0 = slump_profile(Xm) * ytaper(Ym)


# Plot the initial slump shape

# In[12]:


fig,ax = subplots()
pcolormesh(X,Y,dZ0,cmap=colormaps.blue_white_red)
clim(-15,15)
colorbar();
contour(X,Y,dZ0,linspace(-20,20,20),colors='k');
ax.set_aspect(1/cos(y0*pi/180))
xticks(rotation=20)
xlabel('x (longitude)')
ylabel('y (latitude)')
title('Initial slump deformation');


# ## Create dtopo file for progressive slump
# 
# First make dZ at a set of times going from maximum displacment at time 0 to zero displacement at final time (so that final topo agrees with topofile).  As in other models, use decay in time given by $(cos(\pi t/T) + 1)/2$ for $0 \leq t \leq T$, where $T$ is the final time.

# In[13]:


Tfinal = 30 # seconds
ntimes = 16
times = linspace(0,Tfinal,ntimes)
dZshape = (ntimes,X.shape[0],X.shape[1])
dZ = empty(dZshape)
ampl = 0.5*(cos(pi*times/Tfinal) + 1)
for k in range(ntimes):
    dZ[k,:,:] = ampl[k] * dZ0

figure(figsize=(5,3))
plot(times,ampl)
xlim(0,Tfinal)
title('Decay of slump in time')
xlabel('t (seconds)')
grid(True);


# In[14]:


# create dtopo object:
dtopo_slump = dtopotools.DTopography()
dtopo_slump.dZ = dZ
dtopo_slump.X = X
dtopo_slump.Y = Y
dtopo_slump.times = times
dtopo_slump.delta = (dx,dy)


# ### Plot displacment on topo:
# 
# Download topography from the [Hawaii 6 arcsecond coastal DEM](https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.ngdc.mgg.dem:460/html) and plot the slump contours on the topography:

# In[15]:


url = 'https://www.ngdc.noaa.gov/thredds/dodsC/regional/hawaii_6_mllw_2005.nc'
extent = [-156.28,-156.12,20.55,20.7]
topo6s = topotools.read_netcdf(url, extent=extent, coarsen=1, verbose=True)


# In[16]:


fig,ax = subplots(figsize=(8,7))
topo6s.plot(axes=ax, contour_levels=arange(-1500,200,100))
dtopo_slump.plot_dZ_contours(5,dZ_interval = 2,axes=ax)
plot([-156.187],[20.6275],'ro',markersize=5)
title("Initial slump deformation and Nu'u (red dot)\n" \
      + "with 2m contours (100m topo contours)");
#savefig('slump_topo.png')


# ### Write dtopo file for this slump and scaled versions:

# In[17]:


dZ10m = dtopo_slump.dZ.copy()  # make a copy of the 10m amplitude dZ array
for ampl in [5,15,10]:
    dtopo_slump.dZ = (ampl/10.) * dZ10m
    fname_dtopo = 'slump_%sm.dtt3' % str(int(ampl)).zfill(2)
    dtopo_slump.write(fname_dtopo, dtopo_type=3)
    print('Created ',fname_dtopo)


# ### Make animation of slump behavior:

# In[18]:


figs = []
for k,tplot in enumerate(dtopo_slump.times[:-1]):
    fig,ax = subplots()
    dZ_interval = 0.1*A
    dtopo_slump.plot_dZ_colors(t=tplot,axes=ax,cmax_dZ=A, dZ_interval=dZ_interval)
    ticklabel_format(style='plain',useOffset=False)
    xticks(rotation=20);
    ax.set_aspect(1/cos(y0*pi/180))
    title('Slump deformation at time t = %.1f seconds\n Contour interval %.1f m' % (tplot,dZ_interval))
    figs.append(fig)
    close(fig)

images = animation_tools.make_images(figs)
anim = animation_tools.animate_images(images, figsize=(6,5))
animation_tools.make_mp4(anim, file_name='slump_10m.mp4')


# The resulting animation is archived with supplementary information for the paper, as the mpeg file `slump_10m.mp4`.
# 
# You can also view the animation here if you are running the notebook live, or viewing an html rendered version:

# In[19]:


HTML(anim.to_jshtml())


# In[ ]:




