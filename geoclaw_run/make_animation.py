"""
Plot fgout frames and also motion of particles using velocities
interpolated from fgout.
"""

import sys
if 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('Agg')  # Use an image backend

from pylab import *
import os
from clawpack.visclaw import plottools, geoplot, gridtools
from clawpack.visclaw import animation_tools, colormaps
from matplotlib import animation, colors
from clawpack.geoclaw import fgout_tools
from datetime import timedelta 
from scipy.interpolate import interp1d
from nuu_debris import compute_debris_paths


# local module:
import debris_tools
#from clawpack.geoclaw import debris_tools


mstr = '05'  # slump amplitude (rerun code for 05, 10, 15)
outdir = '_output_%sm' % mstr
#outdir = '_output'

print('Looking for output in ',outdir)

output_format = 'binary32'
qmap = 'geoclaw-bouss'  # defines mapping for fgout variables

graphics_dir = os.path.abspath('../topo')
hillshade_image = load(graphics_dir + '/NuuHillshadeOffshore.npy')
hillshade_extent = loadtxt(graphics_dir + '/NuuHillshadeOffshore_extent.txt')

# Instantiate object for reading fgout frames:
fgout_grid1 = fgout_tools.FGoutGrid(2, outdir, output_format, qmap)
fgout_grid1.read_fgout_grids_data(2)
fgout_grid2 = fgout_tools.FGoutGrid(5, outdir, output_format, qmap)
fgout_grid2.read_fgout_grids_data(5)


# transect:
x1trans = -156.17812; y1trans = 20.6257
x2trans = -156.1766; y2trans = 20.6280
xtrans = linspace(x1trans, x2trans, 1000)
ytrans = linspace(y1trans, y2trans, 1000)

# fgout frames to include in animation:

tfinal = 20*60.  # create animation up to this time (or end of computation)

# use all frames for fgout grid 2 up to time tfinal (spaced every 10 sec)
fgframes1 = [n+1 for n in range(len(fgout_grid1.times)) \
               if fgout_grid1.times[n] <= tfinal]
print('+++ fgframes1 = ',fgframes1)
print('For fgout grid 2, using %i frames up to frame %i, t = %i sec' \
        % (len(fgframes1), fgframes1[-1], fgout_grid1.times[fgframes1[-1]-1]))
        

# select frames for fgout grid 5 every 10 seconds up to tfinal:
fgframes2 = [n+1 for n in range(len(fgout_grid2.times)) \
               if ((fgout_grid2.times[n] <= tfinal) and \
                  (mod(fgout_grid2.times[n], 10) == 0))]

#print('+++ fgframes2 = ',fgframes2)
print('For fgout grid 5, using %i frames up to frame %i, t = %i sec' \
        % (len(fgframes2), fgframes2[-1], fgout_grid2.times[fgframes2[-1]-1]))
        
fgframe1 = fgframes1[0] # start with first frame
fgframe2 = fgframes2[0] # for fgout grid 2 

fgout1 = fgout_grid1.read_frame(fgframe1) # fgout grid offshore region
fgout2 = fgout_grid2.read_frame(fgframe2) # finer fgout grid over Nu'u

# compute debris paths based on fgout2:
debris_paths, dbnos_water, dbnos_land = compute_debris_paths(mstr, tfinal)



# ----------
# Plotting:

# Plot one frame of fgout data and define the Artists that will need to
# be updated in subsequent frames:

fig = figure(1, figsize=(13,8))
clf()

# =========
# Ocean:
# =========

#plot_extent = fgout.extent_edges
#plot_extent = [-160.7, -154.54, 18.5, 22.6]
plot_extent = [-156.4, -156., 20.4, 20.7]

ax = axes([.1,.1,.4,.8])
B_plot1 = ax.imshow(flipud(fgout1.B.T), extent=fgout1.extent_edges,
       cmap=geoplot.land1_colormap)

B_plot1.set_clim(0,1500)

eta_water = where(fgout1.h > 0, fgout1.eta, nan)
eta_plot1 = ax.imshow(flipud(eta_water.T), extent=fgout1.extent_edges,
       cmap=geoplot.tsunami_colormap)

eta_plot1.set_clim(-1,1)
axis(plot_extent)
title_text = ax.set_title('Slump amplitude %s meters\nSurface at time %s' \
        % (mstr,timedelta(seconds=fgout1.t)))


if 1:
    cb = colorbar(eta_plot1, extend='both', shrink=0.5, 
                  #orientation='vertical',anchor=(0,0))
                  orientation='horizontal', anchor=(0.2,1))
    cb.set_label('meters')
    
# =========
# Nu'u region
# =========


#onecolor = [0.3,0.3,1]
onecolor = [0.7,0.7,1]
cmap_onecolor = colormaps.make_colormap({0.:onecolor, 1.:onecolor})

# inset axes....
axins = ax.inset_axes([1.1,0.40,1.0,.7])
#axins.imshow(GE_image, extent=GE_extent)
axins.imshow(hillshade_image, extent=hillshade_extent)

#B = ma.masked_where(fgout2.h < 0.01, fgout2.B)
B = nan*fgout2.B  # don't show topo since using GE_image
B_plot2 = axins.imshow(flipud(B.T), extent=fgout2.extent_edges,
       cmap=geoplot.googleearth_transparent)
B_plot2.set_clim(0,8)

dry_depth = 0.02  # gdepth
eta_water2 = where(fgout2.h > 0, fgout2.eta, nan)
h_water2 = where(fgout2.h > dry_depth, fgout2.h, nan)
zeta_water2 = where(fgout2.B>0, h_water2, eta_water2)
eta_plot2 = axins.imshow(flipud(zeta_water2.T), extent=fgout2.extent_edges,
       cmap=cmap_onecolor)
eta_plot2.set_clim(-5,5)
axins.plot([x1trans,x2trans], [y1trans,y2trans], 'b', linewidth=0.9)
axins.text(x2trans-0.0001, y2trans-0.0004, 'Transect', color='b', fontsize=10)
#xticks(rotation=20)

# sub region of the original image
x1, x2, y1, y2 = [-156.1825,-156.175,20.625,20.629]
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xticks(arange(-156.182,-156.175,.002))
axins.ticklabel_format(useOffset=False)

axins.set_yticklabels([])
axins_yaxis2 = axins.secondary_yaxis('right')
axins_yaxis2.ticklabel_format(useOffset=False)
axins.set_title("Nu'u Refuge Pond")
axins.contour(fgout2.X,fgout2.Y,fgout2.B,[1],colors='g')

# show where this region is on the Ocean plot:
ax.indicate_inset_zoom(axins, edgecolor="black")

# -----------
# add debris:
            
xd_water,yd_water = debris_tools.get_debris_xy(fgout2.t, debris_paths, dbnos_water)
xd_land,yd_land = debris_tools.get_debris_xy(fgout2.t, debris_paths, dbnos_land)

if 0:
    points_water, = axins.plot(xd_water, yd_water, color='k',
                            linestyle='', marker='.',markersize=5)
    points_land, = axins.plot(xd_land, yd_land, color='k',
                            linestyle='', marker='s',markersize=2)

#color_rock = 'yellow'; color_coral = 'w'
color_rock = 'k'; color_coral = 'r'
points_water, = axins.plot(xd_water, yd_water,
                        color=color_coral,linestyle='',
                        marker='.',markersize=5)
points_land, = axins.plot(xd_land, yd_land,
                        color=color_rock,linestyle='',
                        marker='s',markersize=2)
                                                    

# =========
# transect:
# =========


def extract_transect(fgout_soln):

    eta1d = gridtools.grid_eval_2d(fgout_soln.X, fgout_soln.Y,
                                   fgout_soln.eta, xtrans, ytrans)
    B1d = gridtools.grid_eval_2d(fgout_soln.X, fgout_soln.Y,
                                 fgout_soln.B, xtrans, ytrans)
    return B1d, eta1d

                                 

axtrans = axes([.55,.1,.4,.3])
axtrans.set_title('Transect')

axtrans.set_xlim(x1trans,x2trans)
axtrans.set_ylim(-10,15)

Btrans1, etatrans1 = extract_transect(fgout1)
Btrans2, etatrans2 = extract_transect(fgout2)

Btrans, etatrans = Btrans2, etatrans2


# filled regions:
Bfill_plot = axtrans.fill_between(xtrans, Btrans-1e4, Btrans, 
                                  color=[.5,1,.5,1])
etafill_plot = axtrans.fill_between(xtrans, Btrans, etatrans, 
                                  color=[.5,.5,1,1])

# surface and topo plots:
etatrans_plot, = axtrans.plot(xtrans, etatrans, 'b')
Btrans_plot, = axtrans.plot(xtrans, Btrans, 'g')


axtrans.grid(True)
axtrans.ticklabel_format(useOffset=False)



# The artists that will be updated for subsequent frames:
update_artists = (B_plot1, eta_plot1, B_plot2, eta_plot2,
                      Bfill_plot, etafill_plot, Btrans_plot, etatrans_plot,
                      title_text, points_water, points_land)
            
figdummy,axdummy = subplots()

def update(k, *update_artists):
    """
    Update an exisiting plot with solution from fgout frame fgframeno.
    The artists in update_artists must have new data assigned.
    """

    #fgout1 = fgout_grid1.read_frame(fgframeno)
    #fgout2 = fgout_grid2.read_frame(2*fgframeno-1)

    fgout1 = fgout_grid1.read_frame(fgframes1[k])
    fgout2 = fgout_grid2.read_frame(fgframes2[k])

    print('Updating plot at time %s and %s' \
            % (timedelta(seconds=fgout1.t), timedelta(seconds=fgout2.t)))
    
    # unpack update_artists (must agree with definition above):
    B_plot1, eta_plot1, B_plot2, eta_plot2, \
          Bfill_plot, etafill_plot, Btrans_plot, etatrans_plot, \
          title_text, points_water, points_land = update_artists
        
    # reset title to current time:
    title_text.set_text('Surface at time %s' % timedelta(seconds=fgout1.t))

    # reset eta and B in plan-view plots to current state:

    eta_water = where(fgout1.h > 0, fgout1.eta, nan)
    eta_plot1.set_array(flipud(eta_water.T))
    B_plot1.set_array(flipud(fgout1.B.T))

    eta_water2 = where(fgout2.h > 0, fgout2.eta, nan)
    #h_water2 = where(fgout2.h > 0.02, fgout2.h, nan)
    h_water2 = where(fgout2.h > dry_depth, fgout2.h, nan)
    zeta_water2 = where(fgout2.B>0, h_water2, eta_water2)
    #eta_water = where(fgout2.X > -156.184, eta_water, nan)
    eta_plot2.set_array(flipud(zeta_water2.T))
    #B_plot2.set_array(flipud(fgout2.B.T))

    # update debris points:
    xd_water,yd_water = debris_tools.get_debris_xy(fgout2.t, debris_paths,
                                                   dbnos_water)
    xd_land,yd_land = debris_tools.get_debris_xy(fgout2.t, debris_paths,
                                                   dbnos_land)
    points_water.set_data(xd_water,yd_water)
    points_land.set_data(xd_land,yd_land)
    

    # update transects:
    
    Btrans1, etatrans1 = extract_transect(fgout1)
    Btrans2, etatrans2 = extract_transect(fgout2)

    Btrans, etatrans = Btrans2, etatrans2

    Btrans_plot.set_data(xtrans,Btrans)
    etatrans_plot.set_data(xtrans,etatrans)
    

    #update the PolyCollections for fill_between plots:             
    dummy = axdummy.fill_between(xtrans, Btrans-1e4, Btrans, 
                                      color=[.5,1,.5,1])
    dp = dummy.get_paths()[0]
    dummy.remove()
    Bfill_plot.set_paths([dp.vertices])

    dummy = axdummy.fill_between(xtrans, Btrans, etatrans, 
                                      color=[.5,.5,1,1])
    dp = dummy.get_paths()[0]
    dummy.remove()
    etafill_plot.set_paths([dp.vertices])
    
    update_artists = (B_plot1, eta_plot1, B_plot2, eta_plot2,
                      Bfill_plot, etafill_plot, Btrans_plot, etatrans_plot,
                      title_text, points_water, points_land)
    return update_artists

def plot_fgframe(fgframeno, save_png=False):
    """
    Convenience function for plotting one frame.
    But if you use this function in IPython and then try to make the animation,
    it may get into an infinite loop (not sure why).  Close the figure to abort.
    """
    update(fgframeno, *update_artists)
    
    if save_png:
        fname = 'fgout_frame%s.png' % str(fgframeno).zfill(4)
        savefig(fname, bbox_inches='tight')
        print('Created ',fname)

                

def make_anim():
    print('Making anim...')
    anim = animation.FuncAnimation(fig, update,
                                   frames=range(len(fgframes1)), 
                                   fargs=update_artists,
                                   interval=200, blit=True)
    return anim

if __name__ == '__main__':
    
    anim = make_anim()
    
    rundir = os.getcwd()
    if 'mmfs1' in rundir:
        # on hyak:
        scrdir = rundir.replace('mmfs1/home','gscratch/tsunami')
        outdir = os.path.join(scrdir, '_output')
    plotdir = outdir.replace('_output','_plots')

    if 0:
        fgout_plotdir = plotdir + '/animations'
        os.system('mkdir -p %s' % fgout_plotdir)
    else:
        fgout_plotdir = '.'

    # Output files:
    name = 'animation_pond%sm' % mstr

    fname_mp4 = os.path.join(fgout_plotdir, name + '.mp4')
    #fname_html = os.path.join(fgout_plotdir, name + '.html')
    fname_html = None
    
    
    if fname_mp4:
        fps = 5
        print('Making mp4...')
        writer = animation.writers['ffmpeg'](fps=fps)
        anim.save(fname_mp4, writer=writer)
        print("Created %s" % fname_mp4)

    if fname_html:
        # html version:
        fname_html = name + '.html'
        animation_tools.make_html(anim, file_name=fname_html, title=name)
