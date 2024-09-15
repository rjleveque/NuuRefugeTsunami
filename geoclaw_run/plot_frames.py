
import sys
if 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('Agg')  # Use an image backend

from pylab import *
import os
from clawpack.visclaw import plottools, geoplot, gridtools
from clawpack.visclaw import colormaps
from matplotlib import animation, colors
from clawpack.geoclaw import fgout_tools
from datetime import timedelta 
import debris_tools
from scipy.interpolate import interp1d
from nuu_debris import compute_debris_paths

graphics_dir = os.path.abspath('../topo')
hillshade_image = load(graphics_dir + '/NuuHillshadeOffshore.npy')
hillshade_extent = loadtxt(graphics_dir + '/NuuHillshadeOffshore_extent.txt')
figure_format = 'png'

mstr = '05'  # slump amplitude (rerun code for 05, 10, 15)
outdir = '_output_%sm' % mstr
output_format = 'binary32'
qmap = 'geoclaw-bouss'  # defines mapping for fgout variables

# Instantiate object for reading fgout frames:
fgout_grid1 = fgout_tools.FGoutGrid(2, outdir, output_format, qmap)
fgout_grid1.read_fgout_grids_data(2)
fgout_grid2 = fgout_tools.FGoutGrid(5, outdir, output_format, qmap)
fgout_grid2.read_fgout_grids_data(5)

#fgframes1 = range(1,122,1) # 20 minutes
#fgframes2 = fgframes1

tfinal = 20*60.  # use fgout frames up to this time (or end of computation)

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

fgout1 = fgout_grid1.read_frame(fgframe1)
fgout2 = fgout_grid2.read_frame(fgframe2)

    
def plot_nuu_fgframe(k, debris_paths=None, dbnos_water=None, dbnos_land=None):

    if debris_paths is None:
        print('Need to compute debris_paths via:')
        print('   debris_paths, dbnos_water, dbnos_land = compute_debris_paths(mstr)')

    #onecolor = [0.3,0.3,1]
    onecolor = [0.7,0.7,1]
    cmap_onecolor = colormaps.make_colormap({0.:onecolor, 1.:onecolor})

    fgout2 = fgout_grid2.read_frame(k)
    print('Read fgout2 at time %.1f' % fgout2.t)
        
    fig = figure(2, figsize=(13,8))
    clf()
        
    axins = axes()
    #axins.imshow(GE_image, extent=GE_extent)
    axins.imshow(hillshade_image, extent=hillshade_extent)

    #B = ma.masked_where(fgout2.h < 0.01, fgout2.B)
    B = nan*fgout2.B  # don't show topo since using GE_image
    B_plot2 = axins.imshow(flipud(B.T), extent=fgout2.extent_edges,
           cmap=geoplot.googleearth_transparent)
    B_plot2.set_clim(0,8)

    eta_water2 = where(fgout2.h > 0, fgout2.eta, nan)
    h_water2 = where(fgout2.h > 0.02, fgout2.h, nan)
    zeta_water2 = where(fgout2.B>0, h_water2, eta_water2)
    eta_plot2 = axins.imshow(flipud(zeta_water2.T), extent=fgout2.extent_edges,
           cmap=cmap_onecolor)
           #cmap=geoplot.tsunami_colormap)
    eta_plot2.set_clim(-5,5)

    # sub region of the original image
    x1, x2, y1, y2 = [-156.1825,-156.175,20.625,20.629]
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    axins.set_xticks(arange(-156.182,-156.175,.001),rotation=20)
    axins.ticklabel_format(useOffset=False)
    axins.set_xlim(-156.181, -156.175)
    
    axins.set_title("Inundation from %s m slump at time %s" \
                    % (int(mstr),timedelta(seconds=fgout2.t)), fontsize=15)
    axins.contour(fgout2.X,fgout2.Y,fgout2.B,[1],colors='g')

    # show where this region is on the Ocean plot:
    #ax.indicate_inset_zoom(axins, edgecolor="black")
    
    xlabel('Longitude', fontsize=12)
    ylabel('Latitude', fontsize=12)
    xticks(fontsize=11,rotation=20)
    yticks(fontsize=11)

    # -----------
    # add debris:
    
    if debris_paths:
                
        xd_water,yd_water = debris_tools.get_debris_xy(fgout2.t, 
                                            debris_paths, dbnos_water)
        xd_land,yd_land = debris_tools.get_debris_xy(fgout2.t, 
                                            debris_paths, dbnos_land)

        #color_rock = 'yellow'; color_coral = 'w'
        color_rock = 'k'; color_coral = 'r'
        points_water, = axins.plot(xd_water, yd_water,
                                color=color_coral,linestyle='',
                                marker='.',markersize=8)
        points_land, = axins.plot(xd_land, yd_land,
                                color=color_rock,linestyle='',
                                marker='s',markersize=5)
                            


def plot_ocean_fgframe(k):     
    fig = figure(1, figsize=(13,8))
    clf()

    # =========
    # Ocean:
    # =========

    fgout1 = fgout_grid1.read_frame(k)
    print('Read fgout1 at time %.1f' % fgout1.t)
    
    #plot_extent = [-156.4, -156., 20.4, 20.7] # for animation
    plot_extent = [-156.26, -156.12, 20.54, 20.64]

    ax = axes()
    B_plot1 = ax.imshow(flipud(fgout1.B.T), extent=fgout1.extent_edges,
           cmap=geoplot.land1_colormap)

    B_plot1.set_clim(0,1500)

    eta_water = where(fgout1.h > 0, fgout1.eta, nan)
    eta_plot1 = ax.imshow(flipud(eta_water.T), extent=fgout1.extent_edges,
           cmap=geoplot.tsunami_colormap)

    eta_plot1.set_clim(-1,1)
    axis(plot_extent)
    title_text = ax.set_title('Surface at time %s' % timedelta(seconds=fgout1.t),
                    fontsize=15)


    if 1:
        cb = colorbar(eta_plot1, extend='both', shrink=0.5)
        cb.set_label('meters')

def annotate_plot():
    # plot transect and gauges
    x1trans = -156.1766; y1trans = 20.6280
    x2trans = -156.17812; y2trans = 20.6257
    xg1,yg1 = -1.5617660000e+02,   2.0627200000e+01
    xg2,yg2 = -1.5617800000e+02,   2.0626400000e+01
    xg102,yg102 = -156.17713,   20.62722
    xg101,yg101 = -156.17777,   20.62622
    #xwedge = -156.176917; ywedge = 20.62752
    xg103,yg103 = -156.176917, 20.62752
    plot([x1trans,x2trans],[y1trans,y2trans],'k')
    #plot([xg102],[yg102],'bo',markersize=5)
    xeps = 0.00005
    yeps = 0.0001
    #text(xg102+xeps,yg102-yeps,'Gauge 102',color='b',fontsize=8,
    #     backgroundcolor='w', alpha=1)
    plot([xg101],[yg101],'bo',markersize=5)
    #text(xg101+xeps,yg101-yeps,'Gauge 101',color='b',fontsize=8,
    #     backgroundcolor='w', alpha=1)
    eps = 0.0001
    annotate("Gauge 101",
            xy=(xg101-eps,yg101-0.1*eps), xycoords='data',
            xytext=(xg101-10*eps, yg101-2*eps), textcoords='data',
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),
            fontsize=10, alpha=1, #backgroundcolor='w', 
            bbox={'facecolor':'w','edgecolor':'k','lw':0.5})
    plot([xg103],[yg103],'bo',markersize=5)
    #text(xg103+xeps,yg103-yeps,'Gauge 103',color='b',fontsize=8,
    #     backgroundcolor='w', alpha=1)
    annotate("Gauge 103",
            xy=(xg103,yg103-0.4*eps), xycoords='data',
            xytext=(xg103-eps, yg103-4*eps), textcoords='data',
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),
            fontsize=10, alpha=1,
            bbox={'facecolor':'w','edgecolor':'k','lw':0.5})

    xw = -156.17656 - 2.8e-4; yw = 20.627657
    plot([xw],[yw],color='w',linestyle='',
         marker='d',fillstyle='top',markersize=9,
         label='Wedge clast 11')
        
    xw = -156.17656; yw = 20.627657
    plot([xw],[yw],color='w',linestyle='',
         marker='d',fillstyle='bottom',markersize=9,
         label='Wedge clast 14')

    legend(loc='upper right',facecolor=[.2,.2,.2],labelcolor='w',framealpha=1)


def make_initial_particle_plot(debris_paths, dbnos_water, dbnos_land):
    # initial particles
    plot_nuu_fgframe(10,debris_paths, dbnos_water, dbnos_land)
    title("Initial particle distribution",fontsize=15)
    fname = 'nuu_initial_particles.%s' % figure_format
    savefig(fname, bbox_inches='tight')
    print('Created ',fname)
                    
def make_gauge_transect_plot():
    # annotated at initial time 
    plot_nuu_fgframe(10) # with no particles
    annotate_plot()
    title("Nu'u Refuge with pond, gauges, and transect",fontsize=15)
    fname = 'nuu_gauges_transect.%s' % figure_format
    savefig(fname, bbox_inches='tight')
    print('Created ',fname)
     
def make_plots(debris_paths, dbnos_water, dbnos_land):

    make_gauge_transect_plot()
    
    for k in [10,22,121]:
        if 0:
            plot_ocean_fgframe(k)
            fname = 'ocean_frame%s_%sm.%s' \
                    % (str(k).zfill(3), mstr, figure_format)
            savefig(fname, bbox_inches='tight')
            print('Created ',fname)
        
        k2 = fgframes2[k-1]
        print('+++ k2 = %i, t = %.1f' % (k2,fgout_grid2.times[k2-1]))
        plot_nuu_fgframe(k2, debris_paths, dbnos_water, dbnos_land)
            
        fname = 'nuu_frame%s_%sm.%s' \
                % (str(k).zfill(3), mstr, figure_format)
        savefig(fname, bbox_inches='tight')
        print('Created ',fname)
    
    
if __name__=='__main__':

    print('Making plots for slump amplitude %s' % mstr)
    debris_paths, dbnos_water, dbnos_land = compute_debris_paths(mstr, tfinal)
    make_plots(debris_paths, dbnos_water, dbnos_land)
    
    make_gauge_transect_plot()
    make_initial_particle_plot(debris_paths, dbnos_water, dbnos_land)
    
    if 0:
        # for debugging:
        dbkeys = array(list(debris_paths.keys()))
        idbkeys = array([floor_divide(d,1000) for d in dbkeys])
        jdbkeys = array([remainder(d,1000) for d in dbkeys])
        #figure()
        #plot(idbkeys,jdbkeys,'o',markersize=3)
