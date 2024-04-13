
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import numpy as np
import matplotlib.pyplot as plt
import os

from clawpack.geoclaw import topotools

if 0:
    image = plt.imread('GE_PA2.png')

    def background_image(current_data):
        from pylab import imshow
        extent = [-123.5,-123.33,48.10,48.17]
        imshow(image,extent=extent)

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps, geoplot
    from numpy import linspace

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'


    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)
    

    def timeformat(t):
        from numpy import mod
        hours = int(t/3600.)
        tmin = mod(t,3600.)
        min = int(tmin/60.)
        sec = int(mod(tmin,60.))
        timestr = '%s:%s:%s' % (hours,str(min).zfill(2),str(sec).zfill(2))
        return timestr

    def title_hours(current_data):
        from pylab import title
        t = current_data.t
        timestr = timeformat(t)
        title('%s after earthquake' % timestr)

    def surface_or_depth_feet(current_data):
        sdf = geoplot.surface_or_depth(current_data)
        sdf *= 1./0.3048  # convert to feet
        return sdf
    

    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=0)
    plotfigure.show = False
    plotfigure.kwargs = {'figsize':(8,7)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = False
    #plotaxes.xlimits = [-170,-66]
    #plotaxes.ylimits = [-48,32]

    def fixup(current_data):
        from pylab import title, ticklabel_format, gca, cos, pi
        import pylab
        #addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        title('Surface at %4.2f hours' % t, fontsize=10)
        pylab.xticks(fontsize=10)
        pylab.yticks(fontsize=10)
        ticklabel_format(useOffset=False)
        pylab.xticks(rotation=20)
        gca().set_aspect(1./cos(10*pi/180.))
        title_hours(current_data)
    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface
    #plotitem.plot_var = surface_or_depth_feet
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -5. 
    plotitem.pcolor_cmax = 5. 
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.colorbar_kwargs = {'label':'meters'}
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land1_colormap
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 2000.0
    plotitem.add_colorbar = False
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    #plotaxes.xlimits = [-120,-60]
    #plotaxes.ylimits = [-60,0]

    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = linspace(-3000,-3000,1)
    plotitem.amr_contour_colors = ['y']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':2}
    plotitem.amr_contour_show = [1,0,0]  
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


    #-----------------------------------------
    # Figure for 36sec Topo Area
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='RegionTopo36_30sec_1sec', figno=1)
    plotfigure.show = False
    plotfigure.kwargs = {'figsize':(8,7)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = False
    plotaxes.xlimits = [-164.994, -150.006]
    plotaxes.ylimits = [15.006,24.994]

    def fixup(current_data):
        from pylab import title, ticklabel_format, gca, cos, pi
        import pylab
        #addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        title('Surface at %4.2f hours' % t, fontsize=10)
        pylab.xticks(fontsize=10)
        pylab.yticks(fontsize=10)
        ticklabel_format(useOffset=False)
        pylab.xticks(rotation=20)
        gca().set_aspect(1./cos(20*pi/180.))
        title_hours(current_data)
    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface
    #plotitem.plot_var = surface_or_depth_feet
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -5. 
    plotitem.pcolor_cmax = 5. 
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.colorbar_kwargs = {'label':'meters'}
    plotitem.celledges_show = 0
    plotitem.amr_patchedges_show = [0]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land1_colormap
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 2000.0
    plotitem.add_colorbar = False
    plotitem.celledges_show = 0
    plotitem.amr_patchedges_show = [0]
    #plotaxes.xlimits = [-120,-60]
    #plotaxes.ylimits = [-60,0]

    #-----------------------------------------
    # Figure for Greater Islands 30sec
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Greater Islands_30sec', figno=2)
    plotfigure.show = False
    plotfigure.kwargs = {'figsize':(8,7)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = False
    plotaxes.xlimits = [-160.7, -154.54]
    plotaxes.ylimits = [18.5,22.6]

    def fixup(current_data):
        from pylab import title, ticklabel_format, gca, cos, pi
        import pylab
        #addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        title('Surface at %4.2f hours' % t, fontsize=10)
        pylab.xticks(fontsize=10)
        pylab.yticks(fontsize=10)
        ticklabel_format(useOffset=False)
        pylab.xticks(rotation=20)
        gca().set_aspect(1./cos(20*pi/180.))
        title_hours(current_data)
    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface
    #plotitem.plot_var = surface_or_depth_feet
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -5. 
    plotitem.pcolor_cmax = 5. 
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.colorbar_kwargs = {'label':'meters'}
    plotitem.celledges_show = 0
    plotitem.amr_patchedges_show = [0]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land1_colormap
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 2000.0
    plotitem.add_colorbar = False
    plotitem.celledges_show = 0
    plotitem.amr_patchedges_show = [0]
    #plotaxes.xlimits = [-120,-60]
    #plotaxes.ylimits = [-60,0]

    #-----------------------------------------
    # Figure for Maui 30sec to 6sec
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Maui_30sec_6sec', figno=3)
    #plotfigure.show = False
    plotfigure.kwargs = {'figsize':(8,7)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = False
    #plotaxes.xlimits = [-156.85, -155.8]
    #plotaxes.ylimits = [20.5,21.19]
    #plotaxes.xlimits = [-156.4, -156.0]
    #plotaxes.ylimits = [20.4,20.7]
    plotaxes.xlimits = [-156.3, -156.1]
    plotaxes.ylimits = [20.46,20.64]

    def fixup(current_data):
        from pylab import title, ticklabel_format, gca, cos, pi
        import pylab
        #addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        title('Surface at %4.2f hours' % t, fontsize=10)
        pylab.xticks(fontsize=10)
        pylab.yticks(fontsize=10)
        ticklabel_format(useOffset=False)
        pylab.xticks(rotation=20)
        gca().set_aspect(1./cos(20*pi/180.))
        title_hours(current_data)
    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface_or_depth
    #plotitem.plot_var = surface_or_depth_feet
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -2. 
    plotitem.pcolor_cmax = 2. 
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.colorbar_kwargs = {'label':'meters'}
    plotitem.celledges_show = 0
    plotitem.amr_patchedges_show = [0]
    #plotitem.amr_patchedges_show = [0,0,1,1]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land1_colormap
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 2000.0
    plotitem.add_colorbar = False
    plotitem.celledges_show = 0
    plotitem.amr_patchedges_show = [0]
    #plotaxes.xlimits = [-120,-60]
    #plotaxes.ylimits = [-60,0]

    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = [0]  #[0,2,4,6,8]
    plotitem.amr_contour_colors = ['yellow']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':1.0}
    plotitem.amr_contour_show = [0,0,1,0,0,0,0,0,1]  
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


    #-----------------------------------------
    # Figure for Maui 2sec
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Maui_2sec', figno=4)
    plotfigure.show = False
    plotfigure.kwargs = {'figsize':(8,7)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = False
    plotaxes.xlimits = [-156.206, -156.162]
    plotaxes.ylimits = [20.612,20.629]

    def fixup(current_data):
        from pylab import title, ticklabel_format, gca, cos, pi
        import pylab
        #addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        title('Surface at %4.2f hours' % t, fontsize=10)
        pylab.xticks(fontsize=10)
        pylab.yticks(fontsize=10)
        ticklabel_format(useOffset=False)
        pylab.xticks(rotation=20)
        gca().set_aspect(1./cos(20*pi/180.))
        title_hours(current_data)
    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface
    #plotitem.plot_var = surface_or_depth_feet
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -5. 
    plotitem.pcolor_cmax = 5. 
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.colorbar_kwargs = {'label':'meters'}
    plotitem.celledges_show = 0
    plotitem.amr_patchedges_show = [0]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land1_colormap
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 2000.0
    plotitem.add_colorbar = False
    plotitem.celledges_show = 0
    plotitem.amr_patchedges_show = [0]
    #plotaxes.xlimits = [-120,-60]
    #plotaxes.ylimits = [-60,0]


    #-----------------------------------------
    # Figure for Maui_13sec
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Maui_13sec', figno=5)
    #plotfigure.show = False
    plotfigure.kwargs = {'figsize':(10,7)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = False
    plotaxes.xlimits = [-156.189, -156.17]
    #plotaxes.ylimits = [20.618,20.628]
    plotaxes.ylimits = [20.618,20.6295]

    def fixup(current_data):
        from pylab import title, ticklabel_format, gca, cos, pi
        import pylab
        #addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        title('Surface at %4.2f hours' % t, fontsize=10)
        pylab.xticks(fontsize=10)
        pylab.yticks(fontsize=10)
        ticklabel_format(useOffset=False)
        pylab.xticks(rotation=20)
        gca().set_aspect(1./cos(20*pi/180.))
        title_hours(current_data)
    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface_or_depth
    #plotitem.plot_var = surface_or_depth_feet
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -2.
    plotitem.pcolor_cmax = 2. 
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.colorbar_kwargs = {'label':'meters'}
    plotitem.celledges_show = 0
    plotitem.amr_patchedges_show = [0]
    plotitem.amr_data_show = [0,0,1,0,0,0,0,0,1]

    # Land
    if 0:
        plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
        plotitem.plot_var = geoplot.land
        plotitem.pcolor_cmap = geoplot.land1_colormap
        plotitem.pcolor_cmin = 0.0
        plotitem.pcolor_cmax = 40.0
        plotitem.add_colorbar = False
        plotitem.celledges_show = 0
        plotitem.amr_patchedges_show = [0]
        #plotaxes.xlimits = [-120,-60]
        #plotaxes.ylimits = [-60,0]

    if 1:
        plotitem = plotaxes.new_plotitem(plot_type='2d_hillshade')
        plotitem.plot_var = geoplot.land
        plotitem.hillshade_cmap = 'gray'
        plotitem.hillshade_vertical_exaggeration = 2
        plotitem.hillshade_azimuth_degree =  315
        plotitem.hillshade_altitude_degree = 45
        plotitem.hillshade_latlon = True
        #plotitem.add_colorbar = False

    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = [0,2,4,6,8]
    plotitem.amr_contour_colors = ['yellow']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':1.0}
    plotitem.amr_contour_show = [0,0,1,0,0,0,0,0,1]  
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------

    time_scale = 1./60. # minutes
    time_label = ''
    tlimits = [0,20] # minutes

    plotfigure = plotdata.new_plotfigure(name='Gauges',figno=300,type='each_gauge')
    plotfigure.kwargs = {'figsize':(8,8)}
    # plotfigure.kwargs = {'figsize':(12,12)}
    #plotfigure.clf_each_gauge = True
    plotfigure.clf_each_gauge = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(3,1,1)'
    plotaxes.time_scale = time_scale
    plotaxes.time_label = time_label
    plotaxes.xlimits = tlimits
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Quantities of interest'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')

    # Plot eta, the amount above MHW everywhere
    def eta(current_data):
        #from numpy import nan, logical_and
        q = current_data.q
        eta = q[-1,:]  # eta for SWE or Bouss
        return eta
    plotitem.plot_var = eta
    plotitem.plotstyle = 'b-'

    # Plot h:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(3,1,2)'
    plotaxes.time_scale = time_scale
    plotaxes.time_label = ''
    plotaxes.xlimits = tlimits
    plotaxes.ylimits = 'auto'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    def h(current_data):
        q = current_data.q
        h = q[0,:]
        return h
    plotitem.plot_var = h
    plotitem.plotstyle = 'b-'

    # Plot current speed:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(3,1,3)'
    plotaxes.time_scale = time_scale
    plotaxes.time_label = ''
    plotaxes.xlimits = tlimits
    plotaxes.ylimits = 'auto'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    def speed(current_data):
        from numpy import where, sqrt
        h = current_data.q[0,:]
        h = where(h>0.01, h, 1.e6)
        u = current_data.q[1,:] / h
        v = current_data.q[2,:] / h
        s = sqrt(u**2 + v**2)
        return s
    plotitem.plot_var = speed
    plotitem.plotstyle = 'b-'


    def aa(current_data):

        from pylab import clf, subplot,xlabel,ylabel,title,grid,where,ylim
        q = current_data.q
        gaugeno = current_data.gaugeno
        gaugesoln = current_data.plotdata.getgauge(gaugeno)
        level = gaugesoln.level
        max_level = level.max()
        index = where(level == max_level)[0]
        hmax = q[0,index].max()
        hmin = q[0,index].min()
        etamax = q[-1,index].max()
        etamin = q[-1,index].min()

        subplot(311)
        ylabel('Surface (eta=h+B), (m)')
        title('Gauge %i: Max h = %.2f, Max eta = %.2f' \
            % (gaugeno,hmax,etamax))
        grid(color='k', linestyle='dotted', linewidth=0.5)
        buffer = 0.1*(etamax - etamin)
        ylim(etamin-buffer,etamax+buffer)

        subplot(312)
        ylabel('h (m)')
        title('')
        grid(color='k', linestyle='dotted', linewidth=0.5)
        buffer = 0.1*(hmax - hmin)
        ylim(hmin-buffer,hmax+buffer)

        subplot(313)
        xlabel('time (minutes)')
        ylabel('s, speed (m/s)')
        title('')
        grid(color='k', linestyle='dotted', linewidth=0.5)

    plotaxes.afteraxes = aa


    #otherfigure = plotdata.new_otherfigure(name='animations',
    #                fname='animations')

    if 0:
        #-------------------------------
        # Other Figures for this Site 
        #-------------------------------
        other_figures_dir    =  plotdata.plotdir + '/_other_figures'
        print('other_figures_dir = ',other_figures_dir)
        if os.path.isdir(other_figures_dir):
            files = os.listdir(other_figures_dir)
            print('files: ',files,' END of files')
            if len(files) > 0:
                files.sort()
                print('%i files in this directory' % len(files))
                for filename in files:
                    print('\nfilename=',filename)
                    path='_other_figures/'+filename
                    otherfigure = plotdata.new_otherfigure(name=filename,fname=path)
                    print('Added other figure: ',path)
            else:
                print('No files in this directory')
        else:
            print('*** directory not found, will not add to index')

    #-----------------------------------------
    # Figures for fgmax plots
    #-----------------------------------------
    # Note: You need to move fgmax png files into _plots/_other_figures after 
    # creating them, e.g., by running the process_fgmax notebook or script.
    # The lines below just create links to these figures from _PlotIndex.html 

    if 0:
        # included by listdir version above:
        otherfigure = plotdata.new_otherfigure(name='max depth',
                        fname='_other_figures/h_onshore.png')

        otherfigure = plotdata.new_otherfigure(name='max speed',
                        fname='_other_figures/speed.png')
            
    loc   = 'LandTrust'
    event = 'chile-moreno'
    fname_kmz = 'fgmax_results_%s_%s.kmz' % (loc,event)
    otherfigure = plotdata.new_otherfigure(name=fname_kmz,
                    fname='_other_figures/kmlfiles/%s' % fname_kmz)
        

    # Plots of timing (CPU and wall time):

    def make_timing_plots(plotdata):
        import os
        from clawpack.visclaw import plot_timing_stats
        try:
            timing_plotdir = plotdata.plotdir + '/timing_figures'
            os.system('mkdir -p %s' % timing_plotdir)
            units = {'comptime':'hours', 'simtime':'hours', 'cell':'billions'}
            plot_timing_stats.make_plots(outdir=plotdata.outdir, make_pngs=True,
                                          plotdir=timing_plotdir, units=units)
            os.system('cp %s/timing.* %s' % (plotdata.outdir, timing_plotdir))
        except:
            print('*** Error making timing plots')

    # create a link to this webpage from _PlotIndex.html:
    otherfigure = plotdata.new_otherfigure(name='timing',
                    fname='timing_figures/timing.html')
    otherfigure.makefig = make_timing_plots

    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.parallel = True

    return plotdata

