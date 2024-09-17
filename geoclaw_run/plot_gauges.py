
from pylab import *
import numpy
import clawpack.pyclaw.gauges as gauges


gs = 5.  # grounding speed


def plot_gauge(gaugeno,outdir,plot_color):
    gauge = gauges.GaugeSolution(gaugeno, outdir)
    t = gauge.t / 60.   # convert to minutes
    q = gauge.q
    h,hu,hv = q[:3,:]
    u = numpy.divide(hu, h, where=h>0, out=zeros(hu.shape))
    v = numpy.divide(hv, h, where=h>0, out=zeros(hv.shape))
    s = sqrt(u**2 + v**2)
    mflux = h*s
    eta = q[-1,:]

    subplot(211)
    plot(t, eta, color=plot_color,label='Water surface, Gauge %s' % gaugeno)
    xlabel('')
    ylabel('Elevation above MSL (m)', fontsize=12)
    xlim(0,20)
    ylim(0,11)
    xticks(fontsize=12)
    yticks(fontsize=12)
    grid(True)
    #title_str = 'Gauge %s for %s m slump' % (gaugeno,int(mstr))
    title_str = 'Time series at gauges for %s m slump' % int(mstr)
    title(title_str, fontsize=16)
    if plot_color=='r':
        plot([0,20],[8,8],'k-.',label='Wedge clast 14 elevation')
        plot([0,20],[3.2,3.2],'g--',label='Wedge clast 11 elevation')
    legend(loc='upper right',framealpha=1,fontsize=12)

    subplot(212)
    plot(t, s, color=plot_color,label='Speed, Gauge %s' % gaugeno)
    xlabel('time (minutes)', fontsize=12)
    ylabel('speed (m/s)', fontsize=12)
    xlim(0,20)
    ylim(0,13)
    xticks(fontsize=12)
    yticks(fontsize=12)
    grid(True)
    title('')

    if plot_color=='r':
        plot([0,20],[gs,gs],'k--',label='Grounding speed')
    legend(loc='upper right',framealpha=1,fontsize=12)

    tight_layout()



figure(400, figsize=(8,6))

for mstr in ['05','10','15']:  # loop over all 3 amplitudes
#for mstr in ['10']:
    clf()
    name = 'gauges_bouss_%sm' % mstr
    #gaugenos = [101,102,103]; plot_colors = ['b','r','g']
    gaugenos = [101,103]; plot_colors = ['b','r']
    for k,gaugeno in enumerate(gaugenos):
        #name = 'gauge%s_bouss_%sm' % (str(gaugeno).zfill(2), mstr)
        outdir = '_output_%sm' % mstr
        plot_gauge(gaugeno,outdir,plot_colors[k])

    fname = '%s.jpg' % name
    savefig(fname, bbox_inches='tight')
    print('Created %s' % fname)
    
