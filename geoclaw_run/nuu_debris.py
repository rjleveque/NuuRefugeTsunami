from pylab import *

def compute_debris_paths(mstr, tfinal, dxd=0.0001):

    """
    Load in fgout frames up to time tfinal from the tsunami simulation using
    the slump amplitude specified by the string mstr.
    
    Initial debris locations are on a regular grid with spacing dxd in both
    the x and y direction, over a portion of the beach region as specified in
    this routine. 
    dxd = 0.0001 was used for the figures in the paper,
    dxd = 0.0002 was used for making the 3D PyVista animations.
    
    Code from debris_tools is used to move the debris based on the velocities
    loaded from the fgout frames.
    """
    import debris_tools
    #from clawpack.geoclaw import debris_tools
    from clawpack.geoclaw import fgout_tools
    from scipy.interpolate import interp1d

    outdir = '_output_%sm' % mstr
    output_format = 'binary32'
    qmap = 'geoclaw-bouss'  # defines mapping for fgout variables

    fgno = 5
    fgout_grid2 = fgout_tools.FGoutGrid(fgno, outdir, output_format, qmap)
    fgout_grid2.read_fgout_grids_data(fgno)
    #fgout_grid2.read_fgout_grids_data_pre511(fgno)
    print('Looking for output in ',outdir)

    # use all frames up to time tfinal:
    fgframes2 = [n+1 for n in range(len(fgout_grid2.times)) \
                   if fgout_grid2.times[n] <= tfinal]
                   
    print('In nuu_debris, fgout grid 5, using %i frames up to frame %i, t = %i sec' \
        % (len(fgframes2), fgframes2[-1], fgout_grid2.times[fgframes2[-1]-1]))    
                      
    fgframe2 = fgframes2[0] # initial frame for fgout grid 2 
    fgout2 = fgout_grid2.read_frame(fgframe2)

    # Deterime time t of first fgout frame, to initialize particles
    t0 = fgout2.t
    fgout_extent = fgout2.extent_edges

    # initial topography (used to determine land / water points):
    # assumed to be same as B in fgout frames (no subsidence)
    B0_fcn = fgout_tools.make_fgout_fcn_xy(fgout2, 'B')

    # define beach region for placing initial particles only on beach:
    
    xy_beach = array([(-160, 20.62695),
                      (-156.1804431818182, 20.62694805194805),
                      (-156.17894318181817, 20.626974025974025),
                      (-156.17794967532467, 20.626603896103894),
                      (-156.1774237012987, 20.62638961038961),
                      (-156.175, 20.62639)])
     
    ybeach = interp1d(xy_beach[:,0], xy_beach[:,1], 
                       kind='linear', bounds_error=False, fill_value=nan)
                       
    # Initialize debris_paths dictionary and set
    # initial debris particle locations (x0,y0) at time t0.
    # Require a list of dbnos and each array 
    #     debris_paths[dbno]
    # in the dictionary is a 2d array with a single row [t0, x0, y0] to start.

    debris_paths = {}
    dbnos_water = []
    dbnos_land = []

    # same for all debris particles:
    grounding_depth = 0.1
    grounding_speed = 5.
    #drag_factor = None
    drag_factor = 10.
            
    # set initial velocities to 0 (or may want to interpolate from fgout0 if t0>0?)
    u0 = 0.
    v0 = 0.

    x1,x2 = -156.181, -156.177
    y1,y2 = 20.625, 20.6273
    #dxd = 0.0001 # now input argument -- spacing between debris particles

    dyd = dxd * cos(21*pi/180)
    xd = arange(x1, x2, dxd)
    yd = arange(y1, y2, dyd)
    mxd = len(xd)
    myd = len(yd)
    print('will sample array of %i by %i points for particles' \
            % (mxd,myd))

    for i in range(len(xd)):
        x0 = xd[i]
        for j in range(len(yd)):
            y0 = yd[j]
            B0xy = B0_fcn(x0,y0)
            #if B0xy < -10 or B0xy > 5:
            if B0xy < -5 or B0xy > 4:
                continue # skip point unless near shore
            if x0 > -156.1778 and y0 > 20.6265:
                #continue # skip point if in pond
                pass
            if y0 > 20.6273:
                continue # skip points inland
            if y0 > ybeach(x0):
                continue # skip points too far inland
            dbno = 1000*i + j
            db = array([[t0, x0, y0, u0, v0]])
            debris_paths[dbno] = db
            
            if B0xy < 0:
                dbnos_water.append(dbno)
            elif B0xy > 0:
                dbnos_land.append(dbno)

    dbnos = dbnos_water + dbnos_land
    print('Created %i initial debris particles' % len(dbnos))
    print('    %i offshore,  %i onshore' % (len(dbnos_water),len(dbnos_land)))

    # Compute debris path for each particle by using all the fgout frames
    # in the list fgframes (first frame should be one used to set t0 above):

    debris_paths = debris_tools.make_debris_paths(fgframes2, fgout_grid2, 
                                debris_paths, dbnos, drag_factor, grounding_depth, grounding_speed)

    # done computing debris paths
    # ------------------------------
    return debris_paths, dbnos_water, dbnos_land
