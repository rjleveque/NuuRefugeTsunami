r"""
Functions for tracking particles or debris in a velocity field specified
by a set of fgout frames.

Note: The functions clawpack.geoclaw.fgout_tools.make_fgout_fcn_xy and
clawpack.geoclaw.fgout_tools.make_fgout_fcn_xyt are used
to create functions used to interpolate the fluid velocity u,v
from fgout frames to particle locations.

These are experimental tools that are still under development.
Adapted from https://github.com/rjleveque/geoclaw_debris
with the addition of grounding_speed.

"""

from pylab import *
import os
from clawpack.geoclaw.data import Rearth,DEG2RAD



def move_debris(dbnos, debris_paths, fgout1, fgout2, coordinate_system=2,
                drag_factor=None, grounding_depth=0., grounding_speed=0.):
    """
    For each dbno in dbnos: debris_paths[dbno] is a 2D array and it is assumed
    that the last row has the form [t1, x1, y1, u1, v1] with the location and
    velocity of this debris particle at time t1, which should equal fgout1.t.
    
    fgout1 and fgout2 should be two time frames of fgout data, each of class
    fgout_tools.FGoutFrame.  
    
    Compute the location and velocity of the debris particle at time fgout2.t,
    and append [t2, x2, y2, u2, v2] to the bottom of debris_paths[dbno].
    
    coordinate_system = 1 if (x,y) are in meters, 
                        2 if (x,y) are degrees longitude, latitude
    The velocities (u,v) are always in m/s.
    
    drag_factor and grounding_depth can be single values applying to all
    particles, or dictionaries with keys given by the elements of dbnos.
    
    If drag_factor[dbno]==None and grounding_depth[dbno]==0, then we passively
    advect the debris particle with the fluid velocity.  This is
    currently implemented using the 2-step explicit Runge-Kutta method:
    1. Interpolate fgout1.u and fgout1.v to (x1,y1) and move the particle
    over time dt/2 with this velocity, to obtain (xm,ym).
    2. Interpolate (u,v) in space and time to this midpoint location, and then
    use this velocity to move the particle over time dt from (x1,y1) to (x2,y2).
    
    If grounding_depth[dbno] > 0 and/or grounding_speed > 0,
    then the particle is assumed to be grounded
    whenever the local fluid depth h < grounding_depth[dbno]
    or the speed s < grounding_speed, where h and s are
    interpolated first to (x1,y1) and then to (xm,ym).
    
    If drag_factor[dbno] > 0, then 
    """
    
    from clawpack.geoclaw.fgout_tools import make_fgout_fcn_xyt
    
    t1 = fgout1.t
    t2 = fgout2.t
    dt = t2 - t1
    
    print('Moving debris over time dt = %g' % dt)
    print('       from t1 = %s to t2 = %.2f' % (t1,t2))

    h_fcn = make_fgout_fcn_xyt(fgout1, fgout2, 'h')
    u_fcn = make_fgout_fcn_xyt(fgout1, fgout2, 'u')
    v_fcn = make_fgout_fcn_xyt(fgout1, fgout2, 'v')

    for dbno in dbnos:
        debris_path = debris_paths[dbno]
        try:
            t1,xd1,yd1,ud1,vd1 = debris_path[-1,:]
        except:
            errmsg = 'on input, debris_path[%s] should be a 2d array ' % dbno \
                        + 'with at least one row'
            raise ValueError(errmsg)
            
        errmsg = '*** For dbno = %i, expected t1 = %.3f to equal fgout1.t = %.3f' \
                % (dbno, t1, fgout1.t)
        assert t1 == fgout1.t, errmsg
        
        try:
            gd = grounding_depth[dbno]  # if different for each particle
        except:
            gd = grounding_depth  # assume it's a scalar, same for all debris
            
        try:
            gs = grounding_speed[dbno]  # if different for each particle
        except:
            gs = grounding_speed  # assume it's a scalar, same for all debris
             
        try:
            df = drag_factor[dbno]  # if different for each particle
        except:
            df = drag_factor  # assume it's a scalar, same for all debris
        
        if coordinate_system == 2:
            # x,y in degrees, u,v in m/s
            # convert u,v to degrees/second:
            ud1 = ud1 / (Rearth*DEG2RAD * cos(DEG2RAD*yd1))
            vd1 = vd1 / (Rearth*DEG2RAD)
            
        # Assume ud1, vd1 were properly set at time fgout1.t, so no need
        # to check for grounding now.
        
        # Half time step with old velocities:
        xdm = xd1 + 0.5*dt*ud1
        ydm = yd1 + 0.5*dt*vd1
        
        tm = t1 + 0.5*dt  # t at midpoint in time
        
        # depth and fluid velocity at midpoint in time tm:
        hm = h_fcn(xdm,ydm,tm)
        um = u_fcn(xdm,ydm,tm)
        vm = v_fcn(xdm,ydm,tm)
        sm = sqrt(um**2 + vm**2)
    
        if hm < gd or sm < gs:
            # particle is grounded so velocities set to 0:
            udm = 0.
            vdm = 0.        
        elif df is None:
            # no drag factor and debris velocity = fluid velocity:
            udm = um
            vdm = vm
        else:
            # debris velocity (ud,vd) relaxes toward fluid velocity (u,v):
            # solve 
            #   d/dt ud = C * (u-ud)
            #   d/dt vd = C * (v-vd)
            # where C = df*sqrt((u-ud)**2 + (v-vd)**2).
            # Approximate this by setting C to value at t1, so simple
            # exponential decay with rate C:
            u1 = u_fcn(xd1,yd1,t1)
            v1 = v_fcn(xd1,yd1,t1)
            C = df*sqrt((u1-ud1)**2 + (v1-vd1)**2)
            udm = um - exp(-C*dt/2)*(um-ud1)
            vdm = vm - exp(-C*dt/2)*(vm-vd1)            
            
        if coordinate_system == 2:
            udm = udm / (Rearth*DEG2RAD * cos(DEG2RAD*ydm))
            vdm = vdm / (Rearth*DEG2RAD)
        
        # Take full time step with mid-point debris velocities:
        xd2 = xd1 + dt*udm
        yd2 = yd1 + dt*vdm
        
        # Depth and fluid velocity at final time t2 (to append to debris_path):
        h2 = h_fcn(xd2,yd2,t2)
        u2 = u_fcn(xd2,yd2,t2)
        v2 = v_fcn(xd2,yd2,t2)

        if h2 < gd:
            # particle is grounded so velocities set to 0:
            ud2 = 0.
            vd2 = 0.        
        elif df is None:
            # no drag factor and debris velocity = fluid velocity:
            ud2 = u2
            vd2 = v2
        else:
            # debris velocity (ud,vd) relaxes toward fluid velocity (u,v).
            # Take another half time step of decay from (udm,vdm),
            # now approximating C at the midpoint in time:
            C = df*sqrt((um-udm)**2 + (vm-vdm)**2)
            ud2 = u2 - exp(-C*dt/2)*(u2-udm)
            vd2 = v2 - exp(-C*dt/2)*(v2-vdm)    
        
        debris_paths[dbno] = vstack((debris_path, array([t2,xd2,yd2,ud2,vd2])))
        
    return debris_paths


            
def make_debris_paths(fgframes, fgout_grid, debris_paths, dbnos,
                      drag_factor=None, grounding_depth=0., grounding_speed=0.):
    """
    dbnos is a list of debris particle labels (integers) to operate on.
    debris_paths a dictionary indexed by integer dbno.
    Each element 
        debris_path = debris_paths[dbno] 
    is a 2D numpy array with at least one row and three columns t,x,y.
    
    The last row of debris_path[dbno], as passed in, defines the starting time
    and location of the particle dbno.
    
    This routine loops over all fgout frames in fgframes, reads in the frame
    for fgout grid number fgno as fgout2, 
    and then calls move_debris to move each particle from time
    fgout1.t to fgout2.t, where fgout1 is the previous frame.  The new time
    and location are added as a new row in the debris_path array.
    
    It is assumed that the 
    """
    
    fgout1 = fgout_grid.read_frame(fgframes[0])
    for dbno in dbnos:
        try:
            debris_path = debris_paths[dbno]
            t1,xd1,yd1,ud1,vd1 = debris_path[-1,:]
        except:
            errmsg = 'on input, debris_path[%s] should be a 2d array ' % dbno \
                        + 'with at least one row'
            raise ValueError(errmsg)
        errmsg = 'Time of first fgout frame %i is t = %g ' \
             % (fgframes[0],fgout1.t) \
             + 'This should agree with t value in last row of debris_paths[%s]'\
             % dbno
        assert t1 == fgout1.t, errmsg
    
    for fgframe in fgframes[1:]:
        print('Trying to read fgno=%i, fgframe=%i' % (fgout_grid.fgno,fgframe))
        try:
            fgout2 = fgout_grid.read_frame(fgframe)
        except:
            print('Could not read file, exiting loop')
            break
        debris_paths = move_debris(dbnos, debris_paths, fgout1, fgout2,
                                   drag_factor=drag_factor, 
                                   grounding_depth=grounding_depth,
                                   grounding_speed=grounding_speed)
        fgout1 = fgout2
    return debris_paths


def get_debris_xy(t, debris_paths, dbnos):
    
    """
    Determine the location of each debris particle from the list dbnos at time t,
    assuming debris_paths[dbno] has a row corresponding to this time.
    This assumes debris_paths contains the full debris paths as computed by
    make_debris_paths.
    """
    
    xd = []
    yd = []
    for dbno in dbnos:
        db = debris_paths[dbno]
        try:
            j = where(abs(db[:,0]-t) < 1e-6)[0].max()
        except:
            print('Did not find path for dbno=%i at t = %.3f' % (dbno,t))
            j = -1
        if j > -1:
            xd.append(db[j,1])
            yd.append(db[j,2])
            
    return xd,yd
            


# ==================================
# debris pairs: These functions work on pairs of debris particles and constrain
# their motion so that the distance between them remains constant.
# Useful for tracking debris with some length, and to better visualize rotation
# of particles in a velocity field.

# These functions have not yet been updated to include grounding_depth
# or drag_factor.


def move_debris_pairs(dbnosA, dbnosB, debris_paths, fgout1, fgout2):
    """
    The lists dbnosA and dbnosB should be the same length and the two particles
    dbnosA[k] and dbnosB[k] should be constrained to maintain constant distance 
    between them.
    """
    
    from clawpack.geoclaw.util import haversine 

    # first move each particle by the unconstrained algorithm:
    dbnosAB = list(dbnosA) + list(dbnosB)
    move_debris(dbnosAB, debris_paths, fgout1, fgout2)
        
    # pdb; pdb.set_trace()
    
    # constrain motion so that adjacent particles remain const distance apart:
    for k in range(len(dbnosA)):
        dbnoA = dbnosA[k]
        dbnoB = dbnosB[k]
        
        # previous positions before move_debris performed above:
        tA_old,xdA_old,ydA_old,udA_old,vdA_old = debris_paths[dbnoA][-2,:]
        tB_old,xdB_old,ydB_old,udB_old,vdB_old = debris_paths[dbnoB][-2,:]
        dist_old = haversine(xdB_old, ydB_old, xdA_old, ydA_old)

        # new positions computed by move_debris:
        tA_new,xdA_new,ydA_new,udA_new,vdA_new = debris_paths[dbnoA][-1,:]
        tB_new,xdB_new,ydB_new,udB_new,vdB_new = debris_paths[dbnoB][-1,:]
        dist_new = haversine(xdB_new, ydB_new, xdA_new, ydA_new)
        
        # now adjust so that new dist agrees with dist_old
        # keep midpoint xdm,ydm fixed and adjust distance to each end:

        ratio = dist_old/dist_new

        xdm = 0.5*(xdA_new + xdB_new)
        ydm = 0.5*(ydA_new + ydB_new)
        dx_new = xdB_new - xdA_new
        dy_new = ydB_new - ydA_new
        xdA_new = xdm - ratio*0.5*dx_new
        ydA_new = ydm - ratio*0.5*dy_new
        xdB_new = xdm + ratio*0.5*dx_new
        ydB_new = ydm + ratio*0.5*dy_new
        dist_adjusted = haversine(xdB_new, ydB_new, xdA_new, ydA_new)
        #print('+++ dbnoA=%i: distances %.1f, %.1f, %.1f' \
        #        % (dbnoA, dist_old, dist_new, dist_adjusted))
        
        # Should reset ud, vd also!?
        
        debris_paths[dbnoA][-1,:] = tA_new,xdA_new,ydA_new,udA_new,vdA_new
        debris_paths[dbnoB][-1,:] = tB_new,xdB_new,ydB_new,udB_new,vdB_new
        
    return debris_paths


def make_debris_paths_pairs(fgno, fgframes, plotdata, debris_paths,
                            dbnosA, dbnosB):
    """
    dbnosA,dbnosB are equal-length lists of debris particle labels (integers)
    to operate on, specifying ed points of long debris particles.
    """
    
    fgout1 = read_fgout_frame(fgno, fgframes[0], plotdata)
    
    for fgframe in fgframes[1:]:
        print('Trying to read fgno=%i, fgframe=%i' % (fgno,fgframe))
        try:
            fgout2 = read_fgout_frame(fgno, fgframe, plotdata)
        except:
            print('Could not read file, exiting loop')
            break
        debris_paths = move_debris_pairs(dbnosA, dbnosB, debris_paths, 
                                         fgout1, fgout2)
        fgout1 = fgout2
    return debris_paths


def plot_debris_pairs(t, debris_paths, dbnosA, dbnosB, ax,
                      color='k', linewidth=2):
    """
    Plot the location of each debris particle pair connected by a line,
    from the lists dbnosA, dbnosB at time t,
    assuming debris_paths[dbno] has a row corresponding to this time.
    This assumes debris_paths contains the full debris paths as computed by
    make_debris_paths.
    """
    
    for k in range(len(dbnosA)):
        dbnoA = dbnosA[k]
        dbnoB = dbnosB[k]
        dbA = debris_paths[dbnoA]
        dbB = debris_paths[dbnoB]
        try:
            j = where(abs(dbA[:,0]-t) < 1e-6)[0].max()
        except:
            print('Did not find path for dbno=%i at t = %.3f' % (dbno,t))
            j = -1
        if j > -1:
            xA = dbA[j,1]
            yA = dbA[j,2]
            xB = dbB[j,1]
            yB = dbB[j,2]
            ax.plot([xA,xB], [yA,yB], color=color, linewidth=linewidth)
            
