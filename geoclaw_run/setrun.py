"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

This setrun.py assumes GeoClaw from Clawpack v5.10.0 (or later) is being
used, since it solves the dispersive SGN equation rather than shallow
water equations.  See  https://www.clawpack.org/bouss2d.html

"""


import os
import numpy as np
from clawpack.amrclaw.data import FlagRegion
from clawpack.geoclaw import fgmax_tools
from clawpack.geoclaw import fgout_tools

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")


rundir = os.getcwd()
print('rundir = %s' % rundir)


# set topodir and dtopodir to directory where topo and dtopo files are found:
topodir = '../topo'
dtopodir = '../dtopo'

#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)


    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------
    
    #probdata = rundata.new_UserData(name='probdata',fname='setprob.data')

    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------
    try:
        geo_data = rundata.geo_data
    except:
        print("*** Error, this rundata has no geo_data attribute")
        raise AttributeError("Missing geo_data attribute")

    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    #------------------------------------------------------------------
    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # Set single grid parameters first.
    # See below for AMR parameters.


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:

    # shift so that cell centers on finest grid align with DEM:

    sec16 = 1./(6*3600.)  # one sixth arcsecond

    #Get etopo1 of -180,-66,-48,32  
    clawdata.lower[0] = -156.4        # west longitude
    clawdata.upper[0] = -156.0        # east longitude
    clawdata.lower[1] = 20.4          # south latitude
    clawdata.upper[1] = 20.7          # north latitude

    clawdata.num_cells[0] = 72
    clawdata.num_cells[1] = 54


    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 5  # Required for Boussinesq solvers

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 3

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2

    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0


    # Restart from checkpoint file of a previous run?
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False         # True to restart from prior results
    clawdata.restart_file = ''


    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 1

    if clawdata.output_style==1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.num_output_times = 40
        clawdata.tfinal = 20*60
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        clawdata.output_times = [0, 20*60]

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 1
        clawdata.output_t0 = True
        

    clawdata.output_format = 'binary'

    clawdata.output_q_components = 'all'   # need all
    clawdata.output_aux_components = 'none'  # eta=h+B is in q
    clawdata.output_aux_onlyonce = False    # output aux arrays each frame



    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 1

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.2

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.75
    clawdata.cfl_max = 1.0

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 5000


    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'
    
    # For unsplit method, transverse_waves can be 
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2

    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    clawdata.limiter = ['mc', 'mc', 'mc']

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms
    
    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 'godunov'


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'

    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'



    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    # negative checkpoint_style means alternate between aaaaa and bbbbb files
    # so that at most 2 checkpoint files exist at any time, useful when
    # doing frequent checkpoints of large problems.

    clawdata.checkpt_style = 0

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif clawdata.checkpt_style == 1:
        # Checkpoint only at tfinal.
        pass

    elif abs(clawdata.checkpt_style) == 2:
        # Specify a list of checkpoint times.  
        clawdata.checkpt_times = 3600.*np.arange(1,16,1)

    elif abs(clawdata.checkpt_style) == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5


    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 4

    # List of refinement ratios at each level (length at least mxnest-1)

    # dx = dy = 20", 2", 1/3", 1/9"
    amrdata.refinement_ratios_x = [5,8,6,3]
    amrdata.refinement_ratios_y = [5,8,6,3]
    amrdata.refinement_ratios_t = [5,8,6,3]
    

    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    amrdata.aux_type = ['center','capacity','yleft']


    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag2refine = True

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.700000

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0  

       
    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2
    geo_data.earth_radius = 6367.5e3

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 1.e-3
    geo_data.friction_forcing = True
    geo_data.manning_coefficient =.025
    geo_data.friction_depth = 1e6

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.5

    # == settopo.data values ==
    topo_data = rundata.topo_data
    # for topography, append lines of the form
    #    [topotype, minlevel, maxlevel, t1, t2, fname]

    topofiles = topo_data.topofiles

    # 6-second topo:
    topofiles.append([3, topodir + '/hawaii6s_nuu.asc'])

    # 1/9 arcsecond topo:
    topofiles.append([3, topodir + '/nuu19s.asc'])
        

    # == setdtopo.data values ==
    dtopo_data = rundata.dtopo_data
    # for moving topography, append lines of the form :   (<= 1 allowed for now!)
    #   [topotype, fname]

    dtopo_data.dtopofiles = [[3, dtopodir + '/slump_05m.dtt3']]

    dtopo_data.dt_max_dtopo = 2.


    # == setqinit.data values ==
    rundata.qinit_data.qinit_type = 0
    rundata.qinit_data.qinitfiles = []
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]
    #rundata.qinit_data.qinitfiles.append(['hump.xyz'])

    rundata.qinit_data.variable_eta_init = True  # for Nu'u Pond
    #rundata.qinit_data.variable_eta_init = False  


    # ---------------
    # Regions:
    # ---------------
    rundata.regiondata.regions = []
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]

    flagregions = rundata.flagregiondata.flagregions  # initialized to []
    
    # dx = dy = 2deg, 24', 4', 1', 30", 6", 2", 1/3"":
    # refine would be [5,6,4,2,5,3,6]
    
    # Computational domain Variable Region - 1degree to 4min to 1min:
    # Level 3 below is 4 min, level 4 is 1 min
    # Note that this is a rectangle specified in the new way
    # (other regions below will force/allow more refinement)
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_domain'
    flagregion.minlevel = 1
    flagregion.maxlevel = 2
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [clawdata.lower[0]-0.1,
                                 clawdata.upper[0]+0.1,
                                 clawdata.lower[1]-0.1,
                                 clawdata.upper[1]+0.1]
    flagregions.append(flagregion)

    # Source Region 
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_dtopo'
    flagregion.minlevel = 3
    flagregion.maxlevel = 3
    flagregion.t1 = 0.
    flagregion.t2 = 60.
    flagregion.spatial_region_type = 1  # Rectangle
    source_region = [-156.24, -156.16, 20.55, 20.63]
    flagregion.spatial_region = source_region
    flagregions.append(flagregion)

    # Nearshore Region 
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_dtopo'
    flagregion.minlevel = 2
    flagregion.maxlevel = 3
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    source_region = [-156.22, -156.17, 20.58, 20.63]
    flagregion.spatial_region = source_region
    flagregions.append(flagregion)


    if 1:
        # Maui Site - require  1/3" grids:
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_Maui_Site_13sec'
        flagregion.minlevel = 3
        flagregion.maxlevel = 3
        flagregion.t1 = 0.
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = [-156.189, -156.17, 20.618, 20.63]
        flagregions.append(flagregion)


    if 1:
        # Maui Site - require  1/9" grids:
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_Maui_Site_19sec'
        flagregion.minlevel = 4
        flagregion.maxlevel = 4
        flagregion.t1 = 0.
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = [-156.183, -156.175, 20.622, 20.63]
        flagregions.append(flagregion)

    # To use Boussinesq solver, add bouss_data parameters here
    # Also make sure to use the correct Makefile pointing to bouss version
    # and set clawdata.num_eqn = 5
    from clawpack.geoclaw.data import BoussData
    rundata.add_data(BoussData(),'bouss_data')

    rundata.bouss_data.bouss_equations = 2    # 0=SWE, 1=MS, 2=SGN
    rundata.bouss_data.bouss_min_level = 1    # coarsest level to apply bouss
    rundata.bouss_data.bouss_max_level = 3   # finest level to apply bouss
    rundata.bouss_data.bouss_min_depth = 5.  # depth to switch to SWE
    rundata.bouss_data.bouss_solver = 3       # 1=GMRES, 2=Pardiso, 3=PETSc
    rundata.bouss_data.bouss_tstart = 0.      # time to switch from SWE


    # ---------------
    # Gauges:
    # ---------------
    rundata.gaugedata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    ###  GAUGES:  Pool(1), beach-pool(2),offshore-bay(3), offshore(4)
    rundata.gaugedata.gauges.append([1,-156.1766, 20.6272, 0, 1.e10])
    rundata.gaugedata.gauges.append([2,-156.1780, 20.6264, 0, 1.e10])
    rundata.gaugedata.gauges.append([3,-156.1805, 20.6235, 0, 1.e10])
    rundata.gaugedata.gauges.append([4,-156.185, 20.619, 0, 1.e10])
    
    # Gauge locations for paper:
    xg102,yg102 = -156.17713,  20.62722
    xg101,yg101 = -156.17777,  20.62622
    xg103,yg103 = -156.176917, 20.62752
    rundata.gaugedata.gauges.append([101, xg101, yg101, 0, 1.e10])
    rundata.gaugedata.gauges.append([102, xg102, yg102, 0, 1.e10])
    rundata.gaugedata.gauges.append([103, xg103, yg103, 0, 1.e10])



    # == setfixedgrids.data values ==
    #fixedgrids = rundata.fixed_grid_data.fixedgrids
    # THIS HAS BEEN DEPRECATED.  Use fgmax and/or fgout instead

    # -----------------------------
    # FGMAX GRIDS:
    # NEW STYLE STARTING IN v5.7.0
    # ------------------------------
    # set num_fgmax_val = 1 to save only max depth,
    #                     2 to also save max speed,
    #                     5 to also save max hs,hss,hmin
    rundata.fgmax_data.num_fgmax_val = 2

    fgmax_grids = rundata.fgmax_data.fgmax_grids  # empty list to start

    # Now append to this list objects of class fgmax_tools.FGmaxGrid
    # specifying any fgmax grids.

    #Set fgmax_extent: want the boundaries to be cell centers, so
    #can be an even integer plus some multiple of 1/3sec.  This will be true
    #if the decimal part is a multiple of .1, or .01, .005, or .0025 for
    #example based on how the computational domain was shifted.  This

    #fgmax_extent could be set in params.py if one is using that.
    fgmax_extent=[-156.1825,-156.175,20.625,20.628] 

    if 1:
        # Points on a uniform 2d grid:
        dx_fine = 1./(9*3600.)  # grid resolution at finest level
        fg = fgmax_tools.FGmaxGrid()
        fg.point_style = 2  # uniform rectangular x-y grid
        fg.x1 = fgmax_extent[0] #+ dx_fine/2.
        fg.x2 = fgmax_extent[1] #- dx_fine/2.
        fg.y1 = fgmax_extent[2] #+ dx_fine/2.
        fg.y2 = fgmax_extent[3] #- dx_fine/2.
        fg.dx = dx_fine
        fg.tstart_max =  2100. # when to start monitoring max values
        fg.tend_max = 1.e10                       # when to stop monitoring max values
        fg.dt_check = 5.            # target time (sec) increment between updating
                                    # max values
                                    # which levels to monitor max on
        fg.min_level_check = amrdata.amr_levels_max
        fg.arrival_tol = 1.e-2      # tolerance for flagging arrival
        fg.interp_method = 0        # 0 ==> pw const in cells, recommended
        fgmax_grids.append(fg)      # written to fgmax_grids.data

    # == fgout_grids.data values ==
    # NEW IN v5.9.0
    # Set rundata.fgout_data.fgout_grids to be a list of
    # objects of class clawpack.geoclaw.fgout_tools.FGoutGrid:

    fgout_grids = rundata.fgout_data.fgout_grids  # empty list initially


    if 1:
        # grid south Maui
        fgout_dx = 1./1200   # target resolution
        fgout_dt = 10.      # time increment (sec)
        fgout = fgout_tools.FGoutGrid()
        fgout.fgno = 2
        fgout.point_style = 2       # will specify a 2d grid of points
        fgout.output_format = 'binary32'
        fgout.x1 = -156.4  # specify edges (fgout pts will be cell centers)
        fgout.x2 = -156.0 
        fgout.y1 = 20.4
        fgout.y2 = 20.7
        fgout.nx = int(round((fgout.x2 - fgout.x1) / fgout_dx))
        fgout.ny = int(round((fgout.y2 - fgout.y1) / fgout_dx))
        fgout.tstart = 0.
        fgout.tend = 2*3600.
        fgout.nout = int(round((fgout.tend - fgout.tstart) / fgout_dt)) + 1
        fgout_grids.append(fgout)    # written to fgout_grids.data

    if 1:
        # grid around Nu'u Refuge pond with 1/9" resolution
        fgout_dx = 1./(9*3600)   # target resolution
        fgout_dt = 10.      # time increment (sec)
        fgout = fgout_tools.FGoutGrid()
        fgout.fgno = 5
        fgout.point_style = 2       # will specify a 2d grid of points
        fgout.output_format = 'binary32'
        #fgmax_extent=[-156.1825,-156.175,20.625,20.628] 
        fgout.x1 = -156.1825  # specify edges (fgout pts will be cell centers)
        fgout.x2 = -156.175
        fgout.y1 = 20.622
        fgout.y2 = 20.63
        fgout.nx = int(round((fgout.x2 - fgout.x1) / fgout_dx))
        fgout.ny = int(round((fgout.y2 - fgout.y1) / fgout_dx))
        fgout.tstart = 0.
        fgout.tend = 2*3600.
        fgout.nout = int(round((fgout.tend - fgout.tstart) / fgout_dt)) + 1
        fgout_grids.append(fgout)    # written to fgout_grids.data



    #  ----- For developers ----- 
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting
    
    return rundata

if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    from clawpack.geoclaw import kmltools
    rundata = setrun(*sys.argv[1:])
    rundata.write()
    
    kmltools.make_input_data_kmls(rundata)
