Run the GeoClaw code in this directory to simulate the tsunami generated
by a submarine slump.

Clawpack v5.10.0 (or later) is required, since the dispersive SGN equation
is used rather than shallow water equations.  See
    https://www.clawpack.org/installing.html
    https://www.clawpack.org/bouss2d.html

For general information on running GeoClaw, see the documentation at
    https://www.clawpack.org/contents.html#using-the-fortran-codes
    https://www.clawpack.org/geoclaw.html

Before running this code, topo files and dtopo files must be created.
See the README.txt files in the topo and dtopo directories.

The dtopo file to be used is specified in setrun.py, in the line

    dtopo_data.dtopofiles = [[3, dtopodir + '/slump_10m.dtt3']]

and the code must be run separately for each amplitude, using the dtopo files
    slump_05m.dtt3, slump_10m.dtt3, slump_15m.dtt3
    
In order to run the post-processing scripts to transport particles or
create plots of time frames or gauges, redirect the output to an appropriate
output directory, e.g. for the 10 m slump:

    make .output OUTDIR=_output_10m
    
or you can rename _output after running the code.

--------------------------------------------------------------------
Make gauge plots:

After running the code, the gauge plots can be created by running

    python plot_gauges.py
    
This is set up to loop over 1 or more output directories, making gauge plots
for each.  (Adjust the loop over mstr in the main program to adjust.)

--------------------------------------------------------------------
Make frame plots:

The frame plots showing inundation and particle positions at various times
are created by running:

    python plot_frames.py

Adjust mstr in this script to adjust which output directory is being used.

--------------------------------------------------------------------
Make animation:

An animation can be made via:

    python make_animation.py

Adjust mstr in this script to adjust which output directory is being used.
Creating the animation requires that ffmpeg is installed.

