# 2D Heat Conduction Code

This repository contains the Python code to solve the 2D Heat conduction equations.

# Features:
-Vertex-centred, finite volume method, Explicit

-2nd order central differences for diffusion flux, harmonic or linear interpolation at control surfaces

-Planar or axisymmetric geometries

# Current state:
-Uniform heat generation source term

-Combustion source term from Kim

-two species model with Darcy's law, species tracker (gas and solid only); TESTING PHASE

-can be run from command prompt; must be in parallel with equal nodes in x OR y in each process

-can restart a simulation using variable data from previous run

-post-processing script outputs Temperature, reaction progress, reaction rate, species and pressure contours

# Run code from cmd:
cd [directory]

mpiexec -n [proc] python main.py [input file name] [Output directory]

where:

[proc]-number of processors used

[input file name]-name of input file including extension in name (.txt files have been tested); based on current directory

[Output directory]-directory to save data files to; based on current directory; will create if non-existent

# Post-processing data:
python Post.py [Output directory] [1D graphs]

where:

[Output directory]-directory where data files are stored

[1D graphs] indicates whether 1D graphs should be output (1 or 0)