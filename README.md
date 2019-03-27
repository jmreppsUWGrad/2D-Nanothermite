# 2D Heat Conduction Code

This repository contains the Python code to solve the 2D Heat conduction equations.

# Current state:
-solve 2D Planar heat conduction with uniform heat generation (Explicit)

-solve Axisymmetric heat conduction with uniform heat generation (Explicit)

-Combustion source term from Kim (Explicit)

-can be run from command prompt

-can restart a simulation using Temperature data from previous run

-post-processing script outputs Temperature, reaction progress and reaction rate contours

# Run code from cmd:
cd [directory]

python main.py [input file name] [Output directory]

where:

[input file name]-name of input file including extension in name (.txt files have been tested); based on current directory

[Output directory]-directory to save data files to; based on current directory; will create if non-existent

# Post-processing data:
python Post-processing.py [Output directory] [1D graphs]

where:

[Output directory]-directory where data files are stored
[1D graphs] indicates whether 1D graphs should be output (1 or 0); default is 0