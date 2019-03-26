# 2D Heat Conduction Code

This repository contains the Python code to solve the 2D Heat conduction equations.

# Current state:
-solve 2D planar heat conduction with no source terms (Explicit)

-solve 2D axisymmetric heat conduction with no source terms (Explicit)

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
python Post-processing.py [Output directory]

where:

[Output directory]-directory where data files are stored