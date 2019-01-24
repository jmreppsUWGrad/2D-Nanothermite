# 2D Heat Conduction Code

This repository contains the Python code to solve the 2D Heat conduction equations.

# Current state:
-solve 2D heat conduction with no source terms (Explicit or implicit)

-Combustion source term from Kim (Explicit only)

-can be run from command prompt

-post-processing script outputs Temperature, reaction progress and reaction rate contours

# Run code from cmd:
cd [directory]

python main.py [input file name] [Output directory]

where:

[input file name]-name of input file; must be txt file; do NOT include extension in name

[Output directory]-directory to save data files to; based on current directory; will create if non-existent

# Post-processing data:
python Post-processing.py [Output directory]

where:

[Output directory]-directory where data files are stored