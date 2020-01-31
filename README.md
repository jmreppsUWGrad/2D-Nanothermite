# 2D Heat Conduction Code

This repository contains the Python code to solve the 2D Heat conduction equations with a new model for nano-thermite combustion.

## Solver Details
- Planar or axisymmetric geometries
- Vertex-centred, finite volume method, Explicit time advancement
- 2nd order central differences for diffusion flux; harmonic or linear interpolation at control surfaces
- Solve heat conduction equations (Heat model) or nano-thermite model (Species model)
- Customizable specific heat capacity based on reaction progress (Arrhenius source term), temperature or a constant
- Customizable thermal conductivity models and calculation methods
- Boundary Conditions: Constant, Flux or convective
- Run from command prompt; must be in parallel with equal nodes in x OR y in each process
- Can restart a simulation using variable data from previous run

## Heat Model
- 2D Heat equation
- Uniform heat generation or exponential (Arrhenius with reaction progress variable) source term options

## Species Model (In Development)
- 2D Heat equation
- Two species tracker of solid material with a fraction of it converted to a gas phase
- Exponential (Arrhenius) source term for heat generation and species conversion rate
- Darcy's law implemented in convective terms with permeability calculated via Carmen-Kozeny

## Run code from terminal:
mpiexec -n [proc] python main.py [input file name] [Output directory]

where:

[proc]-number of processors used

[input file name]-name and relative path of input file including extension in name (.txt files have been tested)

[Output directory]-relative path to directory to save data files to; will create if non-existent

## Post-processing data
### Post.py
- outputs Temperature, reaction progress, species and pressure contours
- (Optional) Reaction rate contours, Darcy velocity contours or 1D plots of afformentioned variables along centre

python Post.py [input file name]

where:

[input file name]-relative path to post-processing input file

### Post_timeEvolv.py
- output temperature, pressure or reaction progress evolution in time at one location using the same input file
