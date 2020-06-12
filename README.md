# 2D Heat Conduction Code

This repository contains the Python (2.7) code to solve the 2D Heat conduction equations with a new model for nano-thermite combustion.

## Solver Details
- Planar or axisymmetric geometries
- Vertex-centred, finite volume method, Explicit time advancement
- 2nd order central differences for diffusion fluxes; 1st order harmonic or linear interpolation at control surfaces
- Solve heat conduction equations (Heat model) or nano-thermite model (Species model)
- Customizable specific heat capacity based on reaction progress (Arrhenius source term), temperature or a constant
- Customizable thermal conductivity models and calculation methods
- Run from command prompt, parallel code (MPI)
- Can restart a simulation using variable data from previous run

## Heat Model
- 2D Heat equation
- Uniform heat generation or exponential (Arrhenius form with reaction progress variable) source term options
- Boundary Conditions: Constant, Flux or convective

## Species Model
- 2D Heat equation
- Two phase tracker of solid pellet with a fraction of it converted to a gas phase
- Exponential (Arrhenius) source term for heat generation and species conversion rate
- Darcy's law implemented in convective terms with permeability calculated via Kozeny-Carmen
- Boundary Conditions: Heat Model + Dirichlet pressure (mass allowed to leave to maintain pressure)

## Run code from terminal:
mpiexec -n [proc] python main.py [input file name] [Output directory]

where:

[proc]-number of processors used

[input file name]-name and relative path of input file including extension in name (.txt files have been tested)

[Output directory]-relative path to directory to save data files to; will create if non-existent

## Post-processing data
### Post.py
- Calculate characteristic values and non-dimensional quantities, percentage mass loss to ambient (via boundaries for Species model)
- Specify directory for data (*Post_Input.txt* file)
- Customize temperature colorbar bounds and bin sizes (*Post_Input.txt* file)
- (Optional contours) Temperature, reaction progress, species, pressure, reaction rate, Darcy velocities
- (Optional) 1D plots of afformentioned variables along centre of geometry

python Post.py [input file name]

where:

[input file name]-relative path to post-processing input file

### Post_timeEvolv.py
- Temperature and one of pressure, densities or reaction progress evolution in time at one location
- Uses same input file and command call as Post.py
