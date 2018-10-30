# -*- coding: utf-8 -*-
"""
2D Heat Conduction solver

Features:
    -Customizable boundary conditions (along each side)
    -Explicit solver (Implicit in the works)
    -

Boundary condition options:
    -'zero_grad' will impose 0 normal gradient of that variable
    
Features to include (across all classes):
    -CoolProp library for material properties (track down needed functions)
    -Fix biasing meshing tools (this script and GeomClasses)
        ->Figure out biasing wrt dx and dy array sizes and mesh griding those (GeomClasses)
    -File reader for settings

Desired:
    -Cantera for adding combustion reactions
    -Implicit solver
    

@author: Joseph
"""

##########################################################################
# ----------------------------------Libraries and classes
##########################################################################
import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime
import os
#import CoolProp.CoolProp as CP

#from GeomClasses import OneDimLine as OneDimLine
from GeomClasses import TwoDimPlanar as TwoDimPlanar
#import MatClasses as Mat
import SolverClasses as Solvers
import FileClasses

##########################################################################
# ------------------------------ Geometry, Domain and BCs Setup
#    Reference directions:
#    left-smallest x coordinate
#    right-largest x value
#    north-largest y coordinate
#    south-smallest y coordinate
##########################################################################
settings={} # Dictionary of problem settings
BCs={} # Dictionary of boundary conditions
# Geometry details
settings['Length']                  = 4.0
settings['Width']                   = 4.0
settings['Nodes_x']                 = 13
settings['Nodes_y']                 = 13
settings['k']                       = 10
settings['Cp']                      = 800
settings['rho']                     = 8000

# Meshing details
settings['bias_type_x']             = None
settings['bias_size_x']             = 0.003 # Smallest element size (IN PROGRESS)
settings['bias_type_y']             = None
settings['bias_size_y']             = 10**(-6) # Smallest element size (IN PROGRESS)

# Boundary conditions
#['C',(10,300),(1,-2)]
#['F',4*10**8,(1,-299),'C',(10,300),(2,-2)]
BCs['bc_left']                      = ['T',600,(0,-1)]
#BCs['bc_left_flux']                 = None
#BCs['bc_left_conv']                 = None
# numpy.linspace(400, 900, settings['Nodes_y'])
BCs['bc_right']                     = ['T',300,(0,-1)]
#BCs['bc_right_flux']                = None
#BCs['bc_right_conv']                = None
# numpy.linspace(400, 900, settings['Nodes_y'])
BCs['bc_south']                     = ['T',600,(0,-1)]
#BCs['bc_south_flux']                = None
#BCs['bc_south_conv']                = None
# numpy.linspace(400, 900, settings['Nodes_x'])
BCs['bc_north']                     = ['T',300,(0,-1)]
#BCs['bc_north_flux']                = None
#BCs['bc_north_conv']                = None
# numpy.linspace(400, 900, settings['Nodes_x'])

# Time advancement
settings['Fo']                      = 0.2
settings['total_time_steps']        = 10
settings['Time_Scheme']             = 'Implicit'
settings['Convergence']             = 0.001
settings['Max_iterations']          = 50


print 'Initializing geometry package...'
#domain=OneDimLine(L,Nx)
domain=TwoDimPlanar(settings, 'Solid')
domain.mesh()

##########################################################################
# -------------------------------------Initialize solver and domain
##########################################################################

print 'Initializing solver package...'
solver=Solvers.TwoDimPlanarSolve(domain, settings, BCs)

print 'Initializing domain...'
domain.T[:,:]=300
#T[:2,:]=600

#domain.T[1:-1,1:-1]=domain.p[1:-1,1:-1]/domain.rho[1:-1,1:-1]/domain.R

##########################################################################
# -------------------------------------File setups
##########################################################################
print 'Initializing files...'
os.chdir('Tests')
datTime=str(datetime.date(datetime.now()))+'_'+'{:%H%M}'.format(datetime.time(datetime.now()))
isBinFile=False

#output_file=FileClasses.FileOut('Output_'+datTime, isBinFile)
input_file=FileClasses.FileOut('Input_'+datTime, isBinFile)

# Write headers to files
input_file.header_cond('INPUT')
#output_file.header('OUTPUT')

# Write input file with settings
input_file.input_writer_cond(settings, BCs, domain.T)
input_file.close()

##########################################################################
# -------------------------------------Solve
##########################################################################
print('######################################################')
print('#             2D Heat Conduction Solver              #')
print('#              Created by J. Mark Epps               #')
print('#          Part of Masters Thesis at UW 2018-2020    #')
print('######################################################\n')

print 'Solving:'
for nt in range(settings['total_time_steps']):
    print 'Time step %i of %i'%(nt+1, settings['total_time_steps'])
    err=solver.Advance_Soln_Cond()
    if err==1:
        print '#################### Solver aborted #######################'
        break

#output_file.close()

##########################################################################
# ------------------------------------Post-processing
##########################################################################
# 2D plot
#fig=pyplot.figure(figsize=(7, 7))
#ax = fig.gca(projection='3d')
#ax.plot_surface(domain.X, domain.Y, T, rstride=1, cstride=1, cmap=cm.viridis,linewidth=0, antialiased=True)
##ax.set_xlim(0,0.001)
##ax.set_ylim(0.005,0.006)
#ax.set_zlim(300, BCs['bc_south_T'])
#ax.set_xlabel('$x$ (m)')
#ax.set_ylabel('$y$ (m)')
#ax.set_zlabel('T (K)');
#fig.savefig('Plot1.png',dpi=300)

# 1D Plot
#fig2=pyplot.figure(figsize=(7,7))
#pyplot.plot(domain.Y[:,1]*1000, domain.T[:,1],marker='x')
#pyplot.xlabel('$y$ (mm)')
#pyplot.ylabel('T (K)')
#pyplot.title('Temperature distribution at 2nd x')
#pyplot.xlim(5,6);
#fig2.savefig('Plot2.png',dpi=300)

# Temperature contour
fig4=pyplot.figure(figsize=(7, 7))
pyplot.contourf(domain.X, domain.Y, domain.T, alpha=0.5, cmap=cm.viridis)  
pyplot.colorbar()
pyplot.xlabel('$x$ (m)')
pyplot.ylabel('$y$ (m)')
pyplot.title('Temperature distribution');
#fig4.savefig(datTime+'_Temp.png',dpi=300)

print('Solver has finished its run')