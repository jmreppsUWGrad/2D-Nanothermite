# -*- coding: utf-8 -*-
"""
2D Heat Conduction solver

Features:
    -Customizable boundary conditions (along each side) including distributions
    -Explicit and implicit solver
    -Built in Fourrier number stability check (not incl convective check)

Desired:
    -Combustion reaction modelling (own class)
    -meshing tool with biasing [DONE; need to test with each solver]
        [BC function NOT reflecting biased mesh]
    -Option to apply either flux BC to corners; -2 index to reflect this?
    -Radiation BC [Done]
    -Cantera use via Source_Comb class
    -Read an input file; BC settings need re-coding

NANOTHERMITE TESTING:
    settings['Length']                  = 10**(-3)
    settings['Width']                   = 6.0*10**(-3)
    settings['Nodes_x']                 = 101
    settings['Nodes_y']                 = 601
    settings['k']                       = 10
    settings['Cp']                      = 800
    settings['rho']                     = 8000
    BCs['bc_left']                      = ['F',0,(0,-1)]
    BCs['bc_right']                     = ['C',(30,300),(0,-1)]
    BCs['bc_south']                     = ['F',0,(0,-1)]
    BCs['bc_north']                     = ['F',4*10**8,(1,10-settings['Nodes_x']),'C',(30,300),(10,-1)]
    
fig2=pyplot.figure(figsize=(7,7))
pyplot.plot(domain.Y[:,1]*1000, domain.T[:,1],marker='x')
pyplot.xlabel('$y$ (mm)')
pyplot.ylabel('T (K)')
pyplot.title('Temperature distribution at 2nd x')
pyplot.xlim(5,6);

fig4=pyplot.figure(figsize=(7, 7))
pyplot.contourf(domain.X*1000, domain.Y*1000, domain.T, alpha=0.5, cmap=cm.viridis)  
pyplot.colorbar()
pyplot.xlabel('$x$ (mm)')
pyplot.ylabel('$y$ (mm)')
pyplot.title('Temperature distribution')
pyplot.xlim(0,4)
pyplot.ylim(5,6);

@author: Joseph
"""

##########################################################################
# ----------------------------------Libraries and classes
##########################################################################
import numpy as np
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
settings['Length']                  = 1.0
settings['Width']                   = 1.0
settings['Nodes_x']                 = 31
settings['Nodes_y']                 = 31
settings['k']                       = 10 #0.026384465709828872
settings['Cp']                      = 800 # 714.602
settings['rho']                     = 8000 #1.2

# Source terms
settings['Source_Uniform']          = 'None'
settings['Source_Kim']              = 'None'

# Meshing details
"""
Biasing options:
    -'OneWayUp'   for linearly increasing element sizes with increasing x/y
    -'OneWayDown' for linearly decreasing element sizes with increasing x/y
    -'TwoWayEnd'  for linearly increasing sizes till middle, then decrease again
    -'TwoWayMid'  for linearly decreasing sizes till middle, then increase again
    -size         is the smallest element size based on above selection
"""
settings['bias_type_x']             = 'None'
settings['bias_size_x']             = 0.01 # Smallest element size
settings['bias_type_y']             = 'None'
settings['bias_size_y']             = 0.01 # Smallest element size

# Boundary conditions
"""
Boundary condition options:
    -'T' for constant temperature boundary (specify T)
    -'F' for constant flux boundary (specify q")
    -'C' for convective boundary (specify h and Tinf)
    -format: ['T', 200, (0,-1), ...]
        First index: type of BC
        Second index: Numbers associated with BC (can be an array)
        Third index: Node index range BC is valid for (second number must be negative)
        this pattern repeats...
    -radiation options: None or [emissivity, surrounding_Temp]
"""
#['C',(30,300),(0,-1)]
#['F',4*10**8,(1,-299),'C',(10,300),(2,-2)]
BCs['bc_left']                      = ['T',600,(0,-1)]
BCs['bc_left_rad']                  = 'None'
# numpy.linspace(400, 900, settings['Nodes_y'])
BCs['bc_right']                     = ['T',300,(0,-1)]
BCs['bc_right_rad']                 = 'None'
# numpy.linspace(400, 900, settings['Nodes_y'])
BCs['bc_south']                     = ['T',600,(0,-1)]
BCs['bc_south_rad']                 = 'None'
# numpy.linspace(400, 900, settings['Nodes_x'])
BCs['bc_north']                     = ['T',300,(0,-1)]
BCs['bc_north_rad']                 = 'None'
# numpy.linspace(400, 900, settings['Nodes_x'])

# Time advancement
settings['Fo']                      = 0.1
settings['dt']                      = 'None' # Time step
settings['total_time_steps']        = 1000
settings['total_time']              = 'None'
settings['Time_Scheme']             = 'Explicit' # Explicit or Implicit
settings['Convergence']             = 0.0001 # implicit solver only
settings['Max_iterations']          = 100 #    implicit solver only

print('######################################################')
print('#             2D Heat Conduction Solver              #')
print('#              Created by J. Mark Epps               #')
print('#          Part of Masters Thesis at UW 2018-2020    #')
print('######################################################\n')

##########################################################################
# -------------------------------------Read input file
##########################################################################
print 'Reading input file...'
settings={}
BCs={}
#fin=FileClasses.FileIn('Input_File', 0)
fin=FileClasses.FileIn('Input_File_nt', 0)
fin.Read_Input(settings, BCs)
os.chdir(settings['Output_directory'])

print '################################'

# Initial conditions from previous run/already in memory
#Use_inital_values                   = False


print 'Initializing geometry package...'
#domain=OneDimLine(L,Nx)
domain=TwoDimPlanar(settings, 'Solid')
domain.mesh()
print '################################'

##########################################################################
# -------------------------------------Initialize solver and domain
##########################################################################

print 'Initializing solver package...'
solver=Solvers.TwoDimPlanarSolve(domain, settings, BCs, 'Solid')
print '################################'

print 'Initializing domain...'
domain.T[:,:]=300
#T[:2,:]=600
print '################################'
##########################################################################
# -------------------------------------File setups
##########################################################################
#print 'Initializing files...'
#datTime=str(datetime.date(datetime.now()))+'_'+'{:%H%M}'.format(datetime.time(datetime.now()))
#isBinFile=False
#
##output_file=FileClasses.FileOut('Output_'+datTime, isBinFile)
#input_file=FileClasses.FileOut('Input_'+datTime, isBinFile)
#
## Write headers to files
#input_file.header_cond('INPUT')
##output_file.header('OUTPUT')
#
## Write input file with settings
#input_file.input_writer_cond(settings, BCs, domain.T)
#input_file.close()
#print '################################\n'

print 'Saving data to numpy array files...'
np.save('T_'+'0.000000', domain.T, False)
np.save('eta_'+'0.000000', domain.eta, False)
np.save('X', domain.X, False)
np.save('Y', domain.Y, False)

##########################################################################
# -------------------------------------Solve
##########################################################################
t,nt=0,0
output_data_t,output_data_nt=0,0
if settings['total_time_steps']=='None':
    settings['total_time_steps']=settings['total_time']*10**9
    output_data_t=settings['total_time']/settings['Number_Data_Output']
elif settings['total_time']=='None':
    settings['total_time']=settings['total_time_steps']*10**9
    output_data_nt=int(settings['total_time_steps']/settings['Number_Data_Output'])

BCs_changed=False
print 'Solving:'
while nt<settings['total_time_steps'] and t<settings['total_time']:
#for nt in range(settings['total_time_steps']):
    err,dt=solver.Advance_Soln_Cond(nt, t)
    t+=dt
    nt+=1
#    print 'Time step %i, Step size=%.7f, Time elapsed=%f;'%(nt+1,dt, t)
    if err==1:
        print '#################### Solver aborted #######################'
        break
    
    # Output data to numpy files
    if output_data_nt!=0 and nt%output_data_nt==0:
        print 'Saving data to numpy array files...'
        np.save('T_'+'{:f}'.format(t), domain.T, False)
        np.save('eta_'+'{:f}'.format(t), domain.eta, False)
        
    # Change boundary conditions
    if np.amax(domain.eta)>=1.0 and not BCs_changed:
        solver.BCs['bc_north']=['C',(5,300),(0,-1)]
        BCs_changed=True
#        break
    
#output_file.close()

##########################################################################
# ------------------------------------Post-processing
##########################################################################
T, eta=domain.T, domain.eta
#fig2=pyplot.figure(figsize=(7,7))
#pyplot.plot(domain.Y[:,1]*1000, domain.T[:,1],marker='x')
#pyplot.xlabel('$y$ (mm)')
#pyplot.ylabel('T (K)')
#pyplot.title('Temperature distribution at 2nd x')
#pyplot.xlim(5,6);

# Nano thermite testing figures
fig4=pyplot.figure(figsize=(7, 7))
pyplot.contourf(domain.X*1000, domain.Y*1000, T, alpha=0.5, cmap=cm.viridis)  
pyplot.colorbar()
pyplot.xlabel('$x$ (mm)')
pyplot.ylabel('$y$ (mm)')
pyplot.title('Temperature distribution, t=%.7f'%t)
#pyplot.xlim(0,0.4)
#pyplot.ylim(5,6);

fig4=pyplot.figure(figsize=(7, 7))
pyplot.contourf(domain.X*1000, domain.Y*1000, eta, alpha=0.5, cmap=cm.viridis)  
pyplot.colorbar()
pyplot.xlabel('$x$ (mm)')
pyplot.ylabel('$y$ (mm)')
pyplot.title('Reaction progress, t=%.7f'%t)
#pyplot.xlim(0,0.4)
#pyplot.ylim(5,6);
pyplot.close(fig4)

# 2D plot
#fig=pyplot.figure(figsize=(7, 7))
#ax = fig.gca(projection='3d')
#ax.plot_surface(domain.X, domain.Y, T, rstride=1, cstride=1, cmap=cm.viridis,linewidth=0, antialiased=True)
##ax.set_xlim(0,0.001)
##ax.set_ylim(0.005,0.006)
#ax.set_zlim(300, 700)
#ax.set_xlabel('$x$ (m)')
#ax.set_ylabel('$y$ (m)')
#ax.set_zlabel('T (K)');
#fig.savefig(datTime+'_2DPlot.png',dpi=300)

#fig2=pyplot.figure(figsize=(7,7))
#pyplot.plot(numpy.zeros(len(domain.Y[:,1])), domain.Y[:,1], marker='x')
#pyplot.plot(domain.X[1,:], numpy.zeros(len(domain.X[1,:])), marker='x')
#pyplot.ylabel('Space')
#pyplot.title('Discretizations');

# 1D Plot
#fig2=pyplot.figure(figsize=(7,7))
#pyplot.plot(domain.Y[:,1], T[:,1],marker='x')
#pyplot.xlabel('$y$ (m)')
#pyplot.ylabel('T (K)')
#pyplot.title('Temperature distribution at 2nd x')
##pyplot.xlim(5,6);
#fig2.savefig(datTime+'_Plot2.png',dpi=300)

# Temperature contour
#fig4=pyplot.figure(figsize=(7, 7))
#pyplot.contourf(domain.X, domain.Y, T, alpha=0.5, cmap=cm.viridis)  
#pyplot.colorbar()
#pyplot.xlabel('$x$ (m)')
#pyplot.ylabel('$y$ (m)')
#pyplot.title('Temperature distribution');
#fig4.savefig(datTime+'_Temp.png',dpi=300)

print('Solver has finished its run')