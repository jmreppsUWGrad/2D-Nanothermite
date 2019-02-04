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
#from matplotlib import pyplot, cm
#from mpl_toolkits.mplot3d import Axes3D
#from datetime import datetime
import os
import sys
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
print('######################################################')
print('#             2D Heat Conduction Solver              #')
print('#              Created by J. Mark Epps               #')
print('#          Part of Masters Thesis at UW 2018-2020    #')
print('######################################################\n')

# Get arguments to script execution
settings={}
BCs={}
Sources={}
inputargs=sys.argv
if len(inputargs)>2:
    input_file=inputargs[1]
    settings['Output_directory']=inputargs[2]
else:
    print 'Usage is: python main.py [Input file] [Output directory]\n'
    print 'where\n'
    print '[Input file] is the name of the input file with extension; must be in current directory'
    print '[Output directory] is the directory to output the data; will create relative to current directory if it does not exist'
    print '***********************************'
    sys.exit('Solver not started')
##########################################################################
# -------------------------------------Read input file
##########################################################################
print 'Reading input file...'
#fin=FileClasses.FileIn('Input_File', 0)
#fin=FileClasses.FileIn('Input_File_nt', 0)
fin=FileClasses.FileIn(input_file, 0)
fin.Read_Input(settings, Sources, BCs)
try:
    os.chdir(settings['Output_directory'])
except:
    os.makedirs(settings['Output_directory'])
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
solver=Solvers.TwoDimPlanarSolve(domain, settings, Sources, BCs, 'Solid')
print '################################'

print 'Initializing domain...'
domain.T[:,:]=300
#T[:2,:]=600
print '################################'
##########################################################################
# -------------------------------------File setups
##########################################################################
print 'Saving input file to output directory...'
#datTime=str(datetime.date(datetime.now()))+'_'+'{:%H%M}'.format(datetime.time(datetime.now()))
isBinFile=False

#output_file=FileClasses.FileOut('Output_'+datTime, isBinFile)
input_file=FileClasses.FileOut('Input_file', isBinFile)

# Write headers to files
input_file.header_cond('INPUT')
#output_file.header('OUTPUT')

# Write input file with settings
input_file.input_writer_cond(settings, Sources, BCs, domain.T)
input_file.close()
print '################################\n'

print 'Saving data to numpy array files...'
np.save('T_'+'0.000000', domain.T, False)
np.save('eta_'+'0.000000', domain.eta, False)
np.save('X', domain.X, False)
np.save('Y', domain.Y, False)

##########################################################################
# -------------------------------------Solve
##########################################################################
t,nt,tign=0,0,0
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
        print 'Saving data to numpy array files...'
        np.save('T_'+'{:f}'.format(t), domain.T, False)
        np.save('eta_'+'{:f}'.format(t), domain.eta, False)
        break
    
    # Output data to numpy files
    if output_data_nt!=0 and nt%output_data_nt==0:
        print 'Saving data to numpy array files...'
        np.save('T_'+'{:f}'.format(t), domain.T, False)
        np.save('eta_'+'{:f}'.format(t), domain.eta, False)
        
    # Change boundary conditions
    if np.amax(domain.T)>=1200 and not BCs_changed:
#    if np.amax(domain.eta)>=0.50 and not BCs_changed:
        solver.BCs['bc_north']=['C',(30,300),(0,-1)]
        BCs_changed=True
        tign=t
#        break
    
print 'Ignition time: %f ms'%(tign*1000)
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
#fig4=pyplot.figure(figsize=(7, 7))
#pyplot.contourf(domain.X*1000, domain.Y*1000, T, alpha=0.5, cmap=cm.viridis)  
#pyplot.colorbar()
#pyplot.xlabel('$x$ (mm)')
#pyplot.ylabel('$y$ (mm)')
#pyplot.title('Temperature distribution, t=%.7f'%t)
#pyplot.xlim(0,0.4)
#pyplot.ylim(5,6);

#fig4=pyplot.figure(figsize=(7, 7))
#pyplot.contourf(domain.X*1000, domain.Y*1000, eta, alpha=0.5, cmap=cm.viridis)  
#pyplot.colorbar()
#pyplot.xlabel('$x$ (mm)')
#pyplot.ylabel('$y$ (mm)')
#pyplot.title('Reaction progress, t=%.7f'%t)
#pyplot.xlim(0,0.4)
#pyplot.ylim(5,6);
#pyplot.close(fig4)

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