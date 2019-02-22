# -*- coding: utf-8 -*-
"""
######################################################
#             2D Heat Conduction Solver              #
#              Created by J. Mark Epps               #
#          Part of Masters Thesis at UW 2018-2020    #
######################################################

This file contains the main executable script for solving 2D conduction:
    -Uses FileClasses.py to read and write input files to get settings for solver
    and geometry
    -Creates a domain class from GeomClasses.py
    -Creates solver class from SolverClasses.py with reference to domain class
    -Can be called from command line with: 
        python main.py [Input file name+extension] [Output directory relative to current directory]
    -Calculates the time taken to run solver
    -Changes boundary conditions based on ignition criteria
    -Saves temperature data (.npy) at intervals defined in input file
    -Saves x,y meshgrid arrays (.npy) to output directory

Features:
    -Ignition condition met, will change north BC to that of right BC
    -Saves temperature and reaction data (.npy) depending on input file 
    settings

"""

##########################################################################
# ----------------------------------Libraries and classes
##########################################################################
import numpy as np
import string as st
#from datetime import datetime
import os
import sys
import time

from GeomClasses import TwoDimPlanar as TwoDimPlanar
import SolverClasses as Solvers
import FileClasses

print('######################################################')
print('#             2D Heat Conduction Solver              #')
print('#              Created by J. Mark Epps               #')
print('#          Part of Masters Thesis at UW 2018-2020    #')
print('######################################################\n')

# Start timer
time_begin=time.time()

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
domain.E[:,:]=300*5643*599*domain.CV_vol()
print '################################'
##########################################################################
# ------------------------Write Input File settings to output directory
##########################################################################
print 'Saving input file to output directory...'
#datTime=str(datetime.date(datetime.now()))+'_'+'{:%H%M}'.format(datetime.time(datetime.now()))
isBinFile=False

input_file=FileClasses.FileOut('Input_file', isBinFile)

# Write header to file
input_file.header_cond('INPUT')

# Write input file with settings
input_file.input_writer_cond(settings, Sources, BCs)
print '################################\n'

print 'Saving data to numpy array files...'
T=domain.TempFromConserv()
np.save('T_'+'0.000000', T, False)
if st.find(Sources['Source_Kim'],'True')>=0:
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
Sources['Ignition']=st.split(Sources['Ignition'], ',')
Sources['Ignition'][1]=float(Sources['Ignition'][1])
BCs_changed=False

print 'Solving:'
while nt<settings['total_time_steps'] and t<settings['total_time']:
    err,dt=solver.Advance_Soln_Cond(nt, t)
    t+=dt
    nt+=1
    if err==1:
        print '#################### Solver aborted #######################'
        print 'Saving data to numpy array files...'
        T=domain.TempFromConserv()
        np.save('T_'+'{:f}'.format(t), T, False)
        if st.find(Sources['Source_Kim'],'True')>=0:
            np.save('eta_'+'{:f}'.format(t), domain.eta, False)
        break
    
    # Output data to numpy files
    if output_data_nt!=0 and nt%output_data_nt==0:
        print 'Saving data to numpy array files...'
        T=domain.TempFromConserv()
        np.save('T_'+'{:f}'.format(t), T, False)
        if st.find(Sources['Source_Kim'],'True')>=0:
            np.save('eta_'+'{:f}'.format(t), domain.eta, False)
        
    # Change boundary conditions
    T=domain.TempFromConserv()
    if ((Sources['Ignition'][0]=='eta' and np.amax(domain.eta)>=Sources['Ignition'][1])\
        or (Sources['Ignition'][0]=='Temp' and np.amax(T)>=Sources['Ignition'][1]))\
        and not BCs_changed:
        solver.BCs['bc_north']=solver.BCs['bc_right']
        BCs_changed=True
        tign=t
#        break
    
time_end=time.time()
print 'Ignition time: %f ms'%(tign*1000)
print 'Solver time: %f min'%((time_end-time_begin)/60.0)

T, eta=domain.TempFromConserv(), domain.eta

print('Solver has finished its run')