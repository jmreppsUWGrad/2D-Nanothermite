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

import GeomClasses as Geom
import SolverClasses as Solvers
import FileClasses

def save_data(Domain, Sources, Species, time):
    T=Domain.TempFromConserv()
    P=Domain.P
    np.save('T_'+time, T, False)
    np.save('P_'+time, P, False)
    if st.find(Sources['Source_Kim'],'True')>=0:
        np.save('eta_'+time, Domain.eta, False)
    if bool(Species):
        for i in Species['keys']:
            np.save('m_'+i+'_'+time, Domain.m_species[i], False)

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
Species={}
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
fin.Read_Input(settings, Sources, Species, BCs)
try:
    os.chdir(settings['Output_directory'])
except:
    os.makedirs(settings['Output_directory'])
    os.chdir(settings['Output_directory'])

print '################################'

# Initial conditions from previous run/already in memory
#Use_inital_values                   = False


print 'Initializing geometry package...'
if settings['Domain']=='Planar':
    domain=Geom.TwoDimPlanar(settings, Species, 'Solid')
elif settings['Domain']=='Axisymmetric':
    domain=Geom.AxisymDomain(settings, Species, 'Solid')

domain.mesh()
print '################################'

##########################################################################
# -------------------------------------Initialize solver and domain
##########################################################################

print 'Initializing solver package...'
if settings['Domain']=='Planar':
    solver=Solvers.TwoDimPlanarSolve(domain, settings, Sources, BCs, 'Solid')
elif settings['Domain']=='Axisymmetric':
    solver=Solvers.AxisymmetricSolve(domain, settings, Sources, BCs, 'Solid')
print '################################'

print 'Initializing domain...'
time_max='0.000000'
T=300*np.ones_like(domain.E)
# Restart from previous data
if type(settings['Restart']) is int:
    times=os.listdir('.')
    i=len(times)
    if i<2:
        sys.exit('Cannot find a file to restart a simulation with')
    j=0
    while i>j:
        if st.find(times[j],'T')==0 and st.find(times[j],'.npy')>0 \
            and st.find(times[j],str(settings['Restart']))>=0:
            times[j]=st.split(st.split(times[j],'_')[1],'.npy')[0]
#            if st.find(times[j],str(settings['Restart']))>=0:
            time_max=times[j]
            j+=1
            break
        else:
            del times[j]
            i-=1
    T=np.load('T_'+time_max+'.npy')
    if st.find(Sources['Source_Kim'],'True')>=0:
        domain.eta=np.load('eta_'+time_max+'.npy')

k,rho,Cv,D=domain.calcProp()
vol=domain.CV_vol()
domain.E[:,:]=rho*Cv*vol*T
del k,rho,Cv,D,T
if bool(domain.m_species):
    for i in range(len(Species['Species'])):
        domain.m_species[Species['Species'][i]][:,:]=Species['Species_IC'][i]
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
input_file.input_writer_cond(settings, Sources, Species, BCs)
print '################################\n'

print 'Saving data to numpy array files...'
save_data(domain, Sources, Species, time_max)
np.save('X', domain.X, False)
np.save('Y', domain.Y, False)

##########################################################################
# -------------------------------------Solve
##########################################################################
t,nt,tign=float(time_max),0,0 # time, number steps and ignition time initializations
v_0,v_1,v,N=0,0,0,0 # combustion wave speed variables initialization

# Setup intervals to save data
output_data_t,output_data_nt=0,0
if settings['total_time_steps']=='None':
    output_data_t=settings['total_time']/settings['Number_Data_Output']
    settings['total_time_steps']=settings['total_time']*10**9
    t_inc=int(t/output_data_t)+1
elif settings['total_time']=='None':
    output_data_nt=int(settings['total_time_steps']/settings['Number_Data_Output'])
    settings['total_time']=settings['total_time_steps']*10**9
    t_inc=0

# Ignition conditions
Sources['Ignition']=st.split(Sources['Ignition'], ',')
Sources['Ignition'][1]=float(Sources['Ignition'][1])
BCs_changed=False

print 'Solving:'
while nt<settings['total_time_steps'] and t<settings['total_time']:
    # First point in calculating combustion propagation speed
    T_0=domain.TempFromConserv()
    if st.find(Sources['Source_Kim'],'True')>=0 and BCs_changed:
#        v_0=np.sum(domain.eta[:,int(len(domain.eta[0,:])/2)]*domain.dy)
        v_0=np.sum(domain.eta*solver.dy)/len(domain.eta[0,:])
    err,dt=solver.Advance_Soln_Cond(nt, t, vol)
    t+=dt
    nt+=1
    if err>0:
        print '#################### Solver aborted #######################'
        print 'Saving data to numpy array files...'
        save_data(domain, Sources, Species, '{:f}'.format(t))
        input_file.Write_single_line('#################### Solver aborted #######################')
        input_file.Write_single_line('Time step %i, Time elapsed=%f, error code=%i;'%(nt,t,err))
        input_file.Write_single_line('Error codes: 1-time step, 2-Energy, 3-reaction progress')
        input_file.Write_single_line('4-Species balance\n')
        break
    
    # Output data to numpy files
    if (output_data_nt!=0 and nt%output_data_nt==0) or \
        (output_data_t!=0 and (t>=output_data_t*t_inc and t-dt<output_data_t*t_inc)):
        print 'Saving data to numpy array files...'
        save_data(domain, Sources, Species, '{:f}'.format(t))
        t_inc+=1
        
    # Change boundary conditions
    T=domain.TempFromConserv()
    if ((Sources['Ignition'][0]=='eta' and np.amax(domain.eta)>=Sources['Ignition'][1])\
        or (Sources['Ignition'][0]=='Temp' and np.amax(T)>=Sources['Ignition'][1]))\
        and not BCs_changed:
        solver.BCs['bc_north']=solver.BCs['bc_right']
        input_file.fout.write('##bc_north_new:')
        input_file.Write_single_line(str(solver.BCs['bc_north']))
        input_file.fout.write('\n')
        BCs_changed=True
        tign=t
        save_data(domain, Sources, Species, '{:f}'.format(t))
#    if not BCs_changed:
#        k,rho,Cv=domain.calcProp()
#        T_theo=300+2*solver.BCs['bc_north'][1]/k[-1,0]\
#            *np.sqrt(k[-1,0]/rho[-1,0]/Cv[-1,0]*t/np.pi)
#        dT_theo=solver.BCs['bc_north'][1]/k[-1,0]\
#            *np.sqrt(k[-1,0]/rho[-1,0]/Cv[-1,0]/np.pi/t)
##        if (T_theo-T[-1,0])/T_theo<-0.3:
#        if (dT_theo-((T[-1,0]-T_0[-1,0])/dt))/dT_theo<-1.0:
#            solver.BCs['bc_north']=solver.BCs['bc_right']
#            BCs_changed=True
#            tign=t
#            save_data(domain, Sources, '{:f}'.format(t))
    
    # Second point in calculating combustion propagation speed
    if st.find(Sources['Source_Kim'],'True')>=0 and BCs_changed:
#        v_1=np.sum(domain.eta[:,int(len(domain.eta[0,:])/2)]*domain.dy)
        v_1=np.sum(domain.eta*solver.dy)/len(domain.eta[0,:])
        if (v_1-v_0)/dt>0.001:
            v+=(v_1-v_0)/dt
            N+=1
        
time_end=time.time()
input_file.Write_single_line('Final time step size: %f ms'%(dt*1000))
print 'Ignition time: %f ms'%(tign*1000)
input_file.Write_single_line('Ignition time: %f ms'%(tign*1000))
print 'Solver time per 1000 time steps: %f min'%((time_end-time_begin)/60.0*1000/nt)
input_file.Write_single_line('Solver time per 1000 time steps: %f min'%((time_end-time_begin)/60.0*1000/nt))
try:
    print 'Average wave speed: %f m/s'%(v/N)
    input_file.Write_single_line('Average wave speed: %f m/s'%(v/N))
    input_file.close()
except:
    print 'Average wave speed: 0 m/s'
    input_file.Write_single_line('Average wave speed: 0 m/s')
    input_file.close()

print('Solver has finished its run')