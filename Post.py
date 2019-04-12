# -*- coding: utf-8 -*-
"""
######################################################
#             2D Heat Conduction Solver              #
#              Created by J. Mark Epps               #
#          Part of Masters Thesis at UW 2018-2020    #
######################################################

This file contains the post-processing script for 2D conduction:
    -Called from command line by:
        python Post-processing.py [Data directory relative to current directory]
    -Reads input file to get necessary parameters
    -Reads x,y meshgrid arrays (.npy) for graph output
    -Reads variable arrays (.npy files) and outputs graphs (.png) for 
    each time step in directory

Features:
    -Graphs of Temperature, reaction progress, reaction rate

Desired:
    -display variables available OR pass in desired graphs as argument for
    executing script
    
"""

import numpy as np
#import CoolProp.CoolProp as CP
import os
import sys
import string as st
from matplotlib import pyplot, cm
#from mpl_toolkits.mplot3d import Axes3D

pyplot.ioff()

print('######################################################')
print('#            2D Conduction Post-processing           #')
print('#              Created by J. Mark Epps               #')
print('#          Part of Masters Thesis at UW 2018-2020    #')
print('######################################################\n')

inputargs=sys.argv
if len(inputargs)>2:
    dir_files=inputargs[1]
    OneD_graphs=int(inputargs[2])
else:
    print 'Usage is: python Post.py [Output directory] [1D graphs]\n'
    print 'where\n'
    print '[Output directory] is the directory where the data is located'
    print '[1D graphs] indicates whether 1D graphs should be output (1 or 0); default is 0'
    print '***********************************'
    sys.exit('Post-processing halted')

try:
    os.chdir(dir_files)
except:
    sys.exit('Directory "'+dir_files+'" not found')
# Get Arrhenius parameters
A0=-1.0
Ea=-1.0
source='False'
try:
    input_file=open('Input_file.txt')
except:
    try:
        input_file=open('Input_file_stats.txt')
    except:
        sys.exit('Input file missing')

titles=[]
while A0<0 or Ea<0 or source=='False':
    line=input_file.readline()
    if st.find(line, 'Ea')==0:
        Ea=float(st.split(line, ':')[1])
    elif st.find(line, 'A0')==0:
        A0=float(st.split(line, ':')[1])
    elif st.find(line, 'Source_Kim')==0:
        source=st.split(line, ':')[1]
    elif st.find(line, 'Species')==0:
        titles=st.split(st.split(st.split(line, ':')[1], '\n')[0], ',')
input_file.close()

# Get times to process
times=os.listdir('.')
i=len(times)
j=0
while i>j:
    if st.find(times[j],'T')==0 and st.find(times[j],'.npy')>0:
        times[j]=st.split(st.split(times[j],'_')[1],'.npy')[0]
        j+=1
    else:
        del times[j]
        i-=1

# Graph spacial limits
xmin,xmax=0,0.001
ymin,ymax=0,0.006
#xmin,xmax=0.0004,0.0006
#ymin,ymax=0.005,0.006

# Generate graphs
X=np.load('X.npy', False)
Y=np.load('Y.npy', False)
for time in times:
    T=np.load('T_'+time+'.npy', False)
    if st.find(source,'True')>=0:
        eta=np.load('eta_'+time+'.npy', False)
        Y_tot=np.zeros_like(Y)
    
    # Temperature contour
    fig=pyplot.figure(figsize=(6, 6))
    pyplot.contourf(X, Y, T, alpha=0.5, cmap=cm.viridis)#, vmin=270, vmax=2000)  
    pyplot.colorbar()
    pyplot.xlabel('$x$ (m)')
    pyplot.ylabel('$y$ (m)')
#    pyplot.clim(300, 10000)
    pyplot.xlim([xmin,xmax])
    pyplot.ylim([ymin,ymax])
    pyplot.title('Temperature distribution t='+time);
    fig.savefig('T_'+time+'.png',dpi=300)
    pyplot.close(fig)
    
    # 1D temperature profile at centreline
    # if OneD_graphs==1:
        # fig=pyplot.figure(figsize=(6, 6))
        # pyplot.plot(Y[:,1], T[:,int(len(T[0,:])/2)])
        # pyplot.xlabel('$y$ (m)')
        # pyplot.ylabel('T (K)')
        # pyplot.title('Centreline Temperature distribution t='+time)
        # fig.savefig('T_1D_'+time+'.png',dpi=300)
        # pyplot.close(fig)
    
    if st.find(source,'True')>=0:
        # Progress contour
        fig=pyplot.figure(figsize=(6, 6))
        pyplot.contourf(X, Y, eta, alpha=0.5, cmap=cm.viridis)#, vmin=0.0, vmax=1.0)  
        pyplot.colorbar()
        pyplot.xlabel('$x$ (m)')
        pyplot.ylabel('$y$ (m)')
    #    pyplot.clim(0.0, 1.0)
        pyplot.xlim([xmin,xmax])
        pyplot.ylim([ymin,ymax])
        pyplot.title('Progress distribution t='+time);
        fig.savefig('eta_'+time+'.png',dpi=300)
        pyplot.close(fig)
        
        # Reaction rate contour
        phi=A0*(1-eta)*np.exp(-Ea/8.314/T)
        fig=pyplot.figure(figsize=(6, 6))
        pyplot.contourf(X, Y, phi, alpha=0.5, cmap=cm.viridis)#, vmin=0.0, vmax=1.0)  
        pyplot.colorbar(format='%.2e')
        pyplot.xlabel('$x$ (m)')
        pyplot.ylabel('$y$ (m)')
    #    pyplot.clim(0.0, 1.0)
        pyplot.xlim([xmin,xmax])
        pyplot.ylim([ymin,ymax])
        pyplot.title('Reaction rate t='+time);
        fig.savefig('Phi_'+time+'.png',dpi=300)
        pyplot.close(fig)
        
        # 1D Reaction rate profile at centreline
        if OneD_graphs==1:
            fig=pyplot.figure(figsize=(6, 6))
            pyplot.plot(Y[:,1], phi[:,int(len(T[0,:])/2)])
            pyplot.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            pyplot.xlabel('$y$ (m)')
            pyplot.ylabel('$d\eta/dt$ ($s^{-1}$)')
            pyplot.title('Centreline Reaction rate t='+time)
            fig.savefig('Phi_1D_'+time+'.png',dpi=300)
            pyplot.close(fig)
    try:
        P=np.load('P_'+time+'.npy', False)
        fig=pyplot.figure(figsize=(6, 6))
        pyplot.contourf(X, Y, P, alpha=0.5, cmap=cm.viridis)#, vmin=270, vmax=2000)  
        pyplot.colorbar()
        pyplot.xlabel('$x$ (m)')
        pyplot.ylabel('$y$ (m)')
    #    pyplot.clim(300, 10000)
        pyplot.xlim([xmin,xmax])
        pyplot.ylim([ymin,ymax])
        pyplot.title('Pressure t='+time);
        fig.savefig('P_'+time+'.png',dpi=300)
        pyplot.close(fig)
    except:
        print 'Processed '+time
        continue
    
        # Mass fraction contours
    for i in range(len(titles)):
        Y_0=np.load('m_'+titles[i]+'_'+time+'.npy', False)
        fig=pyplot.figure(figsize=(6, 6))
        pyplot.contourf(X, Y, Y_0, alpha=0.5, cmap=cm.viridis)#, vmin=0.0, vmax=1.0)  
        pyplot.colorbar()
        pyplot.xlabel('$x$ (m)')
        pyplot.ylabel('$y$ (m)')
    #    pyplot.clim(0.0, 1.0)
        pyplot.xlim([xmin,xmax])
        pyplot.ylim([ymin,ymax])
        pyplot.title('Mass; $'+titles[i]+'$, t='+time);
        fig.savefig('m_'+titles[i]+'_'+time+'.png',dpi=300)
        pyplot.close(fig)
        Y_tot+=Y_0
            
        
    print 'Processed '+time
    print '     Mass balance residual: %8f'%(np.amin(Y_tot)*10**6)

if OneD_graphs==1:
    print 'Creating 1D plots'
    fig=pyplot.figure(figsize=(6, 6))
    for time in times:
        T=np.load('T_'+time+'.npy', False)
        # 1D temperature profile at centreline
        pyplot.plot(Y[:,1], T[:,int(len(T[0,:])/2)], label='t='+time)
    pyplot.xlabel('$y$ (m)')
    pyplot.ylabel('T (K)')
    pyplot.legend()
    pyplot.title('Centreline Temperature Evolution')
    fig.savefig('T_1D.png',dpi=300)
    pyplot.close(fig)
    
    if st.find(source,'True')>=0:
        fig=pyplot.figure(figsize=(6, 6))
        for time in times:
            eta=np.load('eta_'+time+'.npy', False)
            T=np.load('T_'+time+'.npy', False)
            phi=A0*(1-eta)*np.exp(-Ea/8.314/T)
            # 1D Reaction rate profile at centreline
            pyplot.plot(Y[:,1], phi[:,int(len(T[0,:])/2)], label='t='+time)
        pyplot.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        pyplot.xlabel('$y$ (m)')
        pyplot.ylabel('$d\eta/dt$ ($s^{-1}$)')
        pyplot.legend()
        pyplot.title('Centreline Reaction rate Evolution')
        fig.savefig('Phi_1D.png',dpi=300)
        pyplot.close(fig)

print '\nPost-processing complete'