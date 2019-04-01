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
#import string as st
from matplotlib import pyplot, cm
#from mpl_toolkits.mplot3d import Axes3D

pyplot.ioff()

print('######################################################')
print('#            2D Conduction Post-processing           #')
print('#              Created by J. Mark Epps               #')
print('#          Part of Masters Thesis at UW 2018-2020    #')
print('######################################################\n')

inputargs=sys.argv
if len(inputargs)>3:
    dir_files=inputargs[1]
    times=[0]*(len(inputargs)-2)
    for i in range(len(inputargs)-2):
        times[i]=inputargs[i+2]
else:
    print 'Usage is: python Post_zoom.py [Output directory] [Times]\n'
    print 'where\n'
    print '[Output directory] is the directory where the data is located'
    print '[Times] is all the times to be processed (spaced list)'
    print '***********************************'
    sys.exit('Post-processing halted')

try:
    os.chdir(dir_files)
except:
    sys.exit('Directory "'+dir_files+'" not found')
# Get Arrhenius parameters
#A0=-1.0
#Ea=-1.0
#source='False'
#try:
#    input_file=open('Input_file.txt')
#except:
#    try:
#        input_file=open('Input_file_stats.txt')
#    except:
#        sys.exit('Input file missing')

#titles=[]
#while A0<0 or Ea<0 or source=='False':
#    line=input_file.readline()
#    if st.find(line, 'Ea')==0:
#        Ea=float(st.split(line, ':')[1])
#    elif st.find(line, 'A0')==0:
#        A0=float(st.split(line, ':')[1])
#    elif st.find(line, 'Source_Kim')==0:
#        source=st.split(line, ':')[1]
#    elif st.find(line, 'Species')==0:
#        titles=st.split(st.split(st.split(line, ':')[1], '\n')[0], ',')
#input_file.close()

# Generate graphs
X=np.load('X.npy', False)
Y=np.load('Y.npy', False)
for time in times:
    T=np.load('T_'+time+'.npy', False)
    
    # Temperature contour
    fig=pyplot.figure(figsize=(5, 5))
    pyplot.contourf(X*1000, Y*1000, T, alpha=0.5, cmap=cm.viridis)#, vmin=300, vmax=12000)  
    pyplot.colorbar()
    pyplot.xlabel('$x$ (mm)')
    pyplot.ylabel('$y$ (mm)')
    pyplot.xlim([0.8,1.0])
    pyplot.ylim([5.0,6.0])
#    pyplot.clim(300, 1000)
    pyplot.title('Temperature distribution t='+time);
    fig.savefig('T_'+time+'_corner.png',dpi=300)
    pyplot.close(fig)
    print 'Processed '+time
    
print '\nPost-processing complete'