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
import matplotlib as mtplt
from matplotlib import pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

plt.ioff()

print('######################################################')
print('#            2D Conduction Post-processing           #')
print('#              Created by J. Mark Epps               #')
print('#          Part of Masters Thesis at UW 2018-2020    #')
print('######################################################\n')

inputargs=sys.argv
if len(inputargs)>1:
    inp_file=inputargs[1]
else:
    print 'Usage is: python Post.py [Input File]\n'
    print 'where\n'
    print '[Input File] is the Post-processing input file'
    print '***********************************'
    sys.exit('Post-processing halted')

##############################################################
#               Read post-processing file
##############################################################
try:
    fin=open(inp_file, 'r')
except:
    sys.exit('Cannot find post-processing input file')
    
for line in fin:
    if st.find(line, ':')>0 and st.find(line, '#')!=0:
        line=st.split(line, ':')
        if line[0]=='Directory':
            dir_files=st.split(line[1], '\n')[0]
        elif line[0]=='Times':
            if st.find(line[1], ',')>0:
                times=st.split(line[1], ',')
                times[-1]=st.split(times[-1], '\n')[0]
            else:
                times=st.split(line[1], '\n')[0]
        elif line[0]=='Temp_min':
            temp_min=float(line[1])
        elif line[0]=='Temp_max':
            temp_max=float(line[1])
        elif line[0]=='Time_Temp_Pos':
            Phi_graphs=st.split(line[1], ',')
            Phi_graphs=[int(Phi_graphs[0]),int(Phi_graphs[1])]
            Phi_graphs=tuple(Phi_graphs)

fin.close()

try:
    os.chdir(dir_files)
except:
    sys.exit('Directory "'+dir_files+'" not found')

##############################################################
#               Times to process (if ALL is selected)
##############################################################
if type(times) is str:
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

##############################################################
#               Figure details
##############################################################
fig_size=(6, 6)
cmap_choice=mtplt.cm.viridis
##############################################################
#               Generate graphs
##############################################################
X=np.load('X.npy', False)
Y=np.load('Y.npy', False)
fig=plt.figure(figsize=fig_size)
for time in times:
    T=np.load('T_'+time+'.npy', False)
    # 1D temperature profile at centreline
    plt.plot(float(time), T[Phi_graphs], marker='o',color='black')
plt.xlabel('$t$ (ms)')
plt.ylabel('T (K)')
#plt.xlim([xmin,xmax])
#plt.ylim([temp_min,temp_max])
plt.legend()
plt.title('Backface Temperature Evolution')
fig.savefig('Temp_time.png',dpi=300)
plt.close(fig)

print '\nPost-processing complete'