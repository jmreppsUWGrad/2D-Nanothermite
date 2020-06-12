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
from myFigs import set_size

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
        elif line[0]=='Variable':
            var=st.split(line[1], '\n')[0]

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
width = 384
fig_size=set_size(width)
#fig_size=(6, 6)
figType='.pdf'
nice_fonts = {
        # Use LaTex to write all text
#        "text.usetex": True,
        "font.family": "serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 11,
        "font.size": 11,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
}
mtplt.rcParams.update(nice_fonts)
cmap_choice=mtplt.cm.viridis
y_label={'Temperature': 'T [$K$]', 'Pressure': 'P [$Pa$]',\
         'eta': '$\eta$ [-]', 'rho_g': 'Density [$kg/m^3$]',\
         'rho_s': 'Density [$kg/m^3$]'}
var_name={'Temperature': 'T', 'Pressure': 'P',\
         'eta': 'eta', 'rho_g': 'rho_g', 'rho_s': 'rho_s'}
##############################################################
#               Generate graphs
##############################################################
X=np.load('X.npy', False)
Y=np.load('Y.npy', False)
#fig=plt.figure(figsize=fig_size)
fig,ax1=plt.subplots()
ax1.set_xlabel('Time [ms]')
ax1.set_ylabel(y_label[var], color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set_ylabel('T [$K$]', color='red')  # we already handled the x-label with ax1
#ax2.ylim([-100,3000])
ax2.tick_params(axis='y', labelcolor='red')

for time in times:
    T=np.load(var_name[var]+'_'+time+'.npy', False)
    T2=np.load('T_'+time+'.npy', False)
    # 1D temperature profile at centreline
    ax1.plot(float(time), T[Phi_graphs], marker='o',color='black')
    ax2.plot(float(time), T2[Phi_graphs], marker='^',color='red')
#plt.xlabel('Time [ms]')
#plt.ylabel(y_label[var])
#plt.xlim([xmin,xmax])
#plt.ylim([temp_min,temp_max])
#plt.legend()
#plt.title(var+' evolution with time at position '+str(Phi_graphs))
fig.tight_layout()  # otherwise the right y-label is slightly clipped
fig.savefig(var+'_time_'+str(Phi_graphs)+figType,dpi=300)
plt.close(fig)

print '\nPost-processing complete'