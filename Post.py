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
import os
import sys
import string as st
import matplotlib as mtplt
from matplotlib import pyplot as plt
from FileClasses import FileIn

# Interpolation function for Darcy u calculations
def interpolate(k1, k2, func):
    if func=='Linear':
        return 0.5*k1+0.5*k2
    else:
        return 2*k1*k2/(k1+k2)

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
        elif line[0]=='x_min':
            xmin=float(line[1])
        elif line[0]=='y_min':
            ymin=float(line[1])
        elif line[0]=='x_max':
            try:
                xmax=float(line[1])
            except:
                xmax=line[1]
        elif line[0]=='y_max':
            try:
                ymax=float(line[1])
            except:
                ymax=line[1]
        elif line[0]=='1D_Plots':
            OneD_graphs=line[1]
        elif line[0]=='Darcy_vel':
            darcy=line[1]
        elif line[0]=='Temp_min':
            temp_min=float(line[1])
        elif line[0]=='Temp_max':
            temp_max=float(line[1])
        elif line[0]=='Temp_pts':
            temp_pts=int(line[1])
        elif line[0]=='eta_pts':
            eta_pts=int(line[1])
        elif line[0]=='Phi_Plots':
            Phi_graphs=line[1]

fin.close()

try:
    os.chdir(dir_files)
except:
    sys.exit('Directory "'+dir_files+'" not found')

##############################################################
#               Read Solver file
##############################################################
try:
    input_file=FileIn('Input_file.txt', False)
except:
    sys.exit('Input file missing')

titles=['g','s']
settings={}
sources={}
Species={}
BCs={}

input_file.Read_Input(settings, sources, Species, BCs)
xmax=float(settings['Length'])*1000
ymax=float(settings['Width'])*1000
try:
    settings['rho_IC']=st.split(settings['rho_IC'], ',')
except:
    settings['rho_IC']=float(settings['rho_IC'])
#while A0<0 or Ea<0 or source=='False':
#    line=input_file.readline()
#    if st.find(line, 'Domain')==0:
#        domain=st.split(st.split(line, ':')[1], '\n')[0]
#    elif st.find(line, 'Ea')==0:
#        Ea=float(st.split(line, ':')[1])
#    elif st.find(line, 'A0')==0:
#        A0=float(st.split(line, ':')[1])
#    elif st.find(line, 'Source_Kim')==0:
#        source=st.split(line, ':')[1]
##    elif st.find(line, 'Species')==0:
##        titles=st.split(st.split(st.split(line, ':')[1], '\n')[0], ',')
#    elif st.find(line, 'Length')==0 and type(xmax) is str:
#        xmax=float(st.split(line, ':')[1])*1000
#    elif st.find(line, 'Width')==0 and type(ymax) is str:
#        ymax=float(st.split(line, ':')[1])*1000
#input_file.close()

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
lvl_eta=np.linspace(0, 1, 11)
lvl_temp=np.linspace(temp_min, temp_max, temp_pts)
norm_eta=mtplt.colors.Normalize(vmin=0, vmax=1.0)
fig_size=(6, 6)
cmap_choice=mtplt.cm.viridis
##############################################################
#               Generate graphs
##############################################################
X=np.load('X.npy', False)
Y=np.load('Y.npy', False)
for time in times:
    T=np.load('T_'+time+'.npy', False)
    if st.find(sources['Source_Kim'],'True')>=0:
        eta=np.load('eta_'+time+'.npy', False)
        Y_tot=0
    
    # Temperature contour
    fig=plt.figure(figsize=fig_size)
    plt.contourf(X*1000, Y*1000, T, alpha=0.5, cmap=cmap_choice, extend='both',levels=lvl_temp)#, vmin=270, vmax=2000)  
    cb=plt.colorbar()
    cb.locator=mtplt.ticker.MaxNLocator(nbins=temp_pts)
    cb.update_ticks()
    plt.xlabel('$x$ (mm)')
    plt.ylabel('$y$ (mm)')
#    plt.clim(300, 3000)
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    plt.title('Temperature distribution t='+time+' ms');
    fig.savefig('T_'+time+'.png',dpi=300)
    plt.close(fig)
    
    # 1D temperature profile at centreline
    # if st.find(OneD_graphs,'True')>=0:
        # fig=plt.figure(figsize=fig_size)
        # plt.plot(Y[:,1], T[:,int(len(T[0,:])/2)])
        # plt.xlabel('$y$ (m)')
        # plt.ylabel('T (K)')
#        plt.xlim([xmin,xmax])
#        plt.ylim([temp_min,temp_max])
        # plt.title('Centreline Temperature distribution t='+time)
        # fig.savefig('T_1D_'+time+'.png',dpi=300)
        # plt.close(fig)
    
    if st.find(sources['Source_Kim'],'True')>=0:
        # Progress contour
        fig=plt.figure(figsize=fig_size)
        plt.contourf(X*1000, Y*1000, eta, alpha=0.5, cmap=cmap_choice, levels=lvl_eta)#, vmin=0.0, vmax=1.0)  
        cb=plt.colorbar()
        cb.locator=mtplt.ticker.MaxNLocator(nbins=eta_pts)
        cb.update_ticks()
        plt.xlabel('$x$ (mm)')
        plt.ylabel('$y$ (mm)')
    #    plt.clim(0.0, 1.0)
        plt.xlim([xmin,xmax])
        plt.ylim([ymin,ymax])
        plt.title('Progress distribution t='+time+' ms');
        fig.savefig('eta_'+time+'.png',dpi=300)
        plt.close(fig)
        
        # Reaction rate contour
        if st.find(Phi_graphs,'True')>=0:
            phi=sources['A0']*(1-eta)*np.exp(-sources['Ea']/8.314/T)
            fig=plt.figure(figsize=fig_size)
            plt.contourf(X*1000, Y*1000, phi, alpha=0.5, cmap=cmap_choice)#, vmin=0.0, vmax=1.0)  
            plt.colorbar(format='%.2e')
            plt.xlabel('$x$ (mm)')
            plt.ylabel('$y$ (mm)')
            plt.xlim([xmin,xmax])
            plt.ylim([ymin,ymax])
            plt.title('Reaction rate t='+time+' ms');
            fig.savefig('Phi_'+time+'.png',dpi=300)
            plt.close(fig)
        
        # 1D Reaction rate profile at centreline
        if st.find(OneD_graphs,'True')>=0:
            fig=plt.figure(figsize=fig_size)
            plt.plot(Y[:,1]*1000, phi[:,int(len(T[0,:])/2)])
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            plt.xlabel('$y$ (mm)')
            plt.ylabel('$d\eta/dt$ ($s^{-1}$)')
            plt.xlim([xmin,xmax])
            plt.title('Centreline Reaction rate t='+time+' ms')
            fig.savefig('Phi_1D_'+time+'.png',dpi=300)
            plt.close(fig)
    try:
            # Mass fraction contours
        for i in range(len(titles)):
            Y_0=np.load('rho_'+titles[i]+'_'+time+'.npy', False)
            fig=plt.figure(figsize=fig_size)
            plt.contourf(X*1000, Y*1000, Y_0, alpha=0.5, cmap=cmap_choice)#, vmin=0.0, vmax=1.0)  
            plt.colorbar()
            plt.xlabel('$x$ (mm)')
            plt.ylabel('$y$ (mm)')
        #    plt.clim(0.0, 1.0)
            plt.xlim([xmin,xmax])
            plt.ylim([ymin,ymax])
            plt.title('Density; $'+titles[i]+'$, t='+time+' ms');
            fig.savefig('rho_'+titles[i]+'_'+time+'.png',dpi=300)
            plt.close(fig)
            Y_tot+=np.sum(Y_0)/np.size(Y_0)
    except:
        print 'Processed '+time
        continue
    
    # Darcy velocities and pressure contours
    P=np.load('P_'+time+'.npy', False)
    u=np.zeros_like(P)
    v=np.zeros_like(P)
    por=settings['Porosity']+\
        (1-Y_0/(float(settings['rho_IC'][1])*settings['Porosity']))\
        *(1-settings['Porosity'])
    perm=por**3*settings['Carmen_diam']**2/(settings['Kozeny_const']*(1-por)**2)
    u[:,1:]=-interpolate(perm[:,1:], perm[:,:-1], settings['diff_interpolation'])\
        /settings['Darcy_mu']*(P[:,1:]-P[:,:-1])/(X[:,1:]-X[:,:-1])
    v[1:,:]=-interpolate(perm[1:,:], perm[:-1,:], settings['diff_interpolation'])\
        /settings['Darcy_mu']*(P[1:,:]-P[:-1,:])/(Y[1:,:]-Y[:-1,:])
#    pl=25
    fig=plt.figure(figsize=fig_size)
#    plt.quiver(X[::pl, ::pl]*1000, Y[::pl, ::pl]*1000, \
#                  u[::pl, ::pl], v[::pl, ::pl])
    plt.contourf(X*1000, Y*1000, P, alpha=0.5, cmap=cmap_choice)#, vmin=270, vmax=2000)  
    plt.colorbar()
    plt.xlabel('$x$ (mm)')
    plt.ylabel('$y$ (mm)')
#    plt.clim(300, 10000)
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    plt.title('Pressure t='+time+' ms');
    fig.savefig('P_'+time+'.png',dpi=300)
    plt.close(fig)
    
    # Darcy Velocity contours
    if st.find(darcy, 'True')>=0:
        fig=plt.figure(figsize=fig_size)
        plt.contourf(X*1000, Y*1000, u, alpha=0.5, cmap=cmap_choice)#, vmin=270, vmax=2000)  
        plt.colorbar()
        plt.xlabel('$x$ (mm)')
        plt.ylabel('$y$ (mm)')
    #    plt.clim(300, 10000)
        plt.xlim([xmin,xmax])
        plt.ylim([ymin,ymax])
        plt.title('Darcy Velocity u t='+time+' ms');
        fig.savefig('u_'+time+'.png',dpi=300)
        plt.close(fig)
        
        fig=plt.figure(figsize=fig_size)
        plt.contourf(X*1000, Y*1000, v, alpha=0.5, cmap=cmap_choice)#, vmin=270, vmax=2000)  
        plt.colorbar()
        plt.xlabel('$x$ (mm)')
        plt.ylabel('$y$ (mm)')
    #    plt.clim(300, 10000)
        plt.xlim([xmin,xmax])
        plt.ylim([ymin,ymax])
        plt.title('Darcy Velocity v t='+time+' ms');
        fig.savefig('v_'+time+'.png',dpi=300)
        plt.close(fig)
    
    print 'Processed '+time
    print '     Mass balance residual: %.1f'%(Y_tot)

if st.find(OneD_graphs,'True')>=0:
    print 'Creating 1D plots'
    fig=plt.figure(figsize=fig_size)
    for time in times:
        T=np.load('T_'+time+'.npy', False)
        # 1D temperature profile at centreline
        plt.plot(Y[:,1]*1000, T[:,int(len(T[0,:])/2)], label='t='+time)
    plt.xlabel('$y$ (mm)')
    plt.ylabel('T (K)')
    plt.xlim([xmin,xmax])
    plt.ylim([temp_min,temp_max])
    plt.legend()
    plt.title('Centreline Temperature Evolution')
    fig.savefig('T_1D.png',dpi=300)
    plt.close(fig)
    
    if st.find(sources['Source_Kim'],'True')>=0:
        fig=plt.figure(figsize=fig_size)
        for time in times:
            eta=np.load('eta_'+time+'.npy', False)
            T=np.load('T_'+time+'.npy', False)
            phi=sources['A0']*(1-eta)*np.exp(-sources['Ea']/8.314/T)
            # 1D Reaction rate profile at centreline
            plt.plot(Y[:,1]*1000, phi[:,int(len(T[0,:])/2)], label='t='+time)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.xlabel('$y$ (mm)')
        plt.ylabel('$d\eta/dt$ ($s^{-1}$)')
        plt.xlim([xmin,xmax])
        plt.legend()
        plt.title('Centreline Reaction rate Evolution')
        fig.savefig('Phi_1D.png',dpi=300)
        plt.close(fig)

print '\nPost-processing complete'