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
from GeomClasses import TwoDimDomain
from Source_Comb import Source_terms
from myFigs import set_size

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
        elif line[0]=='Contour_Plots':
            contours=line[1]

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

# Get ignition delay and burn rate
t_ign,v_BR='0','0'
while (type(t_ign) is str) or (type(v_BR) is str):
    dat=input_file.fin.readline()
    if st.find(dat, 'Average wave')>=0:
        v_BR=float(st.split(st.split(dat, ':')[1], 'm')[0])
    if st.find(dat, 'Ignition time')>=0:
        t_ign=float(st.split(st.split(dat, ':')[1], 'm')[0])
input_file.fin.seek(0)

# Get settings
input_file.Read_Input(settings, sources, Species, BCs)
xmax=float(settings['Length'])*1000
ymax=float(settings['Width'])*1000
#try:
#    settings['rho_IC']=st.split(settings['rho_IC'], ',')
#except:
#    settings['rho_IC']=float(settings['rho_IC'])
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
#width = 384
#fig_size=set_size(width)
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
x_axis_labels={
        'Planar': '$x$ [mm]',
        'Axisymmetric': '$r$ [mm]'
        }
y_axis_labels={
        'Planar': '$y$ [mm]',
        'Axisymmetric': '$z$ [mm]'
        }
mtplt.rcParams.update(nice_fonts)
cmap_choice=mtplt.cm.viridis
##############################################################
#               Post-processing
##############################################################
X=np.load('X.npy', False)
Y=np.load('Y.npy', False)
# Open post-processing file
fout=open('Post_processing.txt', 'a')
fout.write('Post processing results:\n')
fout.write(dir_files+'\n\n')
# Initialize geometry and source objects
geom=TwoDimDomain(settings, Species, settings['Domain'], 0)
geom.mesh()
geom.create_var(Species)
dx,dy=np.meshgrid(geom.dx, geom.dy)
hx,hy=geom.CV_dim()
source=Source_terms(sources['Ea'],sources['A0'],sources['dH'],sources['gas_gen'])
# Initialize variables
rho_avg=0
u_avg=0
p_max=0
for time in times:
    fout.write('Time = '+str(time)+'\n')
    ##############################################################
    #               Generate graphs
    ##############################################################
    T=np.load('T_'+time+'.npy', False)
    if st.find(sources['Source_Kim'],'True')>=0:
        eta=np.load('eta_'+time+'.npy', False)
        Y_tot=np.zeros_like(eta)
    
    # Temperature contour
    if st.find(contours,'True')>=0:
        fig=plt.figure(figsize=fig_size)
        plt.contourf(X*1000, Y*1000, T, alpha=0.5, cmap=cmap_choice, extend='both',levels=lvl_temp)#, vmin=270, vmax=2000)  
        cb=plt.colorbar()
        cb.locator=mtplt.ticker.MaxNLocator(nbins=temp_pts)
        cb.update_ticks()
        plt.xlabel(x_axis_labels[settings['Domain']])
        plt.ylabel(y_axis_labels[settings['Domain']])
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
    
    if st.find(sources['Source_Kim'],'True')>=0 and st.find(contours,'True')>=0:
        # Progress contour
        fig=plt.figure(figsize=fig_size)
        plt.contourf(X*1000, Y*1000, eta, alpha=0.5, cmap=cmap_choice, levels=lvl_eta)#, vmin=0.0, vmax=1.0)  
        cb=plt.colorbar()
        cb.locator=mtplt.ticker.MaxNLocator(nbins=eta_pts)
        cb.update_ticks()
        plt.xlabel(x_axis_labels[settings['Domain']])
        plt.ylabel(y_axis_labels[settings['Domain']])
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
            plt.xlabel(x_axis_labels[settings['Domain']])
            plt.ylabel(y_axis_labels[settings['Domain']])
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
            plt.xlabel(y_axis_labels[settings['Domain']])
            plt.ylabel('$d\eta/dt$ ($s^{-1}$)')
            plt.xlim([xmin,xmax])
            plt.title('Centreline Reaction rate t='+time+' ms')
            fig.savefig('Phi_1D_'+time+'.png',dpi=300)
            plt.close(fig)
    try:
        Y_0=[]
            # Mass fraction contours
        for i in range(len(titles)):
            Y_0.append(np.load('rho_'+titles[i]+'_'+time+'.npy', False))
            if st.find(contours,'True')>=0:
                fig=plt.figure(figsize=fig_size)
                plt.contourf(X*1000, Y*1000, Y_0[i], alpha=0.5, cmap=cmap_choice)#, vmin=0.0, vmax=1.0)  
                plt.colorbar()
                plt.xlabel(x_axis_labels[settings['Domain']])
                plt.ylabel(y_axis_labels[settings['Domain']])
            #    plt.clim(0.0, 1.0)
                plt.xlim([xmin,xmax])
                plt.ylim([ymin,ymax])
                plt.title('Density; $'+titles[i]+'$, t='+time+' ms');
                fig.savefig('rho_'+titles[i]+'_'+time+'.png',dpi=300)
                plt.close(fig)
            Y_tot+=Y_0[i]
    except:
        print 'Processed '+time
        continue
    
    # Darcy velocities and pressure contours
    P=np.load('P_'+time+'.npy', False)
    p_max=max(np.amax(P),p_max)
    u=np.zeros_like(P)
    v=np.zeros_like(P)
    por=np.ones_like(P)*settings['Porosity']
#    por=settings['Porosity']+\
#        (1-Y_0/(float(settings['rho_IC'][1])*settings['Porosity']))\
#        *(1-settings['Porosity'])
    perm=por**3*settings['Carmen_diam']**2/(settings['Kozeny_const']*(1-por)**2)
    u[:,1:]=-interpolate(perm[:,1:], perm[:,:-1], settings['diff_interpolation'])\
        /settings['Darcy_mu']*(P[:,1:]-P[:,:-1])/(X[:,1:]-X[:,:-1])
    v[1:,:]=-interpolate(perm[1:,:], perm[:-1,:], settings['diff_interpolation'])\
        /settings['Darcy_mu']*(P[1:,:]-P[:-1,:])/(Y[1:,:]-Y[:-1,:])
#    pl=25
#    P_plot=np.zeros_like(P)
#    P_plot[:,1:]=P[:,1:]-P[:,:-1]# Pressure difference
#    P_plot[1:,:]=P[1:,:]-P[:-1,:]
    if st.find(contours,'True')>=0:
        fig=plt.figure(figsize=fig_size)
    #    plt.quiver(X[::pl, ::pl]*1000, Y[::pl, ::pl]*1000, \
    #                  u[::pl, ::pl], v[::pl, ::pl])
        plt.contourf(X*1000, Y*1000, P, alpha=0.5, cmap=cmap_choice)#, vmin=270, vmax=2000)  
        plt.colorbar()
        plt.xlabel(x_axis_labels[settings['Domain']])
        plt.ylabel(y_axis_labels[settings['Domain']])
    #    plt.clim(300, 10000)
        plt.xlim([xmin,xmax])
        plt.ylim([ymin,ymax])
#        plt.title('Pressure t='+time+' ms');
        fig.savefig('P_'+time+'.png',dpi=300)
        plt.close(fig)
    
    # Darcy Velocity contours
    if st.find(darcy, 'True')>=0 and st.find(contours,'True')>=0:
        fig=plt.figure(figsize=fig_size)
        plt.contourf(X*1000, Y*1000, u, alpha=0.5, cmap=cmap_choice)#, vmin=270, vmax=2000)  
        plt.colorbar()
        plt.xlabel(x_axis_labels[settings['Domain']])
        plt.ylabel(y_axis_labels[settings['Domain']])
    #    plt.clim(300, 10000)
        plt.xlim([xmin,xmax])
        plt.ylim([ymin,ymax])
#        plt.title('Darcy Velocity u t='+time+' ms');
        fig.savefig('u_'+time+'.png',dpi=300)
        plt.close(fig)
        
        lvl_v=np.linspace(-25, 25, temp_pts+1)
        fig=plt.figure(figsize=fig_size)
        plt.contourf(X*1000, Y*1000, v, alpha=0.5, cmap=cmap_choice, extend='both',levels=lvl_v)#, vmin=270, vmax=2000)  
        cb=plt.colorbar()
        cb.locator=mtplt.ticker.MaxNLocator(nbins=temp_pts)
        cb.update_ticks()
        plt.xlabel(x_axis_labels[settings['Domain']])
        plt.ylabel(y_axis_labels[settings['Domain']])
    #    plt.clim(300, 10000)
        plt.xlim([xmin,xmax])
        plt.ylim([ymin,ymax])
#        plt.title('Darcy Velocity v t='+time+' ms');
        fig.savefig('v_'+time+'.pdf')
        plt.close(fig)
    
    print 'Processed '+time
    fout.write('     Mass balance residual: %.3f'%(\
                          np.sum(Y_tot-(float(st.split(settings['rho_IC'], ',')[0])*settings['Porosity']+\
                            float(st.split(settings['rho_IC'], ',')[1])*(1-settings['Porosity'])))/\
                            (np.size(Y_tot)*(float(st.split(settings['rho_IC'], ',')[0])*settings['Porosity']+\
                            float(st.split(settings['rho_IC'], ',')[1])*(1-settings['Porosity']))))+'\n')
    if st.find(darcy, 'True')>=0:
        fout.write('     Max velocity u: %.1f'%(np.amax(abs(u)))+'\n')
        fout.write('     Max velocity v: %.1f'%(np.amax(abs(v)))+'\n')
    
    ##############################################################
    #               Non-dimensional numbers
    ##############################################################
    # T is temp, P is press, Y_0 is rho_spec, eta is eta, settings dicts
    # Update geometry object
    geom.eta=eta
    geom.T_guess=T
    rhoC=geom.calcProp(T, True)
    geom.E=rhoC*T
    geom.rho_species['s']=Y_0[1]
    geom.rho_species['g']=Y_0[0]
    conv=np.zeros_like(X)
    cond=np.zeros_like(X)
    T2, k, rhoC, Cp=geom.calcProp(T)
    
    # Axisymmetric domain flux in r
    if geom.type=='Axisymmetric':
        
        # Left face
        cond[:,1:-1]   -= 1.0/hx[:,1:-1]/(X[:,1:-1])\
                    *(X[:,1:-1]-dx[:,:-2]/2)\
                    *interpolate(k[:,:-2],k[:,1:-1], settings['diff_interpolation'])\
                    *(T[:,1:-1]-T[:,:-2])/dx[:,:-2]
        cond[:,-1]   -= 1.0/hx[:,-1]/(X[:,-1]-dx[:,-1]/2)\
                    *(X[:,-1]-dx[:,-1]/2)\
                    *interpolate(k[:,-2],k[:,-1], settings['diff_interpolation'])\
                    *(T[:,-1]-T[:,-2])/dx[:,-1]
        conv[:,1:-1]+=1.0/hx[:,1:-1]/(X[:,1:-1])\
            *(X[:,1:-1]-dx[:,:-2]/2)\
            *interpolate(Y_0[0][:,1:-1],Y_0[0][:,:-2],settings['conv_interpolation'])\
            *(-interpolate(perm[:,1:-1],perm[:,:-2], settings['diff_interpolation'])/geom.mu\
            *(P[:,1:-1]-P[:,:-2])/dx[:,:-2])\
            *interpolate(Cp[:,1:-1],Cp[:,:-2],settings['conv_interpolation'])\
            *interpolate(T[:,1:-1],T[:,:-2],settings['conv_interpolation'])
        conv[:,-1]+=1.0/hx[:,-1]/(X[:,-1]-dx[:,-1]/2)\
            *(X[:,-1]-dx[:,-2]/2)\
            *interpolate(Y_0[0][:,-1],Y_0[0][:,-2],settings['conv_interpolation'])\
            *(-interpolate(perm[:,-1],perm[:,-2],settings['diff_interpolation'])/geom.mu\
            *(P[:,-1]-P[:,-2])/dx[:,-2])\
            *interpolate(Cp[:,-1],Cp[:,-2],settings['conv_interpolation'])\
            *interpolate(T[:,-1],T[:,-2],settings['conv_interpolation'])
        # Right face
        cond[:,1:-1] += 1.0/hx[:,1:-1]/(X[:,1:-1])\
                    *(X[:,1:-1]+dx[:,1:-1]/2)\
                    *interpolate(k[:,1:-1],k[:,2:], settings['diff_interpolation'])\
                    *(T[:,2:]-T[:,1:-1])/dx[:,1:-1]
        cond[:,0] += 1.0/hx[:,0]/(dx[:,0]/2)\
                    *(X[:,0]+dx[:,0]/2)\
                    *interpolate(k[:,0],k[:,1], settings['diff_interpolation'])\
                    *(T[:,1]-T[:,0])/dx[:,0]
        conv[:,1:-1]-=1.0/hx[:,1:-1]/(X[:,1:-1])\
            *(X[:,1:-1]+dx[:,1:-1]/2)\
            *interpolate(Y_0[0][:,2:],Y_0[0][:,1:-1],settings['conv_interpolation'])\
            *(-interpolate(perm[:,2:],perm[:,1:-1],settings['diff_interpolation'])/geom.mu\
            *(P[:,2:]-P[:,1:-1])/dx[:,1:-1])\
            *interpolate(Cp[:,2:],Cp[:,1:-1],settings['conv_interpolation'])\
            *interpolate(T[:,2:],T[:,1:-1],settings['conv_interpolation'])
        conv[:,0]-=1.0/hx[:,0]/(dx[:,0]/2)\
            *(X[:,0]+dx[:,0]/2)\
            *interpolate(Y_0[0][:,1],Y_0[0][:,0],settings['conv_interpolation'])\
            *(-interpolate(perm[:,1],perm[:,0],settings['diff_interpolation'])/geom.mu\
            *(P[:,1]-P[:,0])/dx[:,0])\
            *interpolate(Cp[:,1],Cp[:,0],settings['conv_interpolation'])\
            *interpolate(T[:,1],T[:,0],settings['conv_interpolation'])
    # Planar domain flux in r
    else:
        # Left face
        cond[:,1:]   -= 1.0/hx[:,1:]\
                    *interpolate(k[:,:-1],k[:,1:], settings['diff_interpolation'])\
                    *(T[:,1:]-T[:,:-1])/dx[:,:-1]
        conv[:,1:]+=1.0/hx[:,1:]\
            *interpolate(Y_0[0][:,1:],Y_0[0][:,:-1],settings['conv_interpolation'])\
            *(-interpolate(perm[:,1:],perm[:,:-1],settings['diff_interpolation'])/geom.mu\
            *(P[:,1:]-P[:,:-1])/dx[:,:-1])\
            *interpolate(Cp[:,1:],Cp[:,:-1],settings['conv_interpolation'])\
            *interpolate(T[:,1:],T[:,:-1],settings['conv_interpolation'])
        # Right face
        cond[:,:-1] += 1.0/hx[:,:-1]\
                    *interpolate(k[:,:-1],k[:,1:], settings['diff_interpolation'])\
                    *(T[:,1:]-T[:,:-1])/dx[:,:-1]
        conv[:,:-1]-=1.0/hx[:,:-1]\
            *interpolate(Y_0[0][:,1:],Y_0[0][:,:-1],settings['conv_interpolation'])\
            *(-interpolate(perm[:,1:],perm[:,:-1],settings['diff_interpolation'])/geom.mu\
            *(P[:,1:]-P[:,:-1])/dx[:,:-1])\
            *interpolate(Cp[:,1:],Cp[:,:-1],settings['conv_interpolation'])\
            *interpolate(T[:,1:],T[:,:-1],settings['conv_interpolation'])
    
    # South face
    cond[1:,:]   -= 1.0/hy[1:,:]\
                *interpolate(k[1:,:],k[:-1,:], settings['diff_interpolation'])\
                *(T[1:,:]-T[:-1,:])/dy[:-1,:]
    conv[1:,:]+=1.0/hy[1:,:]\
        *interpolate(Y_0[0][1:,:],Y_0[0][:-1,:],settings['conv_interpolation'])\
        *(-interpolate(perm[1:,:],perm[:-1,:],settings['diff_interpolation'])/geom.mu\
        *(P[1:,:]-P[:-1,:])/dy[:-1,:])\
        *interpolate(Cp[1:,:],Cp[:-1,:],settings['conv_interpolation'])\
        *interpolate(T[1:,:],T[:-1,:],settings['conv_interpolation'])
    # North face
    cond[:-1,:]  += 1.0/hy[:-1,:]\
                *interpolate(k[:-1,:],k[1:,:], settings['diff_interpolation'])\
                *(T[1:,:]-T[:-1,:])/dy[:-1,:]
    conv[:-1,:]-=1.0/hy[:-1,:]\
        *interpolate(Y_0[0][1:,:],Y_0[0][:-1,:],settings['conv_interpolation'])\
        *(-interpolate(perm[1:,:],perm[:-1,:],settings['diff_interpolation'])/geom.mu\
        *(P[1:,:]-P[:-1,:])/dy[:-1,:])\
        *interpolate(Cp[1:,:],Cp[:-1,:],settings['conv_interpolation'])\
        *interpolate(T[1:,:],T[:-1,:],settings['conv_interpolation'])
    
    E_kim,deta=source.Source_Comb_Kim(geom.rho_0, T, eta, 0.0)
    
    # Dimensionless numbers
#    L=settings['Carmen_diam']*settings['Porosity']/(1-settings['Porosity'])
#    L=settings['Carmen_diam']*3 # Length scale
#    Ar=source.Ea/source.R/T
    
    # Properties used for Pe and Da
    rho_avg+=np.sum(Y_0[0])/np.size(Y_0[0])
    u_avg+=np.sum(abs(v[1:,:]))/np.size(v[1:,:])
#    rho_avg=max(np.amax(Y_0[0]),rho_avg)
#    u_avg=max(np.amax(abs(v[1:,:])),u_avg)
    
    # Local quantities
#    L=settings['Length']/settings['Nodes_x']
#    Pe_x=abs(Y_0[0][:,:-1]*u[:,1:]*Cp[:,:-1]*L/k[:,:-1])
#    L=settings['Width']/settings['Nodes_y']
#    Pe_y=abs(Y_0[0][:-1,:]*v[1:,:]*Cp[:-1,:]*L/k[:-1,:])
#
#    L=settings['Carmen_diam']*3
#    Re_x=abs(Y_0[0][:,:-1]*u[:,1:]*L/settings['Darcy_mu'])
#    Re_y=abs(Y_0[0][:-1,:]*v[1:,:]*L/settings['Darcy_mu'])
    
    L=settings['Length']/settings['Nodes_x']
    t_conv_x=abs(u[:,1:]/L/sources['A0'])
#    t_conv_x=abs(u[:,1:]/L/deta[:,:-1])
    L=settings['Width']/settings['Nodes_y']
    t_conv_y=abs(v[1:,:]/L/sources['A0'])
#    t_conv_y=abs(v[1:,:]/L/deta[:-1,:])
    t_diff=abs(k/rhoC/(L)**2/sources['A0'])
#    t_diff=abs(k/rhoC/(L)**2/deta)
    
    # Plot non-dimensional numbers
#    num=[Pe_x,Pe_y,t_conv_x,t_conv_y,t_diff]
#    titles_num=['Pe_x','Pe_y','t_conv_x','t_conv_y','t_diff']
#    for i in range(len(num)):
#        fig=plt.figure(figsize=fig_size)
#        plt.contourf(X[:,:-1]*1000, Y[:,:-1]*1000, num[i], alpha=0.5, cmap=cmap_choice)#, vmin=270, vmax=2000
#        plt.contourf(X[:-1,:]*1000, Y[:-1,:]*1000, num[i], alpha=0.5, cmap=cmap_choice)#, vmin=270, vmax=2000)  
#        plt.colorbar()
#        plt.xlabel('$x$ (mm)')
#        plt.ylabel('$y$ (mm)')
#    #    plt.clim(300, 10000)
#        plt.xlim([xmin,xmax])
#        plt.ylim([ymin,ymax])
#        plt.title(titles_num[i]+' t='+time+' ms');
#        fig.savefig(titles_num[i]+'_'+time+'.png',dpi=300)
#        plt.close(fig)
    
    # Data to post-processing results file
#    fout.write('     T error: '+str((np.amin(1-T2/T),np.amax(1-T2/T)))+'\n')
#    fout.write('     Ea/R/T: '+str((np.amin(Ar),np.amax(Ar)))+'\n')
#    fout.write('     Pe_x: '+str((np.amin(Pe_x),np.amax(Pe_x)))+'\n')
#    fout.write('     Pe_y: '+str((np.amin(Pe_y),np.amax(Pe_y)))+'\n')
##    fout.write('     Re_x: '+str((np.amin(Re_x),np.amax(Re_x)))+'\n')
##    fout.write('     Re_y: '+str((np.amin(Re_y),np.amax(Re_y)))+'\n')
#    fout.write('     t_conv_x: '+str((np.amin(t_conv_x),np.amax(t_conv_x)))+'\n')
#    fout.write('     t_conv_y: '+str((np.amin(t_conv_y),np.amax(t_conv_y)))+'\n')
#    fout.write('     t_diff: '+str((np.amin(t_diff),np.amax(t_diff)))+'\n')
    
#    fout.write('     |conv/cond|: '+str((np.amin(abs(conv/cond)),np.amax(abs(conv/cond))))+'\n')
#    fout.write('     |heat gen/heat losses|: '+str((np.amin(abs(E_kim/(cond+conv))),np.amax(abs(E_kim/(cond+conv)))))+'\n')
rho_avg/=len(times)
u_avg/=len(times)
fout.write('Max pressure: %.0f\n'%(p_max))
fout.write('avg rho_g: %.4f\n'%(rho_avg))
fout.write('avg Darcy v: %.4f\n'%(u_avg))

# For non-dimensional numbers
Cp=geom.Cp_calc.get_Cp(np.ones(3)*2844, geom.Cp_g[0])[0]
#L_ref=settings['Carmen_diam']
L_ref=settings['Width']
#rho_ref=BCs['bc_right_P'][1]/settings['gas_constant']/2844
rho_ref=101325/settings['gas_constant']/2844
#rho_ref=rho_avg
#rho_ref=float(geom.rho[1])
#u_ref=geom.perm[0,0]/geom.mu*BCs['bc_right_P'][1]/L
u_ref=geom.perm[0,0]/geom.mu*101325/L_ref
#u_ref=v_BR
#u_ref=u_avg

#fout.write('Pe_y: %.12f\n'%(rho_ref*u_ref*Cp*L_ref/k[0,0]))
fout.write('Pe_y: %.12f\n'%(rho_ref*u_ref*settings['gas_constant']*L_ref/k[0,0]))

fout.write('Da_y: %f\n'%(L_ref/u_ref/(1/sources['A0'])))
#fout.write('Da_y: %f\n'%(L/v_BR/(1/sources['A0'])))

fout.write('Burn rate: %f\n'%(v_BR/u_ref))
#fout.write('Burn rate: %f\n'%(rho*v_BR*L/geom.mu))
#fout.write('Burn rate: %f\n'%(v_BR/(sources['A0']*L)))

#fout.write('Ignition delay: %f\n'%(geom.mu*t_ign/rho/L**2))
fout.write('Ignition delay: %f\n'%(t_ign*u_ref/L_ref))
#fout.write('Ignition delay: %f\n'%(t_ign*sources['A0']))
##############################################################
#               1D plots
##############################################################
if st.find(OneD_graphs,'True')>=0:
    print 'Creating 1D plots'
    fig=plt.figure(figsize=fig_size)
    for time in times:
        T=np.load('T_'+time+'.npy', False)
        # 1D temperature profile at centreline
        plt.plot(Y[:,1]*1000, T[:,int(len(T[0,:])/2)], label='t='+time)
    plt.xlabel(x_axis_labels[settings['Domain']])
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
        plt.xlabel(y_axis_labels[settings['Domain']])
        plt.ylabel('$d\eta/dt$ ($s^{-1}$)')
        plt.xlim([xmin,xmax])
        plt.legend()
        plt.title('Centreline Reaction rate Evolution')
        fig.savefig('Phi_1D.png',dpi=300)
        plt.close(fig)

fout.close()
print '\nPost-processing complete'