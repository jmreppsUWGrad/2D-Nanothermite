# -*- coding: utf-8 -*-
"""
######################################################
#             2D Heat Conduction Solver              #
#              Created by J. Mark Epps               #
#          Part of Masters Thesis at UW 2018-2020    #
######################################################

This file contains the MPI routines:
    -
    
Features:
    -Ignition condition met, will change north BC to that of right BC
    -Saves temperature and reaction data (.npy) depending on input file 
    settings

"""

import numpy as np
import string as st
import copy

class MPI_comms():
    def __init__(self, comm, rank, size, Sources, Species):
        self.comm=comm
        self.rank=rank
        self.size=size
        self.Sources=Sources
        self.Species=Species
        
    # MPI discretization routine
    def MPI_discretize(self, domain):
        # Determine process arrangement
        ranks=np.linspace(0, self.size-1, self.size).astype(int) # Array of processes
        maxDim=int(np.sqrt(self.size)) # Square root of size (max dimension)
        
        while maxDim>1:
            if self.size%maxDim!=0 or domain.Nx%maxDim!=0 or domain.Ny%(self.size/maxDim)!=0:
                maxDim-=1
            else:
                break
        if maxDim==1:
            return 1
        ranks=ranks.reshape((self.size/maxDim, maxDim))
        
        # Discretize domain to each process
        proc_x=len(ranks[0,:])
        proc_y=len(ranks[:,0])
        
        domain.Nx/=proc_x # Process discretization in x
        domain.Ny/=proc_y # Process discretization in y
        
        # Loop through each process rank in the arrangment
        for i in range(self.size/maxDim):
            # Only take parts of array if process belongs to this part of arrangment
            if self.rank in ranks[i]:
                domain.E=domain.E[(self.rank-i)*domain.Ny:(self.rank-i+1)*domain.Ny,\
                                  (self.rank-i)*domain.Nx:(self.rank-i+1)*domain.Nx]
                # domain.X, domain.Y, solver.dx, solver.dy, vol, Ax, Ay
        
        
        return 0
    # Update ghost nodes for processes
    def update_ghosts(self, domain):
        # Send to the left
        self.comm.send(len(domain.E[:,1]), dest=domain.proc_left)
        sen=domain.E[:,1].copy()
        self.comm.Send(sen, dest=domain.proc_left)
        # Receive from the right
        len_arr=len(domain.E[:,-1])
        len_arr=self.comm.recv(source=domain.proc_right)
        a=np.ones(len_arr)*domain.E[:,-1]
        self.comm.Recv(a, source=domain.proc_right)
        domain.E[:,-1]=a
        
        # Send to the right
        self.comm.send(len(domain.E[:,-2]), dest=domain.proc_right)
        sen=domain.E[:,-2].copy()
        self.comm.Send(sen, dest=domain.proc_right)
        # Receive from the left
        a=np.ones(1)*domain.E[:,0]
        self.comm.Recv(a, source=domain.proc_left)
        domain.E[:,0]=a
        
        if st.find(self.Sources['Source_Kim'],'True')>=0:
            # Send to the left, receive from the right
            a=np.ones(1)*domain.eta[-1]
            self.comm.Send(domain.eta[1], dest=domain.proc_left)
            self.comm.Recv(a, source=domain.proc_right)
            domain.eta[-1]=a
            # Send to the right, receive from the left
            a=np.ones(1)*domain.eta[0]
            self.comm.Send(domain.eta[-2], dest=domain.proc_right)
            self.comm.Recv(a, source=domain.proc_left)
            domain.eta[0]=a
        if bool(self.Species):
            # Send to the left, receive from the right
            self.comm.Send(domain.P[1], dest=domain.proc_left)
            a=np.ones(1)*domain.P[-1]
            self.comm.Recv(a, source=domain.proc_right)
            domain.P[-1]=a
            # Send to the right, receive from the left
            self.comm.Send(domain.P[-2], dest=domain.proc_right)
            a=np.ones(1)*domain.P[0]
            self.comm.Recv(a, source=domain.proc_left)
            domain.P[0]=a
            for i in self.Species['keys']:
                # Send to the left, receive from the right
                self.comm.Send(domain.m_species[i][1], dest=domain.proc_left)
                a=np.ones(1)*domain.m_species[i][-1]
                self.comm.Recv(a, source=domain.proc_right)
                domain.m_species[i][-1]=a
                # Send to the right, receive from the left
                self.comm.Send(domain.m_species[i][-2], dest=domain.proc_right)
                a=np.ones(1)*domain.m_species[i][0]
                self.comm.Recv(a, source=domain.proc_left)
                domain.m_species[i][0]=a
    
    # General function to compile a variable from all processes
    def compile_var(self, var, Domain):
        var_global=var[:-1].copy()
        # Gather at process 0
        if self.rank==0:
            for i in range(self.size-1):
                len_arr=self.comm.recv(source=i+1)
                dat=np.empty(len_arr)
                self.comm.Recv(dat, source=i+1)
                var_global=np.block([var_global, dat])
        # Any purely interior processes
        elif (Domain.proc_left>=0) and (Domain.proc_right>=0)\
            and (Domain.proc_top>=0) and (Domain.proc_bottom>=0):
            len_arr=(np.shape(var)[0]-2,np.shape(var)[1]-2)
            self.comm.send(len_arr, dest=0)
            self.comm.Send(var[1:-1,1:-1], dest=0)
        # Any corner processes
        elif (Domain.proc_left*Domain.proc_top>=0) or (Domain.proc_left*Domain.proc_bottom>=0)\
            or (Domain.proc_right*Domain.proc_top>=0) or (Domain.proc_right*Domain.proc_bottom>=0):
            len_arr=(np.shape(var)[0]-1,np.shape(var)[1]-1)
            self.comm.send(len_arr, dest=0)
            if Domain.proc_left<0:
                if Domain.proc_top<0:
                    self.comm.Send(var[:-1,:-1], dest=0)
                else:
                    self.comm.Send(var[1:,:-1], dest=0)
            else:
                if Domain.proc_top<0:
                    self.comm.Send(var[:-1,1:], dest=0)
                else:
                    self.comm.Send(var[1:,1:], dest=0)
        # Processes with BC on bottom or top only
        elif Domain.proc_bottom<0 or Domain.proc_top<0:
            len_arr=(np.shape(var)[0]-1,np.shape(var)[1]-2)
            self.comm.send(len_arr, dest=0)
            if Domain.proc_bottom<0:
                self.comm.Send(var[1:,1:-1], dest=0)
            else:
                self.comm.Send(var[:-1,1:-1], dest=0)
        # Processes with BC on right or left only
        elif Domain.proc_right<0 or Domain.proc_left<0:
            len_arr=(np.shape(var)[0]-2,np.shape(var)[1]-1)
            self.comm.send(len_arr, dest=0)
            if Domain.proc_right<0:
                self.comm.Send(var[1:-1,1:], dest=0)
            else:
                self.comm.Send(var[1:-1,:-1], dest=0)
        else:
            len_arr=len(var)-1
            self.comm.send(len_arr, dest=0)
            self.comm.Send(var[1:], dest=0)
        len_arr=self.comm.bcast(np.shape(var_global), root=0)
        if self.rank!=0:
            var_global=np.empty(len_arr)
        self.comm.Bcast(var_global, root=0)
        return var_global
    
    # Function to save data to npy files
    def save_data(self, Domain, time):
        T=self.compile_var(Domain.TempFromConserv(), Domain, self.rank)
        np.save('T_'+time, T, False)
        # Kim source term
        if st.find(self.Sources['Source_Kim'],'True')>=0:
            eta=self.compile_var(Domain.eta, Domain)
            np.save('eta_'+time, eta, False)
        if bool(self.Species):
            P=self.compile_var(Domain.P, Domain)
            np.save('P_'+time, P, False)
            for i in self.Species['keys']:
                m_i=self.compile_var(Domain.m_species[i], Domain)
                np.save('m_'+i+'_'+time, m_i, False)