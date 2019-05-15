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

class MPI_comms():
    def __init__(self, comm, rank, size, Sources, Species):
        self.comm=comm
        self.rank=rank
        self.size=size
        self.Sources=Sources
        self.Species=Species
        
    # Function to split global array to processes
    # Use for MPI_discretize and restart
    def split_var(self, var_global, domain):
        var_local=np.zeros(2)
        maxDim=len(domain.proc_arrang[0,:])
        for i in range(self.size/maxDim):
            # Only take parts of array if process belongs to this part of arrangment
            # Unique case where there is only one row of processes
            if self.rank in domain.proc_arrang[i] and maxDim==1:
                # Process at y=0
                if i==0:
                    var_local=var_global[(i)*domain.Ny:(i+1)*domain.Ny+1,:]
                    
                # Process at y=y_max
                elif i==self.size/maxDim-1:
                    var_local=var_global[(i)*domain.Ny-1:(i+1)*domain.Ny,:]
                    
                # Interior
                else:
                    var_local=var_global[(i)*domain.Ny-1:(i+1)*domain.Ny+1,:]
                    
            # Processes at y=0
            elif self.rank in domain.proc_arrang[i] and i==0:
                # Corner x=0
                if self.rank==domain.proc_arrang[i,0]:
                    var_local=var_global[(i)*domain.Ny:(i+1)*domain.Ny+1,\
                          (self.rank-domain.proc_arrang[i,0])*domain.Nx:(self.rank-domain.proc_arrang[i,0]+1)*domain.Nx+1]
                    
                    
                # Corner x=x_max
                elif self.rank==domain.proc_arrang[i,-1]:
                    var_local=var_global[(i)*domain.Ny:(i+1)*domain.Ny+1,\
                          (self.rank-domain.proc_arrang[i,0])*domain.Nx-1:(self.rank-domain.proc_arrang[i,0]+1)*domain.Nx]
                    
                    
                # Interior processes
                else:
                    var_local=var_global[(i)*domain.Ny:(i+1)*domain.Ny+1,\
                          (self.rank-domain.proc_arrang[i,0])*domain.Nx-1:(self.rank-domain.proc_arrang[i,0]+1)*domain.Nx+1]
                    
            # Processes at y=y_max
            elif self.rank in domain.proc_arrang[i] and i==self.size/maxDim-1:
                # Edge x=0
                if self.rank==domain.proc_arrang[i,0]:
                    var_local=var_global[(i)*domain.Ny-1:(i+1)*domain.Ny,\
                          (self.rank-domain.proc_arrang[i,0])*domain.Nx:(self.rank-domain.proc_arrang[i,0]+1)*domain.Nx+1]
                    
                # Edge x=x_max
                elif self.rank==domain.proc_arrang[i,-1]:
                    var_local=var_global[(i)*domain.Ny-1:(i+1)*domain.Ny,\
                          (self.rank-domain.proc_arrang[i,0])*domain.Nx-1:(self.rank-domain.proc_arrang[i,0]+1)*domain.Nx]
                    
                # Interior processes
                else:
                    var_local=var_global[(i)*domain.Ny-1:(i+1)*domain.Ny,\
                          (self.rank-domain.proc_arrang[i,0])*domain.Nx-1:(self.rank-domain.proc_arrang[i,0]+1)*domain.Nx+1]
                    
            
            # Interior y values
            elif self.rank in domain.proc_arrang[i]:
                # Edge x=0
                if self.rank==domain.proc_arrang[i,0]:
                    var_local=var_global[(i)*domain.Ny-1:(i+1)*domain.Ny+1,\
                          (self.rank-domain.proc_arrang[i,0])*domain.Nx:(self.rank-domain.proc_arrang[i,0]+1)*domain.Nx+1]
                    
                    
                # Edge x=x_max
                elif self.rank==domain.proc_arrang[i,-1]:
                    var_local=var_global[(i)*domain.Ny-1:(i+1)*domain.Ny+1,\
                          (self.rank-domain.proc_arrang[i,0])*domain.Nx-1:(self.rank-domain.proc_arrang[i,0]+1)*domain.Nx]
                    
                # Interior processes
                else:
                    var_local=var_global[(i)*domain.Ny-1:(i+1)*domain.Ny+1,\
                          (self.rank-domain.proc_arrang[i,0])*domain.Nx-1:(self.rank-domain.proc_arrang[i,0]+1)*domain.Nx+1]
        
        return var_local
    
    # MPI discretization routine
    def MPI_discretize(self, domain):
        # Determine process arrangement
        ranks=np.linspace(0, self.size-1, self.size).astype(int) # Array of processes
        maxDim=int(np.sqrt(self.size))+1 # Square root of size (max dimension)
        
        while maxDim>=1:
            if self.size%maxDim!=0 or domain.Nx%maxDim!=0 or domain.Ny%(self.size/maxDim)!=0:
                maxDim-=1
            else:
                break
        if maxDim==0:
            return 1
        ranks=ranks.reshape((self.size/maxDim, maxDim))
        domain.proc_arrang=ranks # Save process arrangment to each process domain class
        for i in range(len(domain.proc_arrang[:,0])):
                if self.rank in domain.proc_arrang[i,:]:
                    domain.proc_row=i
                    break
                else:
                    continue
        
        # Discretize domain to each process
        domain.Nx/=len(ranks[0,:])
        domain.Ny/=len(ranks[:,0])
        
        domain.X=self.split_var(domain.X, domain)
        domain.Y=self.split_var(domain.Y, domain)
        domain.dX=self.split_var(domain.dX, domain)
        domain.dY=self.split_var(domain.dY, domain)
        domain.E=self.split_var(domain.E, domain)
        
        # Designate neighboring processes
        for i in range(self.size/maxDim):
            # Unique case where there is only one row of processes
            if self.rank in ranks[i] and maxDim==1:
                domain.proc_left=-1
                domain.proc_right=-1
                # Process at y=0
                if i==0:
                    domain.proc_bottom=-1
                    domain.proc_top=ranks[i+1,0]
                    
                # Process at y=y_max
                elif i==self.size/maxDim-1:
                    domain.proc_top=-1
                    domain.proc_bottom=ranks[i-1,0]
                    
                # Interior
                else:
                    domain.proc_bottom=ranks[i-1,0]
                    domain.proc_top=ranks[i+1,0]
                    
            # Processes at y=0
            elif self.rank in ranks[i] and i==0:
                domain.proc_bottom=-1
                # Corner x=0
                if self.rank==ranks[i,0]:
                    domain.proc_top=ranks[i+1,0]
                    domain.proc_left=-1
                    domain.proc_right=ranks[i,1]
                    
                    
                # Corner x=x_max
                elif self.rank==ranks[i,-1]:
                    domain.proc_top=ranks[i+1,-1]
                    domain.proc_left=ranks[i,-2]
                    domain.proc_right=-1
                    
                    
                # Interior processes
                else:
                    domain.proc_top=ranks[i+1,self.rank-ranks[i,0]]
                    domain.proc_left=ranks[i,self.rank-ranks[i,0]-1]
                    domain.proc_right=ranks[i,self.rank-ranks[i,0]+1]
                    
            # Processes at y=y_max
            elif self.rank in ranks[i] and i==self.size/maxDim-1:
                domain.proc_top=-1
                # Edge x=0
                if self.rank==ranks[i,0]:
                    domain.proc_left=-1
                    domain.proc_right=ranks[i,1]
                    domain.proc_bottom=ranks[i-1,0]
                    
                # Edge x=x_max
                elif self.rank==ranks[i,-1]:
                    domain.proc_left=ranks[i,-2]
                    domain.proc_right=-1
                    domain.proc_bottom=ranks[i-1,-1]
                    
                # Interior processes
                else:
                    domain.proc_left=ranks[i,self.rank-ranks[i,0]-1]
                    domain.proc_right=ranks[i,self.rank-ranks[i,0]+1]
                    domain.proc_bottom=ranks[i-1,self.rank-ranks[i,0]]
            
            # Interior y values
            elif self.rank in ranks[i]:
                domain.proc_top=ranks[i+1,self.rank-ranks[i,0]]
                domain.proc_bottom=ranks[i-1,self.rank-ranks[i,0]]
                # Edge x=0
                if self.rank==ranks[i,0]:
                    domain.proc_left=-1
                    domain.proc_right=ranks[i,self.rank-ranks[i,0]+1]
                    
                    
                # Edge x=x_max
                elif self.rank==ranks[i,-1]:
                    domain.proc_left=ranks[i,self.rank-ranks[i,0]-1]
                    domain.proc_right=-1
                    
                # Interior processes
                else:
                    domain.proc_left=ranks[i,self.rank-ranks[i,0]-1]
                    domain.proc_right=ranks[i,self.rank-ranks[i,0]+1]
            
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
        len_arr=len(domain.E[:,0])
        len_arr=self.comm.recv(source=domain.proc_left)
        a=np.ones(len_arr)*domain.E[:,0]
        self.comm.Recv(a, source=domain.proc_left)
        domain.E[:,0]=a
        
        # Send to the bottom
        self.comm.send(len(domain.E[1,:]), dest=domain.proc_bottom)
        sen=domain.E[1,:].copy()
        self.comm.Send(sen, dest=domain.proc_bottom)
        # Receive from the top
        len_arr=len(domain.E[-1,:])
        len_arr=self.comm.recv(source=domain.proc_top)
        a=np.ones(len_arr)*domain.E[-1,:]
        self.comm.Recv(a, source=domain.proc_top)
        domain.E[-1,:]=a
        
        # Send to the top
        self.comm.send(len(domain.E[-2,:]), dest=domain.proc_top)
        sen=domain.E[-2,:].copy()
        self.comm.Send(sen, dest=domain.proc_top)
        # Receive from the bottom
        len_arr=len(domain.E[0,:])
        len_arr=self.comm.recv(source=domain.proc_bottom)
        a=np.ones(len_arr)*domain.E[0,:]
        self.comm.Recv(a, source=domain.proc_bottom)
        domain.E[0,:]=a
        
        if st.find(self.Sources['Source_Kim'],'True')>=0:
            # Send to the left
            self.comm.send(len(domain.eta[:,1]), dest=domain.proc_left)
            sen=domain.eta[:,1].copy()
            self.comm.Send(sen, dest=domain.proc_left)
            # Receive from the right
            len_arr=len(domain.eta[:,-1])
            len_arr=self.comm.recv(source=domain.proc_right)
            a=np.ones(len_arr)*domain.eta[:,-1]
            self.comm.Recv(a, source=domain.proc_right)
            domain.eta[:,-1]=a
            
            # Send to the right
            self.comm.send(len(domain.eta[:,-2]), dest=domain.proc_right)
            sen=domain.eta[:,-2].copy()
            self.comm.Send(sen, dest=domain.proc_right)
            # Receive from the left
            len_arr=len(domain.eta[:,0])
            len_arr=self.comm.recv(source=domain.proc_left)
            a=np.ones(len_arr)*domain.eta[:,0]
            self.comm.Recv(a, source=domain.proc_left)
            domain.eta[:,0]=a
            
            # Send to the bottom
            self.comm.send(len(domain.eta[1,:]), dest=domain.proc_bottom)
            sen=domain.eta[1,:].copy()
            self.comm.Send(sen, dest=domain.proc_bottom)
            # Receive from the top
            len_arr=len(domain.eta[-1,:])
            len_arr=self.comm.recv(source=domain.proc_top)
            a=np.ones(len_arr)*domain.eta[-1,:]
            self.comm.Recv(a, source=domain.proc_top)
            domain.eta[-1,:]=a
            
            # Send to the top
            self.comm.send(len(domain.eta[-2,:]), dest=domain.proc_top)
            sen=domain.eta[-2,:].copy()
            self.comm.Send(sen, dest=domain.proc_top)
            # Receive from the bottom
            len_arr=len(domain.eta[0,:])
            len_arr=self.comm.recv(source=domain.proc_bottom)
            a=np.ones(len_arr)*domain.eta[0,:]
            self.comm.Recv(a, source=domain.proc_bottom)
            domain.eta[0,:]=a
        
        if bool(self.Species):
            # Send to the left
            self.comm.send(len(domain.P[:,1]), dest=domain.proc_left)
            sen=domain.P[:,1].copy()
            self.comm.Send(sen, dest=domain.proc_left)
            # Receive from the right
            len_arr=len(domain.P[:,-1])
            len_arr=self.comm.recv(source=domain.proc_right)
            a=np.ones(len_arr)*domain.P[:,-1]
            self.comm.Recv(a, source=domain.proc_right)
            domain.P[:,-1]=a
            
            # Send to the right
            self.comm.send(len(domain.P[:,-2]), dest=domain.proc_right)
            sen=domain.P[:,-2].copy()
            self.comm.Send(sen, dest=domain.proc_right)
            # Receive from the left
            len_arr=len(domain.P[:,0])
            len_arr=self.comm.recv(source=domain.proc_left)
            a=np.ones(len_arr)*domain.P[:,0]
            self.comm.Recv(a, source=domain.proc_left)
            domain.P[:,0]=a
            
            # Send to the bottom
            self.comm.send(len(domain.P[1,:]), dest=domain.proc_bottom)
            sen=domain.P[1,:].copy()
            self.comm.Send(sen, dest=domain.proc_bottom)
            # Receive from the top
            len_arr=len(domain.P[-1,:])
            len_arr=self.comm.recv(source=domain.proc_top)
            a=np.ones(len_arr)*domain.P[-1,:]
            self.comm.Recv(a, source=domain.proc_top)
            domain.P[-1,:]=a
            
            # Send to the top
            self.comm.send(len(domain.P[-2,:]), dest=domain.proc_top)
            sen=domain.P[-2,:].copy()
            self.comm.Send(sen, dest=domain.proc_top)
            # Receive from the bottom
            len_arr=len(domain.P[0,:])
            len_arr=self.comm.recv(source=domain.proc_bottom)
            a=np.ones(len_arr)*domain.P[0,:]
            self.comm.Recv(a, source=domain.proc_bottom)
            domain.P[0,:]=a
            
            for i in self.Species['keys']:
                # Send to the left
                self.comm.send(len(domain.rho_species[i][:,1]), dest=domain.proc_left)
                sen=domain.rho_species[i][:,1].copy()
                self.comm.Send(sen, dest=domain.proc_left)
                # Receive from the right
                len_arr=len(domain.rho_species[i][:,-1])
                len_arr=self.comm.recv(source=domain.proc_right)
                a=np.ones(len_arr)*domain.rho_species[i][:,-1]
                self.comm.Recv(a, source=domain.proc_right)
                domain.rho_species[i][:,-1]=a
                
                # Send to the right
                self.comm.send(len(domain.rho_species[i][:,-2]), dest=domain.proc_right)
                sen=domain.rho_species[i][:,-2].copy()
                self.comm.Send(sen, dest=domain.proc_right)
                # Receive from the left
                len_arr=len(domain.rho_species[i][:,0])
                len_arr=self.comm.recv(source=domain.proc_left)
                a=np.ones(len_arr)*domain.rho_species[i][:,0]
                self.comm.Recv(a, source=domain.proc_left)
                domain.rho_species[i][:,0]=a
                
                # Send to the bottom
                self.comm.send(len(domain.rho_species[i][1,:]), dest=domain.proc_bottom)
                sen=domain.rho_species[i][1,:].copy()
                self.comm.Send(sen, dest=domain.proc_bottom)
                # Receive from the top
                len_arr=len(domain.rho_species[i][-1,:])
                len_arr=self.comm.recv(source=domain.proc_top)
                a=np.ones(len_arr)*domain.rho_species[i][-1,:]
                self.comm.Recv(a, source=domain.proc_top)
                domain.rho_species[i][-1,:]=a
                
                # Send to the top
                self.comm.send(len(domain.rho_species[i][-2,:]), dest=domain.proc_top)
                sen=domain.rho_species[i][-2,:].copy()
                self.comm.Send(sen, dest=domain.proc_top)
                # Receive from the bottom
                len_arr=len(domain.rho_species[i][0,:])
                len_arr=self.comm.recv(source=domain.proc_bottom)
                a=np.ones(len_arr)*domain.rho_species[i][0,:]
                self.comm.Recv(a, source=domain.proc_bottom)
                domain.rho_species[i][0,:]=a
                
    # General function to compile a variable from all processes
    def compile_var(self, var, Domain):
        var_global=var[1:-1,1:-1].copy()
        # For one column of processes
        if Domain.proc_arrang.shape[1]==1:
            # Start global variable based on which row process is in
            if Domain.proc_row==0:
                var_global=var[:-1,:].copy()
            elif self.rank==Domain.proc_arrang[-1,0]:
                var_global=var[1:,:].copy()
            else:
                var_global=var[1:-1,:].copy()
            # Send variable to process 0 for final compile
            if self.rank!=0:
                self.comm.send(np.shape(var_global), dest=0)
                self.comm.Send(var_global, dest=0)
            # Process 0 compiles global variable
            else:
                for i in range(len(Domain.proc_arrang[:,0])-1):
                    len_arr=self.comm.recv(source=Domain.proc_arrang[i+1,0])
                    a=np.empty(len_arr)
                    self.comm.Recv(a, source=Domain.proc_arrang[i+1,0])
                    var_global=np.block([[var_global],[a]])
        # Collect data from each row of processes at first process of the row
        elif self.rank in Domain.proc_arrang[:,0]:
            # Start global variable based on which row process is in
            if Domain.proc_row==0:
                var_global=var[:-1,:-1].copy()
            elif self.rank==Domain.proc_arrang[-1,0]:
                var_global=var[1:,:-1].copy()
            else:
                var_global=var[1:-1,:-1].copy()
            # Cycle through each process in x (right)
            for i in range(len(Domain.proc_arrang[0,:])-1):
                len_arr=self.comm.recv(source=Domain.proc_arrang[Domain.proc_row,i+1])
                a=np.empty(len_arr)
                self.comm.Recv(a, source=Domain.proc_arrang[Domain.proc_row,i+1])
                var_global=np.block([var_global, a])
            
            # Send compiled variable to process 0 for final compile
            if self.rank!=0:
                self.comm.send(np.shape(var_global), dest=0)
                self.comm.Send(var_global, dest=0)
            # Process 0 compiles global variable
            else:
                for i in range(len(Domain.proc_arrang[:,0])-1):
                    len_arr=self.comm.recv(source=Domain.proc_arrang[i+1,0])
                    a=np.empty(len_arr)
                    self.comm.Recv(a, source=Domain.proc_arrang[i+1,0])
                    var_global=np.block([[var_global],[a]])
        
        # Processes in each row will send data to first process of row
        else:# self.rank in Domain.proc_arrang[:,1:]:
            # First row (y=0)
            if Domain.proc_row==0:
                # Corner (x=x_max)
                if self.rank==Domain.proc_arrang[Domain.proc_row,-1]:
                    sen=var[:-1,1:].copy()
                # Interior processes
                else:
                    sen=var[:-1,1:-1].copy()
            # Last row (y=y_max)
            elif self.rank in Domain.proc_arrang[-1,:]:
                # Corner (x=x_max)
                if self.rank==Domain.proc_arrang[Domain.proc_row,-1]:
                    sen=var[1:,1:].copy()
                # Interior processes
                else:
                    sen=var[1:,1:-1].copy()
            # Interior row
            else:
                # Edge (x=x_max)
                if self.rank==Domain.proc_arrang[Domain.proc_row,-1]:
                    sen=var[1:-1,1:].copy()
                # Interior processes
                else:
                    sen=var[1:-1,1:-1].copy()
            self.comm.send(np.shape(sen), dest=Domain.proc_arrang[Domain.proc_row,0])
            self.comm.Send(sen, dest=Domain.proc_arrang[Domain.proc_row,0])
            
        len_arr=self.comm.bcast(np.shape(var_global), root=0)
        if self.rank!=0:
            var_global=np.empty(len_arr)
        self.comm.Bcast(var_global, root=0)
        
        return var_global
    
    # Function to save data to npy files
    def save_data(self, Domain, Sources, Species, time):
        T=self.compile_var(Domain.TempFromConserv(), Domain)
        np.save('T_'+time, T, False)
        # Kim source term
        if st.find(self.Sources['Source_Kim'],'True')>=0:
            eta=self.compile_var(Domain.eta, Domain)
            np.save('eta_'+time, eta, False)
        if bool(self.Species):
            P=self.compile_var(Domain.P, Domain)
            np.save('P_'+time, P, False)
            for i in self.Species['keys']:
                m_i=self.compile_var(Domain.rho_species[i], Domain)
                np.save('rho_'+i+'_'+time, m_i, False)