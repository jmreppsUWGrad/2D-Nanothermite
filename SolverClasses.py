# -*- coding: utf-8 -*-
"""
######################################################
#             2D Heat Conduction Solver              #
#              Created by J. Mark Epps               #
#          Part of Masters Thesis at UW 2018-2020    #
######################################################

This file contains the solver classes for 2D planar and axisymmetric 
Heat Conduction:
    -Modifies geometry class variables directly
    -Calculate time step based on Fourier number
    -Compute conduction equations
    -Add source terms as needed
    -Applies boundary conditions (can vary along a side)

Features/assumptions:
    -time step based on Fourrier number and local discretizations in x and y
    -equal node spacing in x or y
    -thermal properties can vary in space (call from geometry object)
    -Radiation boundary conditions

"""

import numpy as np
import copy
#import CoolProp.CoolProp as CP
#import temporal_schemes
import Source_Comb
import BCClasses
from mpi4py import MPI

# 2D solver (Cartesian coordinates)
class TwoDimSolver():
    def __init__(self, geom_obj, settings, Sources, BCs, comm):
        self.Domain=geom_obj # Geometry object
        self.time_scheme=settings['Time_Scheme']
        self.dx,self.dy=geom_obj.dX,geom_obj.dY
        self.Fo=settings['Fo']
        self.CFL=settings['CFL']
        self.dt=settings['dt']
        self.conv=settings['Convergence']
        self.countmax=settings['Max_iterations']
        self.comm=comm
        
        # Define source terms and pointer to source object here
        self.get_source=Source_Comb.Source_terms(Sources['Ea'], Sources['A0'], Sources['dH'])
        self.source_unif=Sources['Source_Uniform']
        self.source_Kim=Sources['Source_Kim']
        
        # BC class
        self.BCs=BCClasses.BCs(BCs, self.dx, self.dy, settings['Domain'])
        # Ensure proper BCs for this process
        self.mult_BCs(BCs)
    
    # Modify BCs based on processes next to current one AND if multiple
    # BCs are specified on a given boundary
    def mult_BCs(self, BC_global):
        # Left boundary
        if self.Domain.proc_left>=0:
            self.BCs.BCs['bc_left_E']=['F', 0.0, (0, -1)]
        # Global boundary with multiple BCs
        elif len(BC_global['bc_left_E'])>3:
            i=len(BC_global['bc_left_E'])/3
            j=0
            while i>j:
                # Lower bound of BC in this process
                if BC_global['bc_left_E'][2+3*j][0]>=self.Domain.proc_row*self.Domain.Ny\
                    and BC_global['bc_left_E'][2+3*j][0]<(self.Domain.proc_row+1)*self.Domain.Ny:
                    st=BC_global['bc_left_E'][2+3*j][0]-self.Domain.proc_row*self.Domain.Ny
                    # upper bound of BC in this process
                    if BC_global['bc_left_E'][2+3*j][1]<=(self.Domain.proc_row+1)*self.Domain.Ny:
                        en=BC_global['bc_left_E'][2+3*j][1]-self.Domain.proc_row*self.Domain.Ny
                    # upper bound outside this process
                    else:
                        en=self.Domain.Ny
                    # Ghost node on bottom
                    if self.Domain.proc_bottom>=0:
                        st+=1
                        en+=1
                    elif self.Domain.proc_top<0:
                        en+=1
                    self.BCs.BCs['bc_left_E'][2+3*j]=(st,en)
                    j+=1
                # Lower bound of BC not in this process, but upper bound is
                elif BC_global['bc_left_E'][2+3*j][1]<=(self.Domain.proc_row+1)*self.Domain.Ny\
                    and BC_global['bc_left_E'][2+3*j][1]>self.Domain.proc_row*self.Domain.Ny:
                    st=0
                    en=BC_global['bc_left_E'][2+3*j][1]-self.Domain.proc_row*self.Domain.Ny
                    # Ghost node on bottom
                    if self.Domain.proc_bottom>=0:
                        st+=1
                        en+=1
                    elif self.Domain.proc_top<0:
                        en+=1
                    self.BCs.BCs['bc_left_E'][2+3*j]=(st,en)
                    j+=1
                    
                # Process lies inside the upper and lower bounds (are outside process)
                elif BC_global['bc_left_E'][2+3*j][0]<self.Domain.proc_row*self.Domain.Ny\
                    and BC_global['bc_left_E'][2+3*j][1]>(self.Domain.proc_row+1)*self.Domain.Ny:
                    self.BCs.BCs['bc_left_E'][2+3*j]=(0,-1)
                    j+=1
                # BC has no effect on this process
                else:
                    del self.BCs.BCs['bc_left_E'][3*j:3+3*j]
                    i-=1

        # Right boundary
        if self.Domain.proc_right>=0:
            self.BCs.BCs['bc_right_E']=['F', 0.0, (0, -1)]
        # Global boundary with multiple BCs
        elif len(BC_global['bc_right_E'])>3:
            i=len(BC_global['bc_right_E'])/3
            j=0
            while i>j:
                # Lower bound of BC in this process
                if BC_global['bc_right_E'][2+3*j][0]>=self.Domain.proc_row*self.Domain.Ny\
                    and BC_global['bc_right_E'][2+3*j][0]<(self.Domain.proc_row+1)*self.Domain.Ny:
                    st=BC_global['bc_right_E'][2+3*j][0]-self.Domain.proc_row*self.Domain.Ny
                    # upper bound of BC in this process
                    if BC_global['bc_right_E'][2+3*j][1]<=(self.Domain.proc_row+1)*self.Domain.Ny:
                        en=BC_global['bc_right_E'][2+3*j][1]-self.Domain.proc_row*self.Domain.Ny
                    # upper bound outside this process
                    else:
                        en=self.Domain.Ny
                    # Ghost node on bottom
                    if self.Domain.proc_bottom>=0:
                        st+=1
                        en+=1
                    elif self.Domain.proc_top<0:
                        en+=1
                    self.BCs.BCs['bc_right_E'][2+3*j]=(st,en)
                    j+=1
                # Lower bound of BC not in this process, but upper bound is
                elif BC_global['bc_right_E'][2+3*j][1]<=(self.Domain.proc_row+1)*self.Domain.Ny\
                    and BC_global['bc_right_E'][2+3*j][1]>self.Domain.proc_row*self.Domain.Ny:
                    st=0
                    en=BC_global['bc_right_E'][2+3*j][1]-self.Domain.proc_row*self.Domain.Ny
                    # Ghost node on bottom
                    if self.Domain.proc_bottom>=0:
                        st+=1
                        en+=1
                    elif self.Domain.proc_top<0:
                        en+=1
                    self.BCs.BCs['bc_right_E'][2+3*j]=(st,en)
                    j+=1
                
                # Process lies inside the upper and lower bounds (are outside process)
                elif BC_global['bc_right_E'][2+3*j][0]<self.Domain.proc_row*self.Domain.Ny\
                    and BC_global['bc_right_E'][2+3*j][1]>(self.Domain.proc_row+1)*self.Domain.Ny:
                    self.BCs.BCs['bc_right_E'][2+3*j]=(0,-1)
                    j+=1
                # BC has no effect on this process
                else:
                    del self.BCs.BCs['bc_right_E'][3*j:3+3*j]
                    i-=1
        
        # Locate column of current process
        for i in range(len(self.Domain.proc_arrang[0,:])):
            if self.Domain.rank in self.Domain.proc_arrang[:,i]:
                coln=i
                break
            else:
                continue
        
        # Top boundary
        if self.Domain.proc_top>=0:
            self.BCs.BCs['bc_north_E']=['F', 0.0, (0, -1)]
        # Global boundary with multiple BCs
        elif len(BC_global['bc_north_E'])>3:
            i=len(BC_global['bc_north_E'])/3
            j=0
            while i>j:
                # Lower bound of BC in this process
                if BC_global['bc_north_E'][2+3*j][0]>=coln*self.Domain.Nx\
                    and BC_global['bc_north_E'][2+3*j][0]<(coln+1)*self.Domain.Nx:
                    st=BC_global['bc_north_E'][2+3*j][0]-coln*self.Domain.Nx
                    # upper bound of BC in this process
                    if BC_global['bc_north_E'][2+3*j][1]<=(coln+1)*self.Domain.Nx:
                        en=BC_global['bc_north_E'][2+3*j][1]-coln*self.Domain.Nx
                    # upper bound outside this process
                    else:
                        en=self.Domain.Nx
                    # Ghost node on left
                    if self.Domain.proc_left>=0:
                        st+=1
                        en+=1
                    elif self.Domain.proc_right<0:
                        en+=1
                    self.BCs.BCs['bc_north_E'][2+3*j]=(st,en)
                    j+=1
                # Lower bound of BC not in this process, but upper bound is
                elif BC_global['bc_north_E'][2+3*j][1]<=(coln+1)*self.Domain.Nx\
                    and BC_global['bc_north_E'][2+3*j][1]>coln*self.Domain.Nx:
                    st=0
                    en=BC_global['bc_north_E'][2+3*j][1]-coln*self.Domain.Nx
                    # Ghost node on left
                    if self.Domain.proc_left>=0:
                        st+=1
                        en+=1
                    elif self.Domain.proc_right<0:
                        en+=1
                    self.BCs.BCs['bc_north_E'][2+3*j]=(st,en)
                    j+=1
                
                # Process lies inside the upper and lower bounds (are outside process)
                elif BC_global['bc_north_E'][2+3*j][0]<coln*self.Domain.Nx\
                    and BC_global['bc_north_E'][2+3*j][1]>(coln+1)*self.Domain.Nx:
                    self.BCs.BCs['bc_north_E'][2+3*j]=(0,-1)
                    j+=1
                # BC has no effect on this process
                else:
                    del self.BCs.BCs['bc_north_E'][3*j:3+3*j]
                    i-=1
        
        # Bottom boundary
        if self.Domain.proc_bottom>=0:
            self.BCs.BCs['bc_south_E']=['F', 0.0, (0, -1)]
        # Global boundary with multiple BCs
        elif len(BC_global['bc_south_E'])>3:
            i=len(BC_global['bc_south_E'])/3
            j=0
            while i>j:
                # Lower bound of BC in this process
                if BC_global['bc_south_E'][2+3*j][0]>=coln*self.Domain.Nx\
                    and BC_global['bc_south_E'][2+3*j][0]<(coln+1)*self.Domain.Nx:
                    st=BC_global['bc_south_E'][2+3*j][0]-coln*self.Domain.Nx
                    # upper bound of BC in this process
                    if BC_global['bc_south_E'][2+3*j][1]<=(coln+1)*self.Domain.Nx:
                        en=BC_global['bc_south_E'][2+3*j][1]-coln*self.Domain.Nx
                    # upper bound outside this process
                    else:
                        en=self.Domain.Nx
                    # Ghost node on left
                    if self.Domain.proc_left>=0:
                        st+=1
                        en+=1
                    elif self.Domain.proc_right<0:
                        en+=1
                    self.BCs.BCs['bc_south_E'][2+3*j]=(st,en)
                    j+=1
                # Lower bound of BC not in this process, but upper bound is
                elif BC_global['bc_south_E'][2+3*j][1]<=(coln+1)*self.Domain.Nx\
                    and BC_global['bc_south_E'][2+3*j][1]>coln*self.Domain.Nx:
                    st=0
                    en=BC_global['bc_south_E'][2+3*j][1]-coln*self.Domain.Nx
                    # Ghost node on left
                    if self.Domain.proc_left>=0:
                        st+=1
                        en+=1
                    elif self.Domain.proc_right<0:
                        en+=1
                    self.BCs.BCs['bc_south_E'][2+3*j]=(st,en)
                    j+=1
                
                # Process lies inside the upper and lower bounds (are outside process)
                elif BC_global['bc_south_E'][2+3*j][0]<coln*self.Domain.Nx\
                    and BC_global['bc_south_E'][2+3*j][1]>(coln+1)*self.Domain.Nx:
                    self.BCs.BCs['bc_south_E'][2+3*j]=(0,-1)
                    j+=1
                # BC has no effect on this process
                else:
                    del self.BCs.BCs['bc_south_E'][3*j:3+3*j]
                    i-=1
        
    # Time step check with dx, dy, Fo number
    def getdt(self, k, rho, Cv, T):
        # Time steps depending on Fo
        dt_1=np.amin(self.Fo*rho*Cv/k*(self.dx)**2)
        dt_2=np.amin(self.Fo*rho*Cv/k*(self.dy)**2)
        
        # Time steps depending on CFL (if flow model used)
        if bool(self.Domain.rho_species):
            dt_3=np.amin(self.CFL*self.dx/np.sqrt(1.4*8.314/102*T))
            dt_4=np.amin(self.CFL*self.dy/np.sqrt(1.4*8.314/102*T))
        else:
            dt_3=10**9
            dt_4=10**9
        
        return min(dt_1,dt_2,dt_3,dt_4)
    
    # Interpolation function
    def interpolate(self, k1, k2, func):
        if func=='Linear':
            return 0.5*k1+0.5*k2
        else:
            return 2*k1*k2/(k1+k2)
        
    # Main solver (1 time step)
    def Advance_Soln_Cond(self, nt, t, hx, hy):
        max_Y,min_Y=0,1
        # Calculate properties
        k, rho, Cv, Cp, D=self.Domain.calcProp()
        
        # Copy needed variables and set pointers to other variables
        T_c=self.Domain.TempFromConserv()
        mu=self.Domain.mu
        perm=self.Domain.perm
        if bool(self.Domain.rho_species):
            rho_spec=copy.deepcopy(self.Domain.rho_species)
            species=self.Domain.species_keys
            Cp_spec=self.Domain.Cp_species
#            mu_c=copy.deepcopy(self.Domain.mu_species)
#            mv_c=copy.deepcopy(self.Domain.mv_species)
            
            # Velocity
#            u=mu_c[species[0]]/m_c[species[0]]
#            v=mv_c[species[0]]/m_c[species[0]]
        
        if self.dt=='None':
            dt=self.getdt(k, rho, Cv, T_c)
            # Collect all dt from other processes and send minimum
            dt=self.comm.reduce(dt, op=MPI.MIN, root=0)
            dt=self.comm.bcast(dt, root=0)
        else:
            dt=min(self.dt,self.getdt(k, rho, Cv, T_c))
            # Collect all dt from other processes and send minimum
            dt=self.comm.reduce(dt, op=MPI.MIN, root=0)
            dt=self.comm.bcast(dt, root=0)
        if (np.isnan(dt)) or (dt<=0):
            return 1, dt
        if self.Domain.rank==0:
            print 'Time step %i, Step size=%.7f, Time elapsed=%f;'%(nt+1,dt, t+dt)
        
        ###################################################################
        # Calculate source and Porous medium terms
        ###################################################################
        # Source terms
        E_unif,E_kim=0,0
        if self.source_unif!='None':
            E_unif      = self.source_unif
        if self.source_Kim=='True':
#            self.Domain.eta=self.Domain.rho_species[:,:,2]/0.25
            E_kim, deta =self.get_source.Source_Comb_Kim(rho, T_c, self.Domain.eta, dt)
#            E_kim, deta =self.get_source.Source_Comb_Umbrajkar(rho, T_c, self.Domain.eta, dt)
        
        ###################################################################
        # Conservation of Mass
        ###################################################################
        if bool(self.Domain.rho_species):
            # Adjust pressure
#            print '     Gas mass: %f, %f'%(np.amax(self.Domain.rho_species['g'])*10**6,np.amin(self.Domain.rho_species['g'])*10**6)
#            print '     Gas density: %f, %f'%(np.amax(rho_spec['g']),np.amin(rho_spec['g']))
            self.Domain.P=self.Domain.rho_species['g']*self.Domain.R*T_c
    #        self.BCs.P(self.Domain.P)
#            print '     Pressure: %f, %f'%(np.amax(self.Domain.P),np.amin(self.Domain.P))
            
            # Use Darcy's law to directly calculate the velocities at the faces
            # Ingoing fluxes
            flx=np.zeros_like(self.Domain.P)
            fly=np.zeros_like(self.Domain.P)
            
            # Left face
            flx[:,1:]+=dt/hx[:,1:]\
                *self.interpolate(rho_spec[species[0]][:,1:],rho_spec[species[0]][:,:-1],'Linear')*\
                (-perm/mu*(self.Domain.P[:,1:]-self.Domain.P[:,:-1])/self.dx[:,:-1])
    #            self.interpolate(u[:,1:], u[:,:-1], 'Linear')
            
            # Right face
            flx[:,:-1]-=dt/hx[:,:-1]\
                *self.interpolate(rho_spec[species[0]][:,1:],rho_spec[species[0]][:,:-1], 'Linear')*\
                (-perm/mu*(self.Domain.P[:,1:]-self.Domain.P[:,:-1])/self.dx[:,:-1])
    #            self.interpolate(u[:,1:], u[:,:-1], 'Linear')
            
            # South face
            fly[1:,:]+=dt/hy[1:,:]\
                *self.interpolate(rho_spec[species[0]][1:,:],rho_spec[species[0]][:-1,:],'Linear')*\
                (-perm/mu*(self.Domain.P[1:,:]-self.Domain.P[:-1,:])/self.dy[:-1,:])
    #            self.interpolate(v[1:,:], v[:-1,:], 'Linear')
            
            # North face
            fly[:-1,:]-=dt/hy[:-1,:]\
                *self.interpolate(rho_spec[species[0]][1:,:], rho_spec[species[0]][:-1,:], 'Linear')*\
                (-perm/mu*(self.Domain.P[1:,:]-self.Domain.P[:-1,:])/self.dy[:-1,:])
    #            self.interpolate(v[1:,:], v[:-1,:], 'Linear')
            
#            print '    Gas fluxes in x: %f, %f'%(np.amax(flx)*10**(9),np.amin(flx)*10**(9))
#            print '    Gas fluxes in y: %f, %f'%(np.amax(fly)*10**(9),np.amin(fly)*10**(9))
            
            self.Domain.rho_species[species[0]]+=flx+fly
            
            # Source terms
    #        dm=deta*dt*(m_c[species[0]]+m_c[species[1]])
#            dm=np.zeros_like(deta)
            dm=deta*dt*(self.Domain.rho_0)
#            dm[dm<10**(-9)]=0
#            print '     Mass generated: %f, %f'%(np.amax(dm)*10**(9),np.amin(dm)*10**(9))
    #        (m_c[species[0]]+m_c[species[1]])
            self.Domain.rho_species[species[0]]+=dm/self.Domain.porosity
            self.Domain.rho_species[species[1]]-=dm/(1-self.Domain.porosity)
                    
            max_Y=max(np.amax(self.Domain.rho_species[species[0]]),\
                      np.amax(self.Domain.rho_species[species[1]]))
            min_Y=min(np.amin(self.Domain.rho_species[species[0]]),\
                      np.amin(self.Domain.rho_species[species[1]]))
            
            # Apply BCs
#            self.BCs.mass(self.Domain.rho_species[species[0]], self.Domain.P, Ax, Ay, vol)
        
        ###################################################################
        # Conservation of Momentum (x direction; gas)
        ###################################################################
        # Fluxes
#        self.Domain.mu_species[species[0]][:,1:]+=Ax[:,1:]*dt\
#            *0.5*(mu_c[species[0]][:,1:]+mu_c[species[0]][:,:-1])*self.interpolate(u[:,1:], u[:,:-1], 'Linear')
#        self.Domain.mu_species[species[0]][1:,:]+=Ay[1:,:]*dt\
#            *0.5*(mu_c[species[0]][1:,:]+mu_c[species[0]][:-1,:])*self.interpolate(v[1:,:], v[:-1,:], 'Linear')
#        self.Domain.mu_species[species[0]][:,:-1]-=Ax[:,:-1]*dt\
#            *0.5*(mu_c[species[0]][:,1:]+mu_c[species[0]][:,:-1])*self.interpolate(u[:,1:], u[:,:-1], 'Linear')
#        self.Domain.mu_species[species[0]][:-1,:]-=Ay[:-1,:]*dt\
#            *0.5*(mu_c[species[0]][1:,:]+mu_c[species[0]][:-1,:])*self.interpolate(v[1:,:], v[:-1,:], 'Linear')
#        
#        # Pressure
#        self.Domain.mu_species[species[0]][:,1:]+=Ax[:,1:]*dt\
#            *0.5*(self.Domain.P[:,1:]+self.Domain.P[:,:-1])
#        self.Domain.mu_species[species[0]][:,:-1]-=Ax[:,:-1]*dt\
#            *0.5*(self.Domain.P[:,1:]+self.Domain.P[:,:-1])
#                
#        # Porous medium losses
#        self.Domain.mu_species[species[0]]-=mu/perm*u*vol*dt
#        
#        ###################################################################
#        # Conservation of Momentum (y direction; gas)
#        ###################################################################
#        # Fluxes
#        self.Domain.mv_species[species[0]][:,1:]+=Ax[:,1:]*dt\
#            *0.5*(mv_c[species[0]][:,1:]+mv_c[species[0]][:,:-1])*0.5*(u[:,1:]+u[:,:-1])
#        self.Domain.mv_species[species[0]][1:,:]+=Ay[1:,:]*dt\
#            *0.5*(mv_c[species[0]][1:,:]+mv_c[species[0]][:-1,:])*0.5*(v[1:,:]+v[:-1,:])
#        self.Domain.mv_species[species[0]][:,:-1]-=Ax[:,:-1]*dt\
#            *0.5*(mv_c[species[0]][:,1:]+mv_c[species[0]][:,:-1])*0.5*(u[:,1:]+u[:,:-1])
#        self.Domain.mv_species[species[0]][:-1,:]-=Ay[:-1,:]*dt\
#            *0.5*(mv_c[species[0]][1:,:]+mv_c[species[0]][:-1,:])*0.5*(v[1:,:]+v[:-1,:])
#        
#        # Pressure
#        self.Domain.mv_species[species[0]][1:,:]+=Ay[1:,:]*dt\
#            *0.5*(self.Domain.P[1:,:]+self.Domain.P[:-1,:])
#        self.Domain.mv_species[species[0]][:-1,:]-=Ay[:-1,:]*dt\
#            *0.5*(self.Domain.P[1:,:]+self.Domain.P[:-1,:])
#                
#        # Porous medium losses
#        self.Domain.mu_species[species[0]]-=mu/perm*v*vol*dt
        
        
        ###################################################################
        # Conservation of species
        ###################################################################
#        if bool(self.Domain.m_species):
#            # Mole ratios
#            mole_ratio={}
#            mole_ratio[species[0]]=-2.0/5 # Al
#            mole_ratio[species[1]]=-3.0/5 # CuO
#            mole_ratio[species[2]]=1.0/4  # Al2O3
#            mole_ratio[species[3]]=3.0/4  # Cu
#            
#            for i in species:
#                # Calculate flux coefficients
#                aW,aE,aS,aN=self.get_Coeff(self.dx,self.dy, dt, rho*D[i], 'Linear')
#                
#                # Diffusion contribution (2nd order central schemes)
#                self.Domain.m_species[i][:,1:]    = aW[:,1:]    * m_c[i][:,:-1]
#                self.Domain.m_species[i][:,0]     = aE[:,0]     * m_c[i][:,1]
#                
#                self.Domain.m_species[i][:,1:-1] += aE[:,1:-1]  * m_c[i][:,2:]
#                self.Domain.m_species[i][1:,:]   += aS[1:,:]    * m_c[i][:-1,:]
#                self.Domain.m_species[i][:-1,:]  += aN[:-1,:]   * m_c[i][1:,:]
#                self.Domain.m_species[i]         -= (aW+aE+aS+aN)*m_c[i]
#            
#                # Species generated/destroyed during reaction
#                self.Domain.m_species[i]+=mole_ratio[i]*deta
#                
#                # Species advected from Porous medium equations [TO BE CONTINUED]
#                
#                
#                # Apply data from previous time step
#                self.Domain.m_species[i]*= dt
#                self.Domain.m_species[i]+= m_c[i]
##                print(self.Domain.m_species[i])
#                # IMPLICITLY MAKING SPECIES FLUX 0 AT BOUNDARIES
#                max_Y=max(np.amax(self.Domain.m_species[i]), max_Y)
#                min_Y=min(np.amin(self.Domain.m_species[i]), min_Y)
#        print(self.Domain.m_species)
        ###################################################################
        # Conservation of Energy
        ###################################################################
        # Heat diffusion
            #left faces
        self.Domain.E[:,1:]   -= dt/hx[:,1:]\
                    *self.interpolate(k[:,:-1],k[:,1:], 'Harmonic')\
                    *(T_c[:,1:]-T_c[:,:-1])/self.dx[:,:-1]
            # Right face
        self.Domain.E[:,:-1] += dt/hx[:,:-1]\
                    *self.interpolate(k[:,:-1],k[:,1:], 'Harmonic')\
                    *(T_c[:,1:]-T_c[:,:-1])/self.dx[:,:-1]
            # South face
        self.Domain.E[1:,:]   -= dt/hy[1:,:]\
                    *self.interpolate(k[1:,:],k[:-1,:], 'Harmonic')\
                    *(T_c[1:,:]-T_c[:-1,:])/self.dy[:-1,:]
            # North face
        self.Domain.E[:-1,:]  += dt/hy[:-1,:]\
                    *self.interpolate(k[:-1,:],k[1:,:], 'Harmonic')\
                    *(T_c[1:,:]-T_c[:-1,:])/self.dy[:-1,:]
        
        # Source terms
        self.Domain.E +=E_unif*dt
        self.Domain.E +=E_kim *dt
        
        # Porous medium advection
        if bool(self.Domain.rho_species):
                # Left face
            self.Domain.E[:,1:]+=dt/hx[:,1:]\
                *self.interpolate(rho_spec[species[0]][:,1:],rho_spec[species[0]][:,:-1],'Linear')*\
                (-perm/mu*(self.Domain.P[:,1:]-self.Domain.P[:,:-1])/self.dx[:,:-1])\
                *self.interpolate(Cp_spec[species[0]][:,1:],Cp_spec[species[0]][:,:-1],'Linear')\
                *self.interpolate(T_c[:,1:],T_c[:,:-1],'Linear')
                # Right face
            self.Domain.E[:,:-1]-=dt/hx[:,:-1]\
                *self.interpolate(rho_spec[species[0]][:,1:],rho_spec[species[0]][:,:-1],'Linear')*\
                (-perm/mu*(self.Domain.P[:,1:]-self.Domain.P[:,:-1])/self.dx[:,:-1])\
                *self.interpolate(Cp_spec[species[0]][:,1:],Cp_spec[species[0]][:,:-1],'Linear')\
                *self.interpolate(T_c[:,1:],T_c[:,:-1],'Linear')
                # South face
            self.Domain.E[1:,:]+=dt/hy[1:,:]\
                *self.interpolate(rho_spec[species[0]][1:,:],rho_spec[species[0]][:-1,:],'Linear')*\
                (-perm/mu*(self.Domain.P[1:,:]-self.Domain.P[:-1,:])/self.dy[:-1,:])\
                *self.interpolate(Cp_spec[species[0]][1:,:],Cp_spec[species[0]][:-1,:],'Linear')\
                *self.interpolate(T_c[1:,:],T_c[:-1,:],'Linear')
                # North face
            self.Domain.E[:-1,:]-=dt/hy[:-1,:]\
                *self.interpolate(rho_spec[species[0]][1:,:],rho_spec[species[0]][:-1,:],'Linear')*\
                (-perm/mu*(self.Domain.P[1:,:]-self.Domain.P[:-1,:])/self.dy[:-1,:])\
                *self.interpolate(Cp_spec[species[0]][1:,:],Cp_spec[species[0]][:-1,:],'Linear')\
                *self.interpolate(T_c[1:,:],T_c[:-1,:],'Linear')

        
#        # Radiation effects
#        self.Domain.T[1:-1,1:-1]+=0.8*5.67*10**(-8)*(T_c[:-2,1:-1]**4+T_c[2:,1:-1]**4+T_c[1:-1,:-2]**4+T_c[1:-1,2:]**4)
        
        # Apply boundary conditions
        self.BCs.Energy(self.Domain.E, T_c, dt, rho, Cv, hx, hy)
        
        ###################################################################
        # Divergence/Convergence checks
        ###################################################################
        if (np.isnan(np.amax(self.Domain.E))) \
        or (np.amin(self.Domain.E)<=0):
            return 2, dt
        elif (np.amax(self.Domain.eta)>1.0) or (np.amin(self.Domain.eta)<-10**(-9)):
            return 3, dt
        elif bool(self.Domain.rho_species) and ((min_Y<-10**(-9))\
                  or np.isnan(max_Y)):
            return 4, dt
        else:
            return 0, dt