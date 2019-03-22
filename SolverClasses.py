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
#import PorousClass

# 2D solver (Cartesian coordinates)
class TwoDimPlanarSolve():
    def __init__(self, geom_obj, settings, Sources, BCs, solver):
        self.Domain=geom_obj # Geometry object
        self.time_scheme=settings['Time_Scheme']
        self.dx,self.dy=np.meshgrid(geom_obj.dx,geom_obj.dy)
        self.BCs=BCs
        self.Fo=settings['Fo']
        self.dt=settings['dt']
        self.conv=settings['Convergence']
        self.countmax=settings['Max_iterations']
        
        # Define source terms and pointer to source object here
        self.get_source=Source_Comb.Source_terms(Sources['Ea'], Sources['A0'], Sources['dH'])
        self.source_unif=Sources['Source_Uniform']
        self.source_Kim=Sources['Source_Kim']
        
        # Porous medium solver
#        self.Porous_Eqns=PorousClass(0,0,0)
        
    # Time step check with dx, dy, Fo number
    def getdt(self, k, rho, Cv):
        # Stability check for Fourrier number
        if self.time_scheme=='Explicit':
            self.Fo=min(self.Fo, 1.0)
        elif self.Fo=='None':
            self.Fo=1.0
        
        dt=self.Fo*rho*Cv/k*self.Domain.CV_vol()
        return np.amin(dt)
    
    # Interpolation function
    def interpolate(self, k1, k2, func):
        if func=='Linear':
            return 0.5*k1+0.5*k2
        else:
            return 2*k1*k2/(k1+k2)
        
    # coefficients for temperature weighting in Advance_Soln_Cond
    def get_Coeff(self, dx, dy, dt, k, inter_type):
        aW=np.zeros_like(dx)
        aE=np.zeros_like(dx)
        aS=np.zeros_like(dx)
        aN=np.zeros_like(dx)
        
        # Left/right face factors
        aW[1:-1,1:-1] =0.5*self.interpolate(k[1:-1,1:-1], k[1:-1,:-2], inter_type)\
                    *(dy[1:-1,1:-1]+dy[:-2,1:-1])/(dx[1:-1,:-2])
        aE[1:-1,1:-1] =0.5*self.interpolate(k[1:-1,1:-1],k[1:-1,2:], inter_type)\
                    *(dy[1:-1,1:-1]+dy[:-2,1:-1])/(dx[1:-1,1:-1])
        # At north/south bondaries
        aW[0,1:-1]    =0.5*self.interpolate(k[0,1:-1],k[0,:-2], inter_type)\
            *(dy[0,1:-1])/(dx[0,:-2])
        aE[0,1:-1]    =0.5*self.interpolate(k[0,1:-1],k[0,2:], inter_type)\
            *(dy[0,1:-1])/(dx[0,1:-1])
        aW[-1,1:-1]   =0.5*self.interpolate(k[-1,1:-1],k[-1,:-2], inter_type)\
            *(dy[-1,1:-1])/(dx[-1,:-2])
        aE[-1,1:-1]   =0.5*self.interpolate(k[-1,1:-1],k[-1,2:], inter_type)\
            *(dy[-1,1:-1])/(dx[-1,1:-1])
        # At east/west boundaries
        aE[0,0]       =0.5*self.interpolate(k[0,0],k[0,1], inter_type)\
            *(dy[0,0])/dx[0,0]
        aE[1:-1,0]    =0.5*self.interpolate(k[1:-1,0],k[1:-1,1], inter_type)\
            *(dy[1:-1,0]+dy[:-2,0])/dx[1:-1,0]
        aE[-1,0]      =0.5*self.interpolate(k[-1,0],k[-1,1], inter_type)\
            *(dy[-1,0])/dx[-1,0]
        aW[0,-1]      =0.5*self.interpolate(k[0,-1],k[0,-2], inter_type)\
            *(dy[0,-1])/dx[0,-1]
        aW[1:-1,-1]   =0.5*self.interpolate(k[1:-1,-1],k[1:-1,-2], inter_type)\
            *(dy[1:-1,-1]+dy[:-2,-1])/dx[1:-1,-1]
        aW[-1,-1]     =0.5*self.interpolate(k[-1,-1],k[-1,-2], inter_type)\
            *(dy[-1,-1])/dx[-1,-1]
        
        # Top/bottom faces
        aS[1:-1,1:-1]=0.5*self.interpolate(k[1:-1,1:-1],k[:-2,1:-1], inter_type)\
            *(dx[1:-1,1:-1]+dx[1:-1,:-2])/dy[:-2,1:-1]
        aN[1:-1,1:-1]=0.5*self.interpolate(k[1:-1,1:-1],k[2:,1:-1], inter_type)\
            *(dx[1:-1,1:-1]+dx[1:-1,:-2])/dy[1:-1,1:-1]
        
        # Area account for east/west boundary nodes
        aS[1:-1,0]    =0.5*self.interpolate(k[1:-1,0],k[:-2,0], inter_type)\
            *(dx[1:-1,0])/(dy[:-2,0])
        aN[1:-1,0]    =0.5*self.interpolate(k[1:-1,0],k[2:,0], inter_type)\
            *(dx[1:-1,0])/(dy[1:-1,0])
        aS[1:-1,-1]   =0.5*self.interpolate(k[1:-1,-1],k[:-2,-1], inter_type)\
            *(dx[1:-1,-1])/(dy[:-2,-1])
        aN[1:-1,-1]   =0.5*self.interpolate(k[1:-1,-1],k[2:,-1], inter_type)\
            *(dx[1:-1,-1])/(dy[1:-1,-1])
        # Forward/backward difference for north/south boundaries
        aN[0,0]       =0.5*self.interpolate(k[0,0],k[1,0], inter_type)\
            *dx[0,0]/dy[0,0]
        aN[0,1:-1]    =0.5*self.interpolate(k[0,1:-1],k[1,1:-1], inter_type)\
            *(dx[0,1:-1]+dx[0,:-2])/dy[0,1:-1]
        aN[0,-1]      =0.5*self.interpolate(k[0,-1],k[1,-1], inter_type)\
            *dx[0,-1]/dy[0,-1]
        aS[-1,0]      =0.5*self.interpolate(k[-1,0],k[-2,0], inter_type)\
            *dx[-1,0]/dy[-1,0]
        aS[-1,1:-1]   =0.5*self.interpolate(k[-1,1:-1],k[-2,1:-1], inter_type)\
            *(dx[0,1:-1]+dx[0,:-2])/dy[-1,1:-1]
        aS[-1,-1]     =0.5*self.interpolate(k[-1,-1],k[-2,-1], inter_type)\
            *dx[-1,-1]/dy[-1,-1]
        
        return aW,aE,aS,aN
    
    # Bondary condition handler
    def Apply_BCs_Cond(self, E, T_prev, dt, rho, Cv, vol):
        # Left face
        for i in range(len(self.BCs['bc_left'])/3):
            st=self.BCs['bc_left'][2+3*i][0]
            en=self.BCs['bc_left'][2+3*i][1]
            if self.BCs['bc_left'][3*i]=='T':
                E[st:en,0]=self.BCs['bc_left'][1+3*i]*rho[st:en,0]*Cv[st:en,0]*vol[st:en,0]
                if len(self.BCs['bc_left'])/3-i==1:
                    E[-1,0]=self.BCs['bc_left'][-2]*rho[-1,0]*Cv[-1,0]*vol[-1,0]
            
            else:
                if self.BCs['bc_left'][3*i]=='F':
                    q=self.BCs['bc_left'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_left'][1+3*i][0]*self.BCs['bc_left'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_left'][1+3*i][0]*T_prev[st:en,0] # h*Tij
                
                E[st:en,0]+=(Bi+q)*self.dy[st:en,0]*dt
                if len(self.BCs['bc_left'])/3-i==1:
                    if self.BCs['bc_left'][3*i]=='C':
                        Bi=-self.BCs['bc_left'][1+3*i][0]*T_prev[-1,0] # h*Tij
                    E[-1,0]+=(Bi+q)*self.dy[-1,0]*dt
        
        # Right face
        for i in range(len(self.BCs['bc_right'])/3):
            st=self.BCs['bc_right'][2+3*i][0]
            en=self.BCs['bc_right'][2+3*i][1]
            if self.BCs['bc_right'][3*i]=='T':
                E[st:en,-1]=self.BCs['bc_right'][1+3*i]*rho[st:en,-1]*Cv[st:en,-1]*vol[st:en,-1]
                if len(self.BCs['bc_right'])/3-i==1:
                    E[-1,-1]=self.BCs['bc_right'][-2]*rho[-1,-1]*Cv[-1,-1]*vol[-1,-1]
            
            else:
                if self.BCs['bc_right'][3*i]=='F':
                    q=self.BCs['bc_right'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_right'][1+3*i][0]*self.BCs['bc_right'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_right'][1+3*i][0]*T_prev[st:en,-1] # h*Tij
                
                E[st:en,-1]+=(Bi+q)*self.dy[st:en,-1]*dt
                if len(self.BCs['bc_right'])/3-i==1:
                    if self.BCs['bc_right'][3*i]=='C':
                        Bi=-self.BCs['bc_right'][1+3*i][0]*T_prev[-1,-1] # h*Tij
                    E[-1,-1]+=(Bi+q)*self.dy[-1,-1]*dt
        
        # South face
        for i in range(len(self.BCs['bc_south'])/3):
            st=self.BCs['bc_south'][2+3*i][0]
            en=self.BCs['bc_south'][2+3*i][1]
            if self.BCs['bc_south'][3*i]=='T':
                E[0,st:en]=self.BCs['bc_south'][1+3*i]*rho[0,st:en]*Cv[0,st:en]*vol[0,st:en]
                if len(self.BCs['bc_south'])/3-i==1:
                    E[0,-1]=self.BCs['bc_south'][-2]*rho[0,-1]*Cv[0,-1]*vol[0,-1]
            
            else:
                if self.BCs['bc_south'][3*i]=='F':
                    q=self.BCs['bc_south'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_south'][1+3*i][0]*self.BCs['bc_south'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_south'][1+3*i][0]*T_prev[0,st:en] # h*Tij
                
                E[0,st:en]+=(Bi+q)*self.dx[0,st:en]*dt
                if len(self.BCs['bc_south'])/3-i==1:
                    if self.BCs['bc_south'][3*i]=='C':
                        Bi=-self.BCs['bc_south'][1+3*i][0]*T_prev[0,-1] # h*Tij
                    E[0,-1]+=(Bi+q)*self.dx[0,-1]*dt
                    
        # North face
        for i in range(len(self.BCs['bc_north'])/3):
            st=self.BCs['bc_north'][2+3*i][0]
            en=self.BCs['bc_north'][2+3*i][1]
            if self.BCs['bc_north'][3*i]=='T':
                E[-1,st:en]=self.BCs['bc_north'][1+3*i]*rho[-1,st:en]*Cv[-1,st:en]*vol[-1,st:en]
                if len(self.BCs['bc_north'])/3-i==1:
                    E[-1,-1]=self.BCs['bc_north'][-2]*rho[-1,-1]*Cv[-1,-1]*vol[-1,-1]
            
            else:
                if self.BCs['bc_north'][3*i]=='F':
                    q=self.BCs['bc_north'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_north'][1+3*i][0]*self.BCs['bc_north'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_north'][1+3*i][0]*T_prev[-1,st:en] # h*Tij
                
                E[-1,st:en]+=(Bi+q)*self.dx[-1,st:en]*dt
                if len(self.BCs['bc_north'])/3-i==1:
                    if self.BCs['bc_north'][3*i]=='C':
                        Bi=-self.BCs['bc_north'][1+3*i][0]*T_prev[-1,-1] # h*Tij
                    E[-1,-1]+=(Bi+q)*self.dx[-1,-1]*dt
        
        # Apply radiation BCs
        if self.BCs['bc_left_rad']!='None':
            E[:,0]+=self.dy[:,0]*dt*\
                self.BCs['bc_left_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_left_rad'][1]**4-T_prev[:,0]**4)
        if self.BCs['bc_right_rad']!='None':
            E[:,-1]+=self.dy[:,-1]*dt*\
                self.BCs['bc_right_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_right_rad'][1]**4-T_prev[:,-1]**4)
        if self.BCs['bc_south_rad']!='None':
            E[0,:]+=self.dx[0,:]*dt*\
                self.BCs['bc_south_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_south_rad'][1]**4-T_prev[0,:]**4)
        if self.BCs['bc_north_rad']!='None':
            E[-1,:]+=self.dx[-1,:]*dt*\
                self.BCs['bc_north_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_north_rad'][1]**4-T_prev[-1,:]**4)
    
    # Main solver (1 time step)
    def Advance_Soln_Cond(self, nt, t, vol):
        max_Y,min_Y=0,1
        # Calculate properties
        k, rho, Cv, D=self.Domain.calcProp()
        
        if self.dt=='None':
            dt=self.getdt(k, rho, Cv)
        else:
            dt=min(self.dt,self.getdt(k, rho, Cv))
            
        if (np.isnan(dt)) or (dt<=0):
            print '*********Diverging time step***********'
            return 1, dt
        print 'Time step %i, Step size=%.7f, Time elapsed=%f;'%(nt+1,dt, t+dt)
        
        # Copy needed variables for conservation equations
        T_c=self.Domain.TempFromConserv()
        Y_c=copy.deepcopy(self.Domain.Y_species)
        E_0=copy.deepcopy(self.Domain.E)
        
        ###################################################################
        # Calculate source and Porous medium terms
        ###################################################################
        # Source terms
        E_unif,E_kim=0,0
        if self.source_unif!='None':
            E_unif      = self.get_source.Source_Uniform(self.source_unif, vol)
        if self.source_Kim=='True':
#            self.Domain.eta=self.Domain.Y_species[:,:,2]/0.25
            E_kim, deta =self.get_source.Source_Comb_Kim(rho, T_c, self.Domain.eta, vol, dt)
#            E_kim, deta =self.get_source.Source_Comb_Umbrajkar(rho, T_c, self.Domain.eta, self.Domain.CV_vol(), dt)
            
        # Porous medium equations [TO BE CONTINUED]
#        self.Porous_Eqns
        
        ###################################################################
        # Conservation of species
        ###################################################################
        if bool(self.Domain.Y_species):
            # Mole ratios
            mole_ratio={}
            mole_ratio[self.Domain.species_keys[0]]=-2.0/5 # Al
            mole_ratio[self.Domain.species_keys[1]]=-3.0/5 # CuO
            mole_ratio[self.Domain.species_keys[2]]=1.0/4  # Al2O3
            mole_ratio[self.Domain.species_keys[3]]=3.0/4  # Cu
            
            for i in self.Domain.species_keys:
                # Calculate flux coefficients
                aW,aE,aS,aN=self.get_Coeff(self.dx,self.dy, dt, rho*D[i], 'Linear')
                
                # Diffusion contribution (2nd order central schemes)
                self.Domain.Y_species[i][:,1:]    = aW[:,1:]    * Y_c[i][:,:-1]
                self.Domain.Y_species[i][:,0]     = aE[:,0]     * Y_c[i][:,1]
                
                self.Domain.Y_species[i][:,1:-1] += aE[:,1:-1]  * Y_c[i][:,2:]
                self.Domain.Y_species[i][1:,:]   += aS[1:,:]    * Y_c[i][:-1,:]
                self.Domain.Y_species[i][:-1,:]  += aN[:-1,:]   * Y_c[i][1:,:]
                self.Domain.Y_species[i]         -= (aW+aE+aS+aN)*Y_c[i]
            
                # Species generated/destroyed during reaction
                self.Domain.Y_species[i]+=mole_ratio[i]*deta
                
                # Species advected from Porous medium equations [TO BE CONTINUED]
                
                
                # Apply data from previous time step
                self.Domain.Y_species[i]*= dt
                self.Domain.Y_species[i]+= Y_c[i]
#                print(self.Domain.Y_species[i])
                # IMPLICITLY MAKING SPECIES FLUX 0 AT BOUNDARIES
                max_Y=max(np.amax(self.Domain.Y_species[i]), max_Y)
                min_Y=max(np.amin(self.Domain.Y_species[i]), min_Y)
#        print(self.Domain.Y_species)
        ###################################################################
        # Conservation of Energy
        ###################################################################
        # Calculate flux coefficients
        aW,aE,aS,aN=self.get_Coeff(self.dx,self.dy, dt, k, 'Harmonic')
        
        # Heat diffusion contribution (2nd order central schemes)
        self.Domain.E[:,1:]    = aW[:,1:]*T_c[:,:-1]
        self.Domain.E[:,0]     = aE[:,0]*T_c[:,1]
        
        self.Domain.E[:,1:-1] += aE[:,1:-1]*T_c[:,2:]
        self.Domain.E[1:,:]   += aS[1:,:]*T_c[:-1,:]
        self.Domain.E[:-1,:]  += aN[:-1,:]*T_c[1:,:]
        self.Domain.E         -= (aW+aE+aS+aN)*T_c
        
        # Source terms
        self.Domain.E +=E_unif
        self.Domain.E +=E_kim
        
        # Porous medium advection [TO BE CONTINUED]
        
        
#        # Radiation effects
#        self.Domain.T[1:-1,1:-1]+=0.8*5.67*10**(-8)*(T_c[:-2,1:-1]**4+T_c[2:,1:-1]**4+T_c[1:-1,:-2]**4+T_c[1:-1,2:]**4)
        
        # Apply energy from previous time step and boundary conditions
        self.Domain.E*= dt
        self.Domain.E+= E_0
        self.Apply_BCs_Cond(self.Domain.E, T_c, dt, rho, Cv, vol)
        
        ###################################################################
        # Divergence/Convergence checks
        ###################################################################
        if (np.isnan(np.amax(self.Domain.E))) \
        or (np.amax(self.Domain.E)>100*np.amax(E_0))\
        or (np.amin(self.Domain.E)<=0):
            print '***********Divergence detected - energy************'
            return 1, dt
        elif (np.amax(self.Domain.eta)>1.0) or (np.amin(self.Domain.eta)<-10**(-9)):
            print '***********Divergence detected - reaction progress************'
            return 1, dt
        elif bool(self.Domain.Y_species) and ((max_Y>1.0) or (min_Y<-10**(-9))):
            print '***********Divergence detected - species mass fraction************'
            return 1, dt
        else:
            return 0, dt
        
# 2D axisymmetric solver
class AxisymmetricSolve():
    def __init__(self, geom_obj, settings, Sources, BCs, solver):
        self.Domain=geom_obj # Geometry object
        self.time_scheme=settings['Time_Scheme']
        self.dx,self.dy=np.meshgrid(geom_obj.dx,geom_obj.dy)
        self.BCs=BCs
        self.Fo=settings['Fo']
        self.dt=settings['dt']
        self.conv=settings['Convergence']
        self.countmax=settings['Max_iterations']
        
        # Define source terms and pointer to source object here
        self.get_source=Source_Comb.Source_terms(Sources['Ea'], Sources['A0'], Sources['dH'])
        self.source_unif=Sources['Source_Uniform']
        self.source_Kim=Sources['Source_Kim']
        
        # Porous medium solver
#        self.Porous_Eqns=PorousClass(0,0,0)
        
    # Time step check with dx, dy, Fo number
    def getdt(self, k, rho, Cv, vol):
        # Stability check for Fourrier number
#        if self.time_scheme=='Explicit':
#            self.Fo=min(self.Fo, 0.85)
#        elif self.Fo=='None':
#            self.Fo=1.0
        
        dt=self.Fo*rho*Cv/k*vol
        return np.amin(dt)
    
    # Interpolation function
    def interpolate(self, k1, k2, func):
        if func=='Linear':
            return 0.5*k1+0.5*k2
        else:
            return 2*k1*k2/(k1+k2)
        
    # coefficients for temperature weighting in Advance_Soln_Cond
    def get_Coeff(self, dx, dy, dt, k, inter_type):
        aW=np.zeros_like(dx)
        aE=np.zeros_like(dx)
        aS=np.zeros_like(dx)
        aN=np.zeros_like(dx)
        
        # Left/right face factors
        aW[1:-1,1:-1] =0.5*self.interpolate(k[1:-1,1:-1], k[1:-1,:-2], inter_type)\
                    *(self.Domain.X[1:-1,1:-1]-dx[1:-1,:-2]/2)\
                    *(dy[1:-1,1:-1]+dy[:-2,1:-1])/(dx[1:-1,:-2])
        aE[1:-1,1:-1] =0.5*self.interpolate(k[1:-1,1:-1],k[1:-1,2:], inter_type)\
                    *(self.Domain.X[1:-1,1:-1]+dx[1:-1,1:-1]/2)*(dy[1:-1,1:-1]+dy[:-2,1:-1])/(dx[1:-1,1:-1])
        # At north/south bondaries
        aW[0,1:-1]    =0.5*self.interpolate(k[0,1:-1],k[0,:-2], inter_type)\
            *(self.Domain.X[0,1:-1]-dx[0,:-2]/2)*(dy[0,1:-1])/(dx[0,:-2])
        aE[0,1:-1]    =0.5*self.interpolate(k[0,1:-1],k[0,2:], inter_type)\
            *(self.Domain.X[0,1:-1]+dx[0,1:-1]/2)*(dy[0,1:-1])/(dx[0,1:-1])
        aW[-1,1:-1]   =0.5*self.interpolate(k[-1,1:-1],k[-1,:-2], inter_type)\
            *(self.Domain.X[-1,1:-1]-dx[-1,:-2]/2)*(dy[-1,1:-1])/(dx[-1,:-2])
        aE[-1,1:-1]   =0.5*self.interpolate(k[-1,1:-1],k[-1,2:], inter_type)\
            *(self.Domain.X[-1,1:-1]+dx[-1,1:-1]/2)*(dy[-1,1:-1])/(dx[-1,1:-1])
        # At west/east boundaries
        aE[0,0]       =0.5*self.interpolate(k[0,0],k[0,1], inter_type)\
            *(self.Domain.X[0,0]+dx[0,0]/2)*(dy[0,0])/dx[0,0]
        aE[1:-1,0]    =0.5*self.interpolate(k[1:-1,0],k[1:-1,1], inter_type)\
            *(self.Domain.X[1:-1,0]+dx[1:-1,0]/2)*(dy[1:-1,0]+dy[:-2,0])/dx[1:-1,0]
        aE[-1,0]      =0.5*self.interpolate(k[-1,0],k[-1,1], inter_type)\
            *(self.Domain.X[-1,0]+dx[-1,0]/2)*(dy[-1,0])/dx[-1,0]
        aW[0,-1]      =0.5*self.interpolate(k[0,-1],k[0,-2], inter_type)\
            *(self.Domain.X[0,-1]-dx[0,-1]/2)*(dy[0,-1])/dx[0,-1]
        aW[1:-1,-1]   =0.5*self.interpolate(k[1:-1,-1],k[1:-1,-2], inter_type)\
            *(self.Domain.X[1:-1,-1]-dx[1:-1,-1]/2)*(dy[1:-1,-1]+dy[:-2,-1])/dx[1:-1,-1]
        aW[-1,-1]     =0.5*self.interpolate(k[-1,-1],k[-1,-2], inter_type)\
            *(self.Domain.X[-1,-1]-dx[-1,-1]/2)*(dy[-1,-1])/dx[-1,-1]
        
        # Top/bottom faces
        aS[1:-1,1:-1]=0.5*self.interpolate(k[1:-1,1:-1],k[:-2,1:-1], inter_type)\
            *(self.Domain.X[1:-1,1:-1]-dx[1:-1,:-2]/2)*(dx[1:-1,1:-1]+dx[1:-1,:-2])/dy[:-2,1:-1]
        aN[1:-1,1:-1]=0.5*self.interpolate(k[1:-1,1:-1],k[2:,1:-1], inter_type)\
            *(self.Domain.X[1:-1,1:-1]-dx[1:-1,:-2]/2)*(dx[1:-1,1:-1]+dx[1:-1,:-2])/dy[1:-1,1:-1]
        
        # Area account for east/west boundary nodes
        aS[1:-1,0]    =0.5*self.interpolate(k[1:-1,0],k[:-2,0], inter_type)\
            *(self.Domain.X[1:-1,0])*(dx[1:-1,0])/(dy[:-2,0])
        aN[1:-1,0]    =0.5*self.interpolate(k[1:-1,0],k[2:,0], inter_type)\
            *(self.Domain.X[1:-1,0])*(dx[1:-1,0])/(dy[1:-1,0])
        aS[1:-1,-1]   =0.5*self.interpolate(k[1:-1,-1],k[:-2,-1], inter_type)\
            *(self.Domain.X[1:-1,-1]-dx[1:-1,-1]/2)*(dx[1:-1,-1])/(dy[:-2,-1])
        aN[1:-1,-1]   =0.5*self.interpolate(k[1:-1,-1],k[2:,-1], inter_type)\
            *(self.Domain.X[1:-1,-1]-dx[1:-1,-1]/2)*(dx[1:-1,-1])/(dy[1:-1,-1])
        # Forward/backward difference for north/south boundaries
        aN[0,0]       =0.5*self.interpolate(k[0,0],k[1,0], inter_type)\
            *(self.Domain.X[0,0])*dx[0,0]/dy[0,0]
        aN[0,1:-1]    =0.5*self.interpolate(k[0,1:-1],k[1,1:-1], inter_type)\
            *(self.Domain.X[0,1:-1]-dx[0,:-2]/2)*(dx[0,1:-1]+dx[0,:-2])/dy[0,1:-1]
        aN[0,-1]      =0.5*self.interpolate(k[0,-1],k[1,-1], inter_type)\
            *(self.Domain.X[0,-1]-dx[0,-1]/2)*dx[0,-1]/dy[0,-1]
        aS[-1,0]      =0.5*self.interpolate(k[-1,0],k[-2,0], inter_type)\
            *(self.Domain.X[-1,0])*dx[-1,0]/dy[-1,0]
        aS[-1,1:-1]   =0.5*self.interpolate(k[-1,1:-1],k[-2,1:-1], inter_type)\
            *(self.Domain.X[-1,1:-1]-dx[-1,:-2]/2)*(dx[0,1:-1]+dx[0,:-2])/dy[-1,1:-1]
        aS[-1,-1]     =0.5*self.interpolate(k[-1,-1],k[-2,-1], inter_type)\
            *(self.Domain.X[-1,-1]-dx[-1,-1]/2)*dx[-1,-1]/dy[-1,-1]
        
        return aW,aE,aS,aN
    
    # Bondary condition handler
    def Apply_BCs_Cond(self, E, T_prev, dt, rho, Cv, vol):
        # Left face (zero flux implied unless prescribed temperature)
        if self.BCs['bc_left'][0]=='T':
            E[:,0]=self.BCs['bc_left'][1]*rho[:,0]*Cv[:,0]*vol[:,0]
            
        # Right face
        for i in range(len(self.BCs['bc_right'])/3):
            st=self.BCs['bc_right'][2+3*i][0]
            en=self.BCs['bc_right'][2+3*i][1]
            if self.BCs['bc_right'][3*i]=='T':
                E[st:en,-1]=self.BCs['bc_right'][1+3*i]*rho[st:en,-1]*Cv[st:en,-1]*vol[st:en,-1]
                if len(self.BCs['bc_right'])/3-i==1:
                    E[-1,-1]=self.BCs['bc_right'][-2]*rho[-1,-1]*Cv[-1,-1]*vol[-1,-1]
            
            else:
                if self.BCs['bc_right'][3*i]=='F':
                    q=self.BCs['bc_right'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_right'][1+3*i][0]*self.BCs['bc_right'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_right'][1+3*i][0]*T_prev[st:en,-1] # h*Tij
                
                E[st:en,-1]+=(Bi+q)*(self.Domain.X[st:en,-1]*self.dy[st:en,-1])*dt
                if len(self.BCs['bc_right'])/3-i==1:
                    if self.BCs['bc_right'][3*i]=='C':
                        Bi=-self.BCs['bc_right'][1+3*i][0]*T_prev[-1,-1] # h*Tij
                    E[-1,-1]+=(Bi+q)*(self.Domain.X[-1,-1]*self.dy[-1,-1])*dt
        
        # South face
        for i in range(len(self.BCs['bc_south'])/3):
            st=self.BCs['bc_south'][2+3*i][0]
            en=self.BCs['bc_south'][2+3*i][1]
            if self.BCs['bc_south'][3*i]=='T':
                E[0,st:en]=self.BCs['bc_south'][1+3*i]*rho[0,st:en]*Cv[0,st:en]*vol[0,st:en]
                if len(self.BCs['bc_south'])/3-i==1:
                    E[0,-1]=self.BCs['bc_south'][-2]*rho[0,-1]*Cv[0,-1]*vol[0,-1]
            
            else:
                if self.BCs['bc_south'][3*i]=='F':
                    q=self.BCs['bc_south'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_south'][1+3*i][0]*self.BCs['bc_south'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_south'][1+3*i][0]*T_prev[0,st:en] # h*Tij
                
                E[0,st:en]+=(Bi+q)*(self.Domain.X[0,st:en]*self.dx[0,st:en])*dt
                if len(self.BCs['bc_south'])/3-i==1:
                    if self.BCs['bc_south'][3*i]=='C':
                        Bi=-self.BCs['bc_south'][1+3*i][0]*T_prev[0,-1] # h*Tij
                    E[0,-1]+=(Bi+q)*(self.Domain.X[0,-1]*self.dx[0,-1])*dt
                    
        # North face
        for i in range(len(self.BCs['bc_north'])/3):
            st=self.BCs['bc_north'][2+3*i][0]
            en=self.BCs['bc_north'][2+3*i][1]
            if self.BCs['bc_north'][3*i]=='T':
                E[-1,st:en]=self.BCs['bc_north'][1+3*i]*rho[-1,st:en]*Cv[-1,st:en]*vol[-1,st:en]
                if len(self.BCs['bc_north'])/3-i==1:
                    E[-1,-1]=self.BCs['bc_north'][-2]*rho[-1,-1]*Cv[-1,-1]*vol[-1,-1]
            
            else:
                if self.BCs['bc_north'][3*i]=='F':
                    q=self.BCs['bc_north'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_north'][1+3*i][0]*self.BCs['bc_north'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_north'][1+3*i][0]*T_prev[-1,st:en] # h*Tij
                
                E[-1,st:en]+=(Bi+q)*(self.Domain.X[-1,st:en]*self.dx[-1,st:en])*dt
                if len(self.BCs['bc_north'])/3-i==1:
                    if self.BCs['bc_north'][3*i]=='C':
                        Bi=-self.BCs['bc_north'][1+3*i][0]*T_prev[-1,-1] # h*Tij
                    E[-1,-1]+=(Bi+q)*(self.Domain.X[-1,-1]*self.dx[-1,-1])*dt
        
        # Apply radiation BCs
        if self.BCs['bc_right_rad']!='None':
            E[:,-1]+=self.Domain.X[:,-1]*self.dy[:,-1]*dt*\
                self.BCs['bc_right_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_right_rad'][1]**4-T_prev[:,-1]**4)
        if self.BCs['bc_south_rad']!='None':
            E[0,:]+=self.Domain.X[0,:]*self.dx[0,:]*dt*\
                self.BCs['bc_south_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_south_rad'][1]**4-T_prev[0,:]**4)
        if self.BCs['bc_north_rad']!='None':
            E[-1,:]+=self.Domain.X[-1,:]*self.dx[-1,:]*dt*\
                self.BCs['bc_north_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_north_rad'][1]**4-T_prev[-1,:]**4)
    
    # Main solver (1 time step)
    def Advance_Soln_Cond(self, nt, t, vol):
        max_Y,min_Y=0,1
        # Calculate properties
        k, rho, Cv, D=self.Domain.calcProp()
        
        if self.dt=='None':
            dt=self.getdt(k, rho, Cv, vol)
        else:
            dt=min(self.dt,self.getdt(k, rho, Cv, vol))
#            dt=self.dt
        if (np.isnan(dt)) or (dt<=0):
            print '*********Diverging time step***********'
            return 1, dt
        print 'Time step %i, Step size=%.7f, Time elapsed=%f;'%(nt+1,dt, t+dt)
        
        # Copy needed variables for conservation equations
        T_c=self.Domain.TempFromConserv()
        Y_c=copy.deepcopy(self.Domain.Y_species)
        E_0=copy.deepcopy(self.Domain.E)
        
        ###################################################################
        # Calculate source and Porous medium terms
        ###################################################################
        # Source terms
        E_unif,E_kim=0,0
        if self.source_unif!='None':
            E_unif      = self.get_source.Source_Uniform(self.source_unif, vol)
        if self.source_Kim=='True':
#            self.Domain.eta=self.Domain.Y_species[:,:,2]/0.25
            E_kim, deta =self.get_source.Source_Comb_Kim(rho, T_c, self.Domain.eta, vol, dt)
#            E_kim, deta =self.get_source.Source_Comb_Umbrajkar(rho, T_c, self.Domain.eta, self.Domain.CV_vol(), dt)
            
        # Porous medium equations [TO BE CONTINUED]
#        self.Porous_Eqns
        
        ###################################################################
        # Conservation of species
        ###################################################################
        if bool(self.Domain.Y_species):
            # Mole ratios
            mole_ratio={}
            mole_ratio[self.Domain.species_keys[0]]=-2.0/5 # Al
            mole_ratio[self.Domain.species_keys[1]]=-3.0/5 # CuO
            mole_ratio[self.Domain.species_keys[2]]=1.0/4  # Al2O3
            mole_ratio[self.Domain.species_keys[3]]=3.0/4  # Cu
            
            for i in self.Domain.species_keys:
                # Calculate flux coefficients
                aW,aE,aS,aN=self.get_Coeff(self.dx,self.dy, dt, rho*D[i], 'Linear')
                
                # Diffusion contribution (2nd order central schemes)
                self.Domain.Y_species[i][:,1:]    = aW[:,1:]    * Y_c[i][:,:-1]
                self.Domain.Y_species[i][:,0]     = aE[:,0]     * Y_c[i][:,1]
                
                self.Domain.Y_species[i][:,1:-1] += aE[:,1:-1]  * Y_c[i][:,2:]
                self.Domain.Y_species[i][1:,:]   += aS[1:,:]    * Y_c[i][:-1,:]
                self.Domain.Y_species[i][:-1,:]  += aN[:-1,:]   * Y_c[i][1:,:]
                self.Domain.Y_species[i]         -= (aW+aE+aS+aN)*Y_c[i]
            
                # Species generated/destroyed during reaction
                self.Domain.Y_species[i]+=mole_ratio[i]*deta
                
                # Species advected from Porous medium equations [TO BE CONTINUED]
                
                
                # Apply data from previous time step
                self.Domain.Y_species[i]*= dt
                self.Domain.Y_species[i]+= Y_c[i]
#                print(self.Domain.Y_species[i])
                # IMPLICITLY MAKING SPECIES FLUX 0 AT BOUNDARIES
                max_Y=max(np.amax(self.Domain.Y_species[i]), max_Y)
                min_Y=max(np.amin(self.Domain.Y_species[i]), min_Y)
#        print(self.Domain.Y_species)
        ###################################################################
        # Conservation of Energy
        ###################################################################
        # Calculate flux coefficients
        aW,aE,aS,aN=self.get_Coeff(self.dx,self.dy, dt, k, 'Harmonic')
        
        # Heat diffusion contribution (2nd order central schemes)
        self.Domain.E[:,1:]    = aW[:,1:]*T_c[:,:-1]
        self.Domain.E[:,0]     = aE[:,0]*T_c[:,1]
        
        self.Domain.E[:,1:-1] += aE[:,1:-1]*T_c[:,2:]
        self.Domain.E[1:,:]   += aS[1:,:]*T_c[:-1,:]
        self.Domain.E[:-1,:]  += aN[:-1,:]*T_c[1:,:]
        self.Domain.E         -= (aW+aE+aS+aN)*T_c
        
        # Source terms
        self.Domain.E +=E_unif
        self.Domain.E +=E_kim
        
        # Porous medium advection [TO BE CONTINUED]
        
        
#        # Radiation effects
#        self.Domain.T[1:-1,1:-1]+=0.8*5.67*10**(-8)*(T_c[:-2,1:-1]**4+T_c[2:,1:-1]**4+T_c[1:-1,:-2]**4+T_c[1:-1,2:]**4)
        
        # Apply energy from previous time step and boundary conditions
        self.Domain.E*= dt
        self.Domain.E+= E_0
        self.Apply_BCs_Cond(self.Domain.E, T_c, dt, rho, Cv, vol)
        
        ###################################################################
        # Divergence/Convergence checks
        ###################################################################
        if (np.isnan(np.amax(self.Domain.E))) \
        or (np.amax(self.Domain.E)>100*np.amax(E_0))\
        or (np.amin(self.Domain.E)<=0):
            print '***********Divergence detected - energy************'
            return 1, dt
        elif (np.amax(self.Domain.eta)>1.0) or (np.amin(self.Domain.eta)<-10**(-9)):
            print '***********Divergence detected - reaction progress************'
            return 1, dt
        elif bool(self.Domain.Y_species) and ((max_Y>1.0) or (min_Y<-10**(-9))):
            print '***********Divergence detected - species mass fraction************'
            return 1, dt
        else:
            return 0, dt