# -*- coding: utf-8 -*-
"""
######################################################
#             2D Heat Conduction Solver              #
#              Created by J. Mark Epps               #
#          Part of Masters Thesis at UW 2018-2020    #
######################################################

This file contains the solver classes for 2D planar Heat Conduction:
    -Modifies geometry class variables directly
    -Calculate time step based on Fourier number
    -Compute 2D conduction equations
    -Add source terms as needed
    -Applies boundary conditions (can vary along a side)

Features/assumptions:
    -time step based on Fourrier number and local discretizations in x and y
    -equal node spacing in x or y
    -thermal properties can vary in space (call from geometry object)
    -Radiation boundary conditions

"""

import numpy as np
#import CoolProp.CoolProp as CP
#import temporal_schemes
import Source_Comb
#import PorousClass

# 2D solver
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
    
    # coefficients for temperature weighting in Advance_Soln_Cond
    def get_Coeff(self, dx, dy, dt, k, rho, Cv):
        aW=np.zeros_like(dx)
        aE=np.zeros_like(dx)
        aS=np.zeros_like(dx)
        aN=np.zeros_like(dx)
        
        # Left/right face factors
        aW[1:-1,1:-1] =0.5*(2*k[1:-1,1:-1]*k[1:-1,:-2])/(k[1:-1,1:-1]+k[1:-1,:-2])\
                    *(dy[1:-1,1:-1]+dy[:-2,1:-1])/(dx[1:-1,:-2])
        aE[1:-1,1:-1] =0.5*(2*k[1:-1,1:-1]*k[1:-1,2:])/(k[1:-1,1:-1]+k[1:-1,2:])\
                    *(dy[1:-1,1:-1]+dy[:-2,1:-1])/(dx[1:-1,1:-1])
        # At north/south bondaries
        aW[0,1:-1]    =0.5*(2*k[0,1:-1]*k[0,:-2])/(k[0,1:-1]+k[0,:-2])\
            *(dy[0,1:-1])/(dx[0,:-2])
        aE[0,1:-1]    =0.5*(2*k[0,1:-1]*k[0,2:])/(k[0,1:-1]+k[0,2:])\
            *(dy[0,1:-1])/(dx[0,1:-1])
        aW[-1,1:-1]   =0.5*(2*k[-1,1:-1]*k[-1,:-2])/(k[-1,1:-1]+k[-1,:-2])\
            *(dy[-1,1:-1])/(dx[-1,:-2])
        aE[-1,1:-1]   =0.5*(2*k[-1,1:-1]*k[-1,2:])/(k[-1,1:-1]+k[-1,2:])\
            *(dy[-1,1:-1])/(dx[-1,1:-1])
        # At Left/right boundaries
        aE[0,0]       =0.5*(2*k[0,0]*k[0,1])/(k[0,0]+k[0,1])\
            *(dy[0,0])/dx[0,0]
        aE[1:-1,0]    =0.5*(2*k[1:-1,0]*k[1:-1,1])/(k[1:-1,0]+k[1:-1,1])\
            *(dy[1:-1,0]+dy[:-2,0])/dx[1:-1,0]
        aE[-1,0]      =0.5*(2*k[-1,0]*k[-1,1])/(k[-1,0]+k[-1,1])\
            *(dy[-1,0])/dx[-1,0]
        aW[0,-1]      =0.5*(2*k[0,-1]*k[0,-2])/(k[0,-1]+k[0,-2])\
            *(dy[0,-1])/dx[0,-1]
        aW[1:-1,-1]   =0.5*(2*k[1:-1,-1]*k[1:-1,-2])/(k[1:-1,-1]+k[1:-1,-2])\
            *(dy[1:-1,-1]+dy[:-2,-1])/dx[1:-1,-1]
        aW[-1,-1]     =0.5*(2*k[-1,-1]*k[-1,-2])/(k[-1,-1]+k[-1,-2])\
            *(dy[-1,-1])/dx[-1,-1]
        
        # South/north faces
        aS[1:-1,1:-1]=0.5*(2*k[1:-1,1:-1]*k[:-2,1:-1])/(k[1:-1,1:-1]+k[:-2,1:-1])\
            *(dx[1:-1,1:-1]+dx[1:-1,:-2])/dy[:-2,1:-1]
        aN[1:-1,1:-1]=0.5*(2*k[1:-1,1:-1]*k[2:,1:-1])/(k[1:-1,1:-1]+k[2:,1:-1])\
            *(dx[1:-1,1:-1]+dx[1:-1,:-2])/dy[1:-1,1:-1]
        
        # Heat conduction in y direction (Central differences)
        aS[1:-1,1:-1] =0.5*(2*k[1:-1,1:-1]*k[:-2,1:-1])/(k[1:-1,1:-1]+k[:-2,1:-1])\
            *(dx[1:-1,1:-1]+dx[1:-1,:-2])/(dy[:-2,1:-1])
        aN[1:-1,1:-1] =0.5*(2*k[1:-1,1:-1]*k[2:,1:-1])/(k[1:-1,1:-1]+k[2:,1:-1])\
            *(dx[1:-1,1:-1]+dx[1:-1,:-2])/(dy[1:-1,1:-1])
        # Area account for left/right boundary nodes
        aS[1:-1,0]    =0.5*(2*k[1:-1,0]*k[:-2,0])/(k[1:-1,0]+k[:-2,0])\
            *(dx[1:-1,0])/(dy[:-2,0])
        aN[1:-1,0]    =0.5*(2*k[1:-1,0]*k[2:,0])/(k[1:-1,0]+k[2:,0])\
            *(dx[1:-1,0])/(dy[1:-1,0])
        aS[1:-1,-1]   =0.5*(2*k[1:-1,-1]*k[:-2,-1])/(k[1:-1,-1]+k[:-2,-1])\
            *(dx[1:-1,-1])/(dy[:-2,-1])
        aN[1:-1,-1]   =0.5*(2*k[1:-1,-1]*k[2:,-1])/(k[1:-1,-1]+k[2:,-1])\
            *(dx[1:-1,-1])/(dy[1:-1,-1])
        # Forward/backward difference for north/south boundaries
        aN[0,0]       =0.5*(2*k[0,0]*k[1,0])/(k[0,0]+k[1,0])\
            *dx[0,0]/dy[0,0]
        aN[0,1:-1]    =0.5*(2*k[0,1:-1]*k[1,1:-1])/(k[0,1:-1]+k[1,1:-1])\
            *(dx[0,1:-1]+dx[0,:-2])/dy[0,1:-1]
        aN[0,-1]      =0.5*(2*k[0,-1]*k[1,-1])/(k[0,-1]+k[1,-1])\
            *dx[0,-1]/dy[0,-1]
        aS[-1,0]      =0.5*(2*k[-1,0]*k[-2,0])/(k[-1,0]+k[-2,0])\
            *dx[-1,0]/dy[-1,0]
        aS[-1,1:-1]   =0.5*(2*k[-1,1:-1]*k[-2,1:-1])/(k[-1,1:-1]+k[-2,1:-1])\
            *(dx[0,1:-1]+dx[0,:-2])/dy[-1,1:-1]
        aS[-1,-1]     =0.5*(2*k[-1,-1]*k[-2,-1])/(k[-1,-1]+k[-2,-1])\
            *dx[-1,-1]/dy[-1,-1]
        
        return aW,aE,aS,aN
    
    # Bondary condition handler
    def Apply_BCs_Cond(self, E, T_prev, dt, rho, Cv):
        # Left face
        for i in range(len(self.BCs['bc_left'])/3):
            st=self.BCs['bc_left'][2+3*i][0]
            en=self.BCs['bc_left'][2+3*i][1]
            if self.BCs['bc_left'][3*i]=='T':
                E[st:en,0]=self.BCs['bc_left'][1+3*i]*rho[st:en,0]*Cv[st:en,0]*self.Domain.CV_vol()[st:en,0]
                if len(self.BCs['bc_left'])/3-i==1:
                    E[-1,0]=self.BCs['bc_left'][-2]*rho[-1,0]*Cv[-1,0]
            
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
                E[st:en,-1]=self.BCs['bc_right'][1+3*i]*rho[st:en,-1]*Cv[st:en,-1]*self.Domain.CV_vol()[st:en,-1]
                if len(self.BCs['bc_right'])/3-i==1:
                    E[-1,-1]=self.BCs['bc_right'][-2]*rho[-1,-1]*Cv[-1,-1]
            
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
                E[0,st:en]=self.BCs['bc_south'][1+3*i]*rho[0,st:en]*Cv[0,st:en]*self.Domain.CV_vol()[0,st:en]
                if len(self.BCs['bc_south'])/3-i==1:
                    E[0,-1]=self.BCs['bc_south'][-2]*rho[0,-1]*Cv[0,-1]
            
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
                E[-1,st:en]=self.BCs['bc_north'][1+3*i]*rho[-1,st:en]*Cv[-1,st:en]*self.Domain.CV_vol()[-1,st:en]
                if len(self.BCs['bc_north'])/3-i==1:
                    E[-1,-1]=self.BCs['bc_north'][-2]*rho[-1,-1]*Cv[-1,-1]
            
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
    def Advance_Soln_Cond(self, nt, t):
        E_0=self.Domain.E.copy()# Data from previous time step
        
        # Calculate properties
        k, rho, Cv=self.Domain.calcProp()
        
        if self.dt=='None':
            dt=self.getdt(k, rho, Cv)
        else:
            dt=min(self.dt,self.getdt(k, rho, Cv))
            
        if (np.isnan(dt)) or (dt<=0):
            print '*********Diverging time step***********'
            return 1, dt
        print 'Time step %i, Step size=%.7f, Time elapsed=%f;'%(nt+1,dt, t+dt)
        
        # Calculate flux coefficients
        aW,aE,aS,aN=self.get_Coeff(self.dx,self.dy, dt, k, rho, Cv)
        
        T_c=self.Domain.TempFromConserv()
        
        ###################################################################
        # Calculate source and Porous medium terms
        ###################################################################
        # Source terms
        E_unif,E_kim=0,0
        if self.source_unif!='None':
            E_unif      = self.get_source.Source_Uniform(self.source_unif, self.dx, self.dy)
        if self.source_Kim=='True':
            E_kim, deta =self.get_source.Source_Comb_Kim(rho, T_c, self.Domain.eta, self.dx, self.dy, dt)
            
        # Porous medium equations [TO BE CONTINUED]
#        self.Porous_Eqns
        
        ###################################################################
        # Conservation of Energy
        ###################################################################
        # Conduction contribution (2nd order central schemes)
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
        self.Apply_BCs_Cond(self.Domain.E, T_c, dt, rho, Cv)
        
        ###################################################################
        # Conservation of species
        ###################################################################
        # Species generated/destroyed during reaction
        self.Domain.Y_species[:,:,0]+=deta*dt
        self.Domain.Y_species[:,:,1]-=deta*dt
        
        # Species advected from Porous medium equations [TO BE CONTINUED]
        
        
        ###################################################################
        # Divergence/Convergence checks
        ###################################################################
        if (np.isnan(np.amax(self.Domain.E))) \
        or (np.amax(self.Domain.E)>100*np.amax(E_0)) \
        or (np.amin(self.Domain.E)<=0) or (np.amax(self.Domain.eta)>1.0):
            print '**************Divergence detected****************'
            return 1, dt
        else:
            return 0, dt