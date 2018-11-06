# -*- coding: utf-8 -*-
"""
Created on Sat Sep 29 13:17:11 2018

@author: Joseph

Solver classes for 2D Heat Conduction. Takes in given object (geometry),
time step and convergence information and alters the object's temperature. 
BCs are applied as appropriate, but must be defined and copied into the 
solver object.

Assumptions:
    -equal discretization spacings in either x or y (want to adjust)
    -constant thermal conductivity

Features:
    -time step based on Fourrier number and discretizations in x and y
    -


"""

import numpy
#import GeomClasses
#import MatClasses
#import CoolProp.CoolProp as CP
#import temporal_schemes

# 1D Solvers (CURRENTLY ONLY FOR CONDUCTION)
class OneDimSolve():
    def __init__(self, geom, timeSize, timeSteps, conv):
        self.Domain=geom # Geometry object
        self.dt=timeSize
        self.Nt=timeSteps
        self.conv=conv
        self.T=self.Domain.T
        self.dx=self.Domain.dx
        self.maxCount=1000
        self.Fo=1.0*self.Domain.mat_prop['k']*self.dt\
        /(self.Domain.mat_prop['rho']*self.Domain.mat_prop['Cp'])
        self.BCs={'BCx1': ('T',600,(0,-1)),\
                 'BCx2': ('T',300,(0,-1)),\
                 'BCy1': ('T',600,(0,-1)),\
                 'BCy2': ('T',300,(0,-1))\
                 }
    
    # Convergence checker
    def CheckConv(self, Tprev, Tnew):
        diff=numpy.sum(numpy.abs(Tnew[:]-Tprev[:]))/numpy.sum(numpy.abs(Tprev[:]))
        print(diff)
        if diff<=self.conv:
            return True
        else:
            return False
    # Solve
    def SolveExpTrans(self):
        Tc=numpy.empty_like(self.T)
        for i in range(self.Nt):
            Tc=self.T.copy()
            self.T[1:-1]=2*self.Fo/(self.dx[:-1]+self.dx[1:])*(Tc[:-2]/self.dx[:-1]+Tc[2:]/self.dx[1:])\
            +(1-2*self.Fo/(self.dx[:-1]+self.dx[1:])*(1/self.dx[:-1]+1/self.dx[1:]))*Tc[1:-1]
        
    def SolveSS(self):
        Tc=numpy.empty_like(self.T)
        count=0
        print 'Residuals:'
        while count<self.maxCount:
            Tc=self.T.copy()
            self.T[1:-1]=(self.dx[1:]*Tc[:-2]+self.dx[:-1]*Tc[2:])\
            /(self.dx[1:]+self.dx[:-1])
            if self.CheckConv(Tc,self.T):
                break

# 2D solver
class TwoDimPlanarSolve():
    def __init__(self, geom_obj, settings, BCs, solver):
        self.Domain=geom_obj # Geometry object
        self.time_scheme=settings['Time_Scheme']
        self.dx,self.dy=numpy.meshgrid(geom_obj.dx,geom_obj.dy)
        self.BCs=BCs
        
        if solver=='Fluid':
            self.CFL=settings['CFL']
            self.gx=settings['Gravity_x']
            self.gy=settings['Gravity_y']
        else:
            self.Fo=settings['Fo']
            self.conv=settings['Convergence']
            self.countmax=settings['Max_iterations']
    
    # Time step check with dx, dy, Fo number
    def getdt(self):
        dt=numpy.zeros_like(self.dx)
        dt[1:-1,1:-1]=0.25*self.Fo*self.Domain.rho*self.Domain.Cv/self.Domain.k*\
            (self.dx[1:-1,1:-1]+self.dx[1:-1,:-2])*(self.dy[1:-1,1:-1]+self.dy[:-2,1:-1])
        dt[0,0]      =0.25*self.Fo*self.Domain.rho*self.Domain.Cv/self.Domain.k*\
            (self.dx[0,0])*(self.dy[0,0])
        dt[0,1:-1]   =0.25*self.Fo*self.Domain.rho*self.Domain.Cv/self.Domain.k*\
            (self.dx[0,1:-1]+self.dx[0,:-2])*(self.dy[0,1:-1])
        dt[1:-1,0]   =0.25*self.Fo*self.Domain.rho*self.Domain.Cv/self.Domain.k*\
            (self.dx[1:-1,0])*(self.dy[1:-1,0]+self.dy[:-2,0])
        dt[0,-1]     =0.25*self.Fo*self.Domain.rho*self.Domain.Cv/self.Domain.k*\
            (self.dx[0,-1])*(self.dy[0,-1])
        dt[-1,0]     =0.25*self.Fo*self.Domain.rho*self.Domain.Cv/self.Domain.k*\
            (self.dx[-1,0])*(self.dy[-1,0])
        dt[-1,1:-1]  =0.25*self.Fo*self.Domain.rho*self.Domain.Cv/self.Domain.k*\
            (self.dx[-1,1:-1]+self.dx[-1,:-2])*(self.dy[-1,1:-1])
        dt[1:-1,-1]   =0.25*self.Fo*self.Domain.rho*self.Domain.Cv/self.Domain.k*\
            (self.dx[1:-1,-1])*(self.dy[1:-1,-1]+self.dy[:-2,-1])
        dt[-1,-1]    =0.25*self.Fo*self.Domain.rho*self.Domain.Cv/self.Domain.k*\
            (self.dx[-1,-1])*(self.dy[-1,-1])
        
        return numpy.amin(dt)

    # Convergence checker
    def CheckConv(self, Tprev, Tnew):
        diff=numpy.sum(numpy.abs(Tnew-Tprev))/numpy.sum(numpy.abs(Tprev))
        print(diff)
        if diff<=self.conv:
            return True
        else:
            return False

    # Heat conduction into/out of each node CV (thermal conductivity not included in calc)
    def get_Cond(self, T, dx, dy):
        qx=numpy.empty_like(T)
        qy=numpy.empty_like(T)
        # Heat conduction in x direction (Central differences)
        qx[1:-1,1:-1] =-0.5*(dy[1:-1,1:-1]+dy[:-2,1:-1])*(T[1:-1,1:-1]-T[1:-1,:-2])/(dx[1:-1,:-2])
        qx[1:-1,1:-1]-=-0.5*(dy[1:-1,1:-1]+dy[:-2,1:-1])*(T[1:-1,2:]-T[1:-1,1:-1])/(dx[1:-1,1:-1])
        # Area account for north/south bondary nodes
        qx[0,1:-1]    =-0.5*(dy[0,1:-1])*(T[0,1:-1]-T[0,:-2])/(dx[0,:-2])
        qx[0,1:-1]   -=-0.5*(dy[0,1:-1])*(T[0,2:]-T[0,1:-1])/(dx[0,1:-1])
        qx[-1,1:-1]   =-0.5*(dy[-1,1:-1])*(T[-1,1:-1]-T[-1,:-2])/(dx[-1,:-2])
        qx[-1,1:-1]  -=-0.5*(dy[-1,1:-1])*(T[-1,2:]-T[-1,1:-1])/(dx[-1,1:-1])
        # Forward/backward difference for left/right boundaries
        qx[0,0]       =0.5*(dy[0,0])*(T[0,1]-T[0,0])/dx[0,0]
        qx[1:-1,0]    =0.5*(dy[1:-1,0]+dy[:-2,0])*(T[1:-1,1]-T[1:-1,0])/dx[1:-1,0]
        qx[-1,0]      =0.5*(dy[-1,0])*(T[-1,1]-T[-1,0])/dx[-1,0]
        qx[0,-1]      =-0.5*(dy[0,-1])*(T[0,-1]-T[0,-2])/dx[0,-1]
        qx[1:-1,-1]   =-0.5*(dy[1:-1,-1]+dy[:-2,-1])*(T[1:-1,-1]-T[1:-1,-2])/dx[1:-1,-1]
        qx[-1,-1]     =-0.5*(dy[-1,-1])*(T[-1,-1]-T[-1,-2])/dx[-1,-1]
        
        # Heat conduction in y direction (Central differences)
        qy[1:-1,1:-1] =-0.5*(dx[1:-1,1:-1]+dx[1:-1,:-2])*(T[1:-1,1:-1]-T[:-2,1:-1])/(dy[:-2,1:-1])
        qy[1:-1,1:-1]-=-0.5*(dx[1:-1,1:-1]+dx[1:-1,:-2])*(T[2:,1:-1]-T[1:-1,1:-1])/(dy[1:-1,1:-1])
        # Area account for left/right boundary nodes
        qy[1:-1,0]    =-0.5*(dx[1:-1,0])*(T[1:-1,0]-T[:-2,0])/(dy[:-2,0])
        qy[1:-1,0]   -=-0.5*(dx[1:-1,0])*(T[2:,0]-T[1:-1,0])/(dy[1:-1,0])
        qy[1:-1,-1]   =-0.5*(dx[1:-1,-1])*(T[1:-1,-1]-T[:-2,-1])/(dy[:-2,-1])
        qy[1:-1,-1]  -=-0.5*(dx[1:-1,-1])*(T[2:,-1]-T[1:-1,-1])/(dy[1:-1,-1])
        # Forward/backward difference for north/south boundaries
        qy[0,0]       =0.5*dx[0,0]*(T[0,1]-T[0,0])/dy[0,0]
        qy[0,1:-1]    =0.5*(dx[0,1:-1]+dx[0,:-2])*(T[1,1:-1]-T[0,1:-1])/dy[0,1:-1]
        qy[0,-1]      =0.5*dx[0,-1]*(T[0,-1]-T[0,-2])/dy[0,-1]
        qy[-1,0]      =-0.5*dx[-1,0]*(T[-1,0]-T[-2,0])/dy[-1,0]
        qy[-1,1:-1]   =-0.5*(dx[0,1:-1]+dx[0,:-2])*(T[-1,1:-1]-T[-2,1:-1])/dy[-1,1:-1]
        qy[-1,-1]     =-0.5*dx[-1,-1]*(T[-1,-1]-T[-2,-1])/dy[-1,-1]
        
        return qx+qy
    
    # Bondary condition handler
    def Apply_BCs_Cond(self, T, Fo, T_prev):
#        BC1x,BC1y='T','T'# BC types at corner 1
#        BC2x,BC2y='T','T'# BC types at corner 2
#        BC3x,BC3y='T','T'# BC types at corner 3
#        BC4x,BC4y='T','T'# BC types at corner 4
        
        # Left face
        for i in range(len(self.BCs['bc_left'])/3):
            st=self.BCs['bc_left'][2+3*i][0]
            en=self.BCs['bc_left'][2+3*i][1]
            if self.BCs['bc_left'][3*i]=='T':
                T[st:en,0]=self.BCs['bc_left'][1+3*i]
                if len(self.BCs['bc_left'])/3-i==1:
                    T[-1,0]=self.BCs['bc_left'][-2]
            
            else:
                if self.BCs['bc_left'][3*i]=='F':
                    q=self.BCs['bc_left'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_left'][1+3*i][0]*self.BCs['bc_left'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_left'][1+3*i][0]*T_prev[st:en,0] # h*Tij
                
                T[st:en,0]+=Fo[st:en,0]*(Bi+q)*self.dy[st:en,0]/self.Domain.k
                if len(self.BCs['bc_left'])/3-i==1:
                    if self.BCs['bc_left'][3*i]=='C':
                        Bi=-self.BCs['bc_left'][1+3*i][0]*T_prev[-1,0] # h*Tij
                    T[-1,0]+=Fo[-1,0]*(Bi+q)*self.dy[-1,0]/self.Domain.k
        
        # Right face
        for i in range(len(self.BCs['bc_right'])/3):
            st=self.BCs['bc_right'][2+3*i][0]
            en=self.BCs['bc_right'][2+3*i][1]
            if self.BCs['bc_right'][3*i]=='T':
                T[st:en,-1]=self.BCs['bc_right'][1+3*i]
                if len(self.BCs['bc_right'])/3-i==1:
                    T[-1,-1]=self.BCs['bc_right'][-2]
            
            else:
                if self.BCs['bc_right'][3*i]=='F':
                    q=self.BCs['bc_right'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_right'][1+3*i][0]*self.BCs['bc_right'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_right'][1+3*i][0]*T_prev[st:en,-1] # h*Tij
                
                T[st:en,-1]+=Fo[st:en,-1]*(Bi+q)*self.dy[st:en,-1]/self.Domain.k
                if len(self.BCs['bc_right'])/3-i==1:
                    if self.BCs['bc_right'][3*i]=='C':
                        Bi=-self.BCs['bc_right'][1+3*i][0]*T_prev[-1,-1] # h*Tij
                    T[-1,-1]+=Fo[-1,-1]*(Bi+q)*self.dy[-1,-1]/self.Domain.k
        
        # South face
        for i in range(len(self.BCs['bc_south'])/3):
            st=self.BCs['bc_south'][2+3*i][0]
            en=self.BCs['bc_south'][2+3*i][1]
            if self.BCs['bc_south'][3*i]=='T':
                T[0,st:en]=self.BCs['bc_south'][1+3*i]
                if len(self.BCs['bc_south'])/3-i==1:
                    T[0,-1]=self.BCs['bc_south'][-2]
            
            else:
                if self.BCs['bc_south'][3*i]=='F':
                    q=self.BCs['bc_south'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_south'][1+3*i][0]*self.BCs['bc_south'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_south'][1+3*i][0]*T_prev[0,st:en] # h*Tij
                
                T[0,st:en]+=Fo[0,st:en]*(Bi+q)*self.dx[0,st:en]/self.Domain.k
                if len(self.BCs['bc_south'])/3-i==1:
                    if self.BCs['bc_south'][3*i]=='C':
                        Bi=-self.BCs['bc_south'][1+3*i][0]*T_prev[0,-1] # h*Tij
                    T[0,-1]+=Fo[0,-1]*(Bi+q)*self.dx[0,-1]/self.Domain.k
                    
        # North face
        for i in range(len(self.BCs['bc_north'])/3):
            st=self.BCs['bc_north'][2+3*i][0]
            en=self.BCs['bc_north'][2+3*i][1]
            if self.BCs['bc_north'][3*i]=='T':
                T[-1,st:en]=self.BCs['bc_north'][1+3*i]
                if len(self.BCs['bc_north'])/3-i==1:
                    T[-1,-1]=self.BCs['bc_north'][-2]
            
            else:
                if self.BCs['bc_north'][3*i]=='F':
                    q=self.BCs['bc_north'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_north'][1+3*i][0]*self.BCs['bc_north'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_north'][1+3*i][0]*T_prev[-1,st:en] # h*Tij
                
                T[-1,st:en]+=Fo[-1,st:en]*(Bi+q)*self.dx[-1,st:en]/self.Domain.k
                if len(self.BCs['bc_north'])/3-i==1:
                    if self.BCs['bc_north'][3*i]=='C':
                        Bi=-self.BCs['bc_north'][1+3*i][0]*T_prev[-1,-1] # h*Tij
                    T[-1,-1]+=Fo[-1,-1]*(Bi+q)*self.dx[-1,-1]/self.Domain.k
    
    # Calculate source term for combustion (Cantera module)
    def Source_Comb(self, T, dx, dy):
        # Terms to include after modelling combustion:
            #divide by k
            #resulting expression should have dimensions of K
            #Be sure to account for area of CV properly
        
        
        q_vol=0 # Volumetric heat generation rate
        q_source=numpy.zeros_like(dx)
        
        q_source[1:-1,1:-1]=q_vol/self.Domain.k*\
            0.25*(self.dx[1:-1,1:-1]+self.dx[1:-1,:-2])*(self.dy[1:-1,1:-1]+self.dy[:-2,1:-1])
        q_source[0,0]      =q_vol/self.Domain.k*\
            0.25*(self.dx[0,0])*(self.dy[0,0])
        q_source[0,1:-1]   =q_vol/self.Domain.k*\
            0.25*(self.dx[0,1:-1]+self.dx[0,:-2])*(self.dy[0,1:-1])
        q_source[1:-1,0]   =q_vol/self.Domain.k*\
            0.25*(self.dx[1:-1,0])*(self.dy[1:-1,0]+self.dy[:-2,0])
        q_source[0,-1]     =q_vol/self.Domain.k*\
            0.25*(self.dx[0,-1])*(self.dy[0,-1])
        q_source[-1,0]     =q_vol/self.Domain.k*\
            0.25*(self.dx[-1,0])*(self.dy[-1,0])
        q_source[-1,1:-1]  =q_vol/self.Domain.k*\
            0.25*(self.dx[-1,1:-1]+self.dx[-1,:-2])*(self.dy[-1,1:-1])
        q_source[1:-1,-1]  =q_vol/self.Domain.k*\
            0.25*(self.dx[1:-1,-1])*(self.dy[1:-1,-1]+self.dy[:-2,-1])
        q_source[-1,-1]    =q_vol/self.Domain.k*\
            0.25*(self.dx[-1,-1])*(self.dy[-1,-1])
        
        return q_source
    
    # Main solver (1 time step)
    def Advance_Soln_Cond(self):
        T_c=self.Domain.T.copy()
        T_0=self.Domain.T.copy()
        
        dt=self.getdt()
        if (numpy.isnan(dt)) or (dt<=0):
            print '*********Diverging time step***********'
            return 1
        Fo=numpy.zeros_like(T_c)
        Fo[1:-1,1:-1]=self.Domain.k*dt/(self.Domain.rho*self.Domain.Cv*\
            0.25*(self.dx[1:-1,1:-1]+self.dx[1:-1,:-2])*(self.dy[1:-1,1:-1]+self.dy[:-2,1:-1]))
        Fo[0,0]      =self.Domain.k*dt/(self.Domain.rho*self.Domain.Cv*\
            0.25*(self.dx[0,0])*(self.dy[0,0]))
        Fo[0,1:-1]   =self.Domain.k*dt/(self.Domain.rho*self.Domain.Cv*\
            0.25*(self.dx[0,1:-1]+self.dx[0,:-2])*(self.dy[0,1:-1]))
        Fo[1:-1,0]   =self.Domain.k*dt/(self.Domain.rho*self.Domain.Cv*\
            0.25*(self.dx[1:-1,0])*(self.dy[1:-1,0]+self.dy[:-2,0]))
        Fo[0,-1]     =self.Domain.k*dt/(self.Domain.rho*self.Domain.Cv*\
            0.25*(self.dx[0,-1])*(self.dy[0,-1]))
        Fo[-1,0]     =self.Domain.k*dt/(self.Domain.rho*self.Domain.Cv*\
            0.25*(self.dx[-1,0])*(self.dy[-1,0]))
        Fo[-1,1:-1]  =self.Domain.k*dt/(self.Domain.rho*self.Domain.Cv*\
            0.25*(self.dx[-1,1:-1]+self.dx[-1,:-2])*(self.dy[-1,1:-1]))
        Fo[1:-1,-1]   =self.Domain.k*dt/(self.Domain.rho*self.Domain.Cv*\
            0.25*(self.dx[1:-1,-1])*(self.dy[1:-1,-1]+self.dy[:-2,-1]))
        Fo[-1,-1]    =self.Domain.k*dt/(self.Domain.rho*self.Domain.Cv*\
            0.25*(self.dx[-1,-1])*(self.dy[-1,-1]))
        
        print 'Time step size: %.7f'%dt
        count=0
        while (count<self.countmax):
            ###################################################################
            # Temperature (2nd order central schemes)
            ###################################################################
            dTdt =self.get_Cond(T_c, self.dx, self.dy)
            dTdt+=self.Source_Comb(T_c, self.dx, self.dy)
            
            self.Domain.T = T_0 + Fo*dTdt
            
            ###################################################################
            # Apply boundary conditions
            ###################################################################
            self.Apply_BCs_Cond(self.Domain.T, Fo, T_c)
            
            ###################################################################
            # Divergence/Convergence checks
            ###################################################################
            if (numpy.isnan(numpy.amax(self.Domain.T))) \
            or (numpy.amax(self.Domain.T)>100*numpy.amax(T_0)) \
            or (numpy.amin(self.Domain.T)<=0):
                print '**************Divergence detected****************'
                return 1
            
            # Break while loop if converged OR is explicit solve
            if (self.time_scheme=='Explicit'):
                break
            elif (self.CheckConv(T_c, self.Domain.T)):
                break
            count+=1
            T_c=self.Domain.T.copy()
        ###################################################################
        # Output data to file?????
        ###################################################################
        
        
        if count==self.countmax:
            print '*************No convergence reached*****************'
            return 1
        else:
            return 0