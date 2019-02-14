# -*- coding: utf-8 -*-
"""
Created on Sat Sep 22 23:14:27 2018

This class is a 2D planar domain for solving temp, vel or density.

This is intended to be used with either 2D Compressible Navier-Stokes solver 
or 2D Heat Conduction solvers.

Requires:
    -length and width of domain
    -number of nodes across length and width
    
Features:
    -Linear biasing one way or two way based on specified smallest element size
    -Creates distance between nodes arrays (dx and dy)

Desired:
    -store CV dimensions (dx and dy) for each point; makes for equal sized
    arrays with temp/vel/press/density

Features:
    -

@author: Joseph
"""
import numpy as np
import string as st

class TwoDimPlanar:
    def __init__(self, settings, solver):
        
        self.L=settings['Length']
        self.W=settings['Width']
        self.Nx=settings['Nodes_x']
        self.Ny=settings['Nodes_y']
        self.x=np.zeros(self.Nx)
        self.y=np.zeros(self.Ny)
        self.dx=np.zeros(self.Nx) # NOTE: SIZE MADE TO MATCH REST OF ARRAYS (FOR NOW)
        self.dy=np.zeros(self.Ny) # NOTE: SIZE MADE TO MATCH REST OF ARRAYS (FOR NOW)
        self.k=settings['k']
        # Fluid solver
        if solver=='Fluid':
            self.fluid=settings['Fluid']
            self.gamma=settings['gamma']
            self.mu=settings['mu']
            self.R=settings['R']
            # Setup variable arrays
            self.tau11=np.zeros((self.Ny, self.Nx)) # Shear stress arrays
            self.tau12=np.zeros((self.Ny, self.Nx))
            self.tau22=np.zeros((self.Ny, self.Nx))
            self.rhoE=np.zeros((self.Ny, self.Nx)) # Conservative arrays
            self.rhou=np.zeros((self.Ny,self.Nx))
            self.rhov=np.zeros((self.Ny,self.Nx))
            self.rho=np.zeros((self.Ny,self.Nx)) # Primitive arrays
            
            self.Cv=self.R/(self.gamma-1)
        else:
            self.E=np.zeros((self.Ny, self.Nx))
            self.eta=np.zeros((self.Ny, self.Nx))
            self.rho=settings['rho']
            self.Cv=settings['Cp']
#            self.Y_species=np.zeros((self.Ny, self.Nx, 15)) # species array
#            self.P=np.zeros((self.Ny, self.Nx))
            if type(self.rho) is str:
                line=st.split(self.rho, ',')
                self.rho0=float(line[1])
                self.rho1=float(line[2])
            if type(self.Cv) is str:
                line=st.split(self.Cv, ',')
                self.Cv0=float(line[1])
                self.Cv1=float(line[2])
            if type(self.k) is str:
                line=st.split(self.k, ',')
                self.k0=float(line[1])
                self.k1=float(line[2])
        
        # Biasing options       
        self.xbias=[settings['bias_type_x'], settings['bias_size_x']]
        self.ybias=[settings['bias_type_y'], settings['bias_size_y']]
        self.isMeshed=False
        # Other useful calculations (put elsewhere??)
        
        
    # Discretize domain and save dx and dy
    def mesh(self):
        # Discretize x
        if self.xbias[0]=='OneWayUp':
            smallest=self.xbias[1]
            self.dx[:-1]=np.linspace(2*self.L/(self.Nx-1)-smallest,smallest,self.Nx-1)
            self.dx[-1]=self.dx[-2]
            print 'One way biasing in x: smallest element at x=%2f'%self.L
            print 'Element size range: %2f, %2f'%(smallest, 2*self.L/(self.Nx-1)-smallest)
        elif self.xbias[0]=='OneWayDown':
            smallest=self.xbias[1]
            self.dx[:-1]=np.linspace(smallest,2*self.L/(self.Nx-1)-smallest,self.Nx-1)
            self.dx[-1]=self.dx[-2]
            print 'One way biasing in x: smallest element at x=0'
            print 'Element size range: %2f, %2f'%(smallest, 2*self.L/(self.Nx-1)-smallest)
        elif self.xbias[0]=='TwoWayEnd':
            smallest=self.xbias[1]
            self.dx[:int(self.Nx/2)]=np.linspace(smallest,2*self.L/(self.Nx-1)-smallest,(self.Nx-1)/2)
            self.dx[int(self.Nx/2):-1]=np.linspace(2*self.L/(self.Nx-1)-smallest,smallest,(self.Nx-1)/2)
            self.dx[-1]=self.dx[-2]
            print 'Two way biasing in x: smallest elements at x=0 and %2f'%self.L
            print 'Element size range: %2f, %2f'%(smallest, 2*self.L/(self.Nx-1)-smallest)
        elif self.xbias[0]=='TwoWayMid':
            smallest=self.xbias[1]
            self.dx[:int(self.Nx/2)]=np.linspace(2*self.L/(self.Nx-1)-smallest,smallest,(self.Nx-1)/2)
            self.dx[int(self.Nx/2):-1]=np.linspace(smallest,2*self.L/(self.Nx-1)-smallest,(self.Nx-1)/2)
            self.dx[-1]=self.dx[-2]
            print 'Two way biasing in x: smallest elements around x=%2f'%(self.L/2)
            print 'Element size range: %2f, %2f'%(smallest, 2*self.L/(self.Nx-1)-smallest)
        else:
            self.dx[:]=self.L/(self.Nx-1)
            print 'No biasing schemes specified in x'
        
        # Discretize y
        if self.ybias[0]=='OneWayUp':
            smallest=self.ybias[1]
            self.dy[:-1]=np.linspace(2*self.W/(self.Ny-1)-smallest,smallest,self.Ny-1)
            self.dy[-1]=self.dy[-2]
            print 'One way biasing in y: smallest element at y=%2f'%self.W
            print 'Element size range: %2f, %2f'%(smallest, 2*self.W/(self.Ny-1)-smallest)
        elif self.ybias[0]=='OneWayDown':
            smallest=self.ybias[1]
            self.dy[:-1]=np.linspace(smallest,2*self.W/(self.Ny-1)-smallest,self.Ny-1)
            self.dy[-1]=self.dy[-2]
            print 'One way biasing in y: smallest element at y=0'
            print 'Element size range: %2f, %2f'%(smallest, 2*self.W/(self.Ny-1)-smallest)
        elif self.ybias[0]=='TwoWayEnd':
            smallest=self.ybias[1]
            self.dy[:int(self.Ny/2)]=np.linspace(smallest,2*self.W/(self.Ny-1)-smallest,(self.Ny-1)/2)
            self.dy[int(self.Ny/2):-1]=np.linspace(2*self.W/(self.Ny-1)-smallest,smallest,(self.Ny-1)/2)
            self.dy[-1]=self.dy[-2]
            print 'Two way biasing in y: smallest elements at y=0 and %2f'%self.W
            print 'Element size range: %2f, %2f'%(smallest, 2*self.W/(self.Ny-1)-smallest)
        elif self.ybias[0]=='TwoWayMid':
            smallest=self.ybias[1]
            self.dy[:int(self.Ny/2)]=np.linspace(2*self.W/(self.Ny-1)-smallest,smallest,(self.Ny-1)/2)
            self.dy[int(self.Ny/2):-1]=np.linspace(smallest,2*self.W/(self.Ny-1)-smallest,(self.Ny-1)/2)
            self.dy[-1]=self.dy[-2]
            print 'Two way biasing in y: smallest elements around y=%2f'%(self.W/2)
            print 'Element size range: %2f, %2f'%(smallest, 2*self.W/(self.Ny-1)-smallest)
        else:
            self.dy[:]=self.W/(self.Ny-1)
            print 'No biasing schemes specified in y'

        for i in range(self.Nx-1):
            self.x[i+1]=self.x[i]+self.dx[i]
        for i in range(self.Ny-1):
            self.y[i+1]=self.y[i]+self.dy[i]
        self.X,self.Y=np.meshgrid(self.x,self.y)
        
        self.isMeshed=True
    
    def primitiveFromConserv(self, rho, rhou, rhov, rhoE):
        u=rhou/rho
        v=rhov/rho
        T=(rhoE/rho-0.5*(u**2+v**2))/self.Cv
        
        # Ideal gas law assumed
        p=rho*self.R*T
        
        return u,v,p,T
    
    # Calculate and return volume of each node
    def CV_vol(self):
        v=np.zeros_like(self.eta)
        dx,dy=np.meshgrid(self.dx, self.dy)
        v[1:-1,1:-1]=0.25*(dx[1:-1,1:-1]+dx[1:-1,:-2])*(dy[1:-1,1:-1]+dy[:-2,1:-1])
        v[0,0]      =0.25*(dx[0,0])*(dy[0,0])
        v[0,1:-1]   =0.25*(dx[0,1:-1]+dx[0,:-2])*(dy[0,1:-1])
        v[1:-1,0]   =0.25*(dx[1:-1,0])*(dy[1:-1,0]+dy[:-2,0])
        v[0,-1]     =0.25*(dx[0,-1])*(dy[0,-1])
        v[-1,0]     =0.25*(dx[-1,0])*(dy[-1,0])
        v[-1,1:-1]  =0.25*(dx[-1,1:-1]+dx[-1,:-2])*(dy[-1,1:-1])
        v[1:-1,-1]   =0.25*(dx[1:-1,-1])*(dy[1:-1,-1]+dy[:-2,-1])
        v[-1,-1]    =0.25*(dx[-1,-1])*(dy[-1,-1])
        return v
    
    # Calculate temperature dependent properties
    def calcProp(self):
        k=np.zeros_like(self.eta)
        rho=np.zeros_like(self.eta)
        Cv=np.zeros_like(self.eta)
        
        # Calculate properties based on eta or constant
        if type(self.k) is str:
            k=(self.eta/self.k1+(1-self.eta)/self.k0)**(-1)
        else:
            k[:,:]=self.k
        if type(self.Cv) is str:
            Cv=self.eta*self.Cv1+(1-self.eta)*self.Cv0
        else:
            Cv[:,:]=self.Cv
        if type(self.rho) is str:
            rho=self.eta*self.rho1+(1-self.eta)*self.rho0
        else:
            rho[:,:]=self.rho
                
        return k, rho, Cv
    
    # Calculate temperature from energy
    def TempFromConserv(self):
        k,rho,Cv=self.calcProp()
        return self.E/Cv/rho/self.CV_vol()
    # Check everything before solving
    def IsReadyToSolve(self):
        if self.isMeshed:
            return True
        else:
            return False
