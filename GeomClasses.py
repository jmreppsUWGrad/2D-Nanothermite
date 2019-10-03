# -*- coding: utf-8 -*-
"""
######################################################
#             2D Heat Conduction Solver              #
#              Created by J. Mark Epps               #
#          Part of Masters Thesis at UW 2018-2020    #
######################################################

This file contains the 2D planar and axisymmetric domain classes:
    -holds conservative variable (energy) at each node
    -holds thermal properties at each node
    -holds x and y coordinate arrays
    -holds dx and dy discretization arrays
    -calculates thermal properties
    -meshing function (biasing feature not functional in solver)
    -function to return temperature given conservative variable (energy)
    -calculate CV 'volume' at each node

Requires:
    -length and width of domain
    -number of nodes across length and width
    -values for thermal properties (depends on what method of calculation)
    
This is intended to be used with 2D Heat Conduction solver with possiblility
of using with 2D Compressible Navier-Stokes solver.

"""
import numpy as np
import string as st
import copy
import MatClasses

class TwoDimDomain():
    def __init__(self, settings, Species, solver, rank):
        
        self.L=settings['Length']
        self.W=settings['Width']
        self.Nx=settings['Nodes_x']
        self.Ny=settings['Nodes_y']
        self.porosity=settings['Porosity']
        self.pore_gas=settings['pore_gas']
        self.type=solver
        self.x=np.zeros(self.Nx)
        self.y=np.zeros(self.Ny)
        self.dx=np.zeros(self.Nx) # NOTE: SIZE MADE TO MATCH REST OF ARRAYS (FOR NOW)
        self.dy=np.zeros(self.Ny) # NOTE: SIZE MADE TO MATCH REST OF ARRAYS (FOR NOW)
        self.rank=rank
        
        # Variables for conservation equations
        self.E=np.zeros((self.Ny, self.Nx)) # Lumped energy
        self.max_iter=settings['Max_iterations']
        self.conv=settings['Convergence']
        
        # Thermal properties
        self.k=settings['k']
        self.rho=settings['rho']
        self.Cv=settings['Cp']
        if type(self.rho) is str and (st.find(self.rho, 'eta')>=0):
            line=st.split(self.rho, ',')
            self.rho0=float(line[1])
            self.rho1=float(line[2])
        if type(self.Cv) is str and (st.find(self.Cv, 'eta')>=0):
            line=st.split(self.Cv, ',')
            self.Cv0=float(line[1])
            self.Cv1=float(line[2])
        if type(self.k) is str and (st.find(self.k, 'eta')>=0):
            line=st.split(self.k, ',')
            self.k0=float(line[1])
            self.k1=float(line[2])
        self.mu=settings['Darcy_mu']
        self.perm=settings['Darcy_perm']
        self.R=settings['gas_constant']
        self.Diff=MatClasses.Diff_Coef()
        self.Cp_calc=MatClasses.Cp()
        self.k_calc=MatClasses.therm_cond()
        
        # Biasing options       
        self.xbias=[settings['bias_type_x'], settings['bias_size_x']]
        self.ybias=[settings['bias_type_y'], settings['bias_size_y']]
        self.isMeshed=False
        
        # MPI information (set by function in main)
        self.proc_left=-1
        self.proc_right=-1
        self.proc_top=-1
        self.proc_bottom=-1
        self.proc_arrang=0 # Array holding process arrangment in domain
        self.proc_row=0 # Row number where rank is in proc_arrang
        
    # Discretize domain and save dx and dy
    def mesh(self):
        # Discretize x
        if self.xbias[0]=='OneWayUp':
            smallest=self.xbias[1]
            self.dx[:-1]=np.linspace(2*self.L/(self.Nx-1)-smallest,smallest,self.Nx-1)
            self.dx[-1]=self.dx[-2]
#            print 'One way biasing in x: smallest element at x=%2f'%self.L
#            print 'Element size range: %2f, %2f'%(smallest, 2*self.L/(self.Nx-1)-smallest)
        elif self.xbias[0]=='OneWayDown':
            smallest=self.xbias[1]
            self.dx[:-1]=np.linspace(smallest,2*self.L/(self.Nx-1)-smallest,self.Nx-1)
            self.dx[-1]=self.dx[-2]
#            print 'One way biasing in x: smallest element at x=0'
#            print 'Element size range: %2f, %2f'%(smallest, 2*self.L/(self.Nx-1)-smallest)
        elif self.xbias[0]=='TwoWayEnd':
            smallest=self.xbias[1]
            self.dx[:int(self.Nx/2)]=np.linspace(smallest,2*self.L/(self.Nx-1)-smallest,(self.Nx-1)/2)
            self.dx[int(self.Nx/2):-1]=np.linspace(2*self.L/(self.Nx-1)-smallest,smallest,(self.Nx-1)/2)
            self.dx[-1]=self.dx[-2]
#            print 'Two way biasing in x: smallest elements at x=0 and %2f'%self.L
#            print 'Element size range: %2f, %2f'%(smallest, 2*self.L/(self.Nx-1)-smallest)
        elif self.xbias[0]=='TwoWayMid':
            smallest=self.xbias[1]
            self.dx[:int(self.Nx/2)]=np.linspace(2*self.L/(self.Nx-1)-smallest,smallest,(self.Nx-1)/2)
            self.dx[int(self.Nx/2):-1]=np.linspace(smallest,2*self.L/(self.Nx-1)-smallest,(self.Nx-1)/2)
            self.dx[-1]=self.dx[-2]
#            print 'Two way biasing in x: smallest elements around x=%2f'%(self.L/2)
#            print 'Element size range: %2f, %2f'%(smallest, 2*self.L/(self.Nx-1)-smallest)
        else:
            self.dx[:]=self.L/(self.Nx-1)
#            print 'No biasing schemes specified in x'
        
        # Discretize y
        if self.ybias[0]=='OneWayUp':
            smallest=self.ybias[1]
            self.dy[:-1]=np.linspace(2*self.W/(self.Ny-1)-smallest,smallest,self.Ny-1)
            self.dy[-1]=self.dy[-2]
#            print 'One way biasing in y: smallest element at y=%2f'%self.W
#            print 'Element size range: %2f, %2f'%(smallest, 2*self.W/(self.Ny-1)-smallest)
        elif self.ybias[0]=='OneWayDown':
            smallest=self.ybias[1]
            self.dy[:-1]=np.linspace(smallest,2*self.W/(self.Ny-1)-smallest,self.Ny-1)
            self.dy[-1]=self.dy[-2]
#            print 'One way biasing in y: smallest element at y=0'
#            print 'Element size range: %2f, %2f'%(smallest, 2*self.W/(self.Ny-1)-smallest)
        elif self.ybias[0]=='TwoWayEnd':
            smallest=self.ybias[1]
            self.dy[:int(self.Ny/2)]=np.linspace(smallest,2*self.W/(self.Ny-1)-smallest,(self.Ny-1)/2)
            self.dy[int(self.Ny/2):-1]=np.linspace(2*self.W/(self.Ny-1)-smallest,smallest,(self.Ny-1)/2)
            self.dy[-1]=self.dy[-2]
#            print 'Two way biasing in y: smallest elements at y=0 and %2f'%self.W
#            print 'Element size range: %2f, %2f'%(smallest, 2*self.W/(self.Ny-1)-smallest)
        elif self.ybias[0]=='TwoWayMid':
            smallest=self.ybias[1]
            self.dy[:int(self.Ny/2)]=np.linspace(2*self.W/(self.Ny-1)-smallest,smallest,(self.Ny-1)/2)
            self.dy[int(self.Ny/2):-1]=np.linspace(smallest,2*self.W/(self.Ny-1)-smallest,(self.Ny-1)/2)
            self.dy[-1]=self.dy[-2]
#            print 'Two way biasing in y: smallest elements around y=%2f'%(self.W/2)
#            print 'Element size range: %2f, %2f'%(smallest, 2*self.W/(self.Ny-1)-smallest)
        else:
            self.dy[:]=self.W/(self.Ny-1)
#            print 'No biasing schemes specified in y'

        for i in range(self.Nx-1):
            self.x[i+1]=self.x[i]+self.dx[i]
        for i in range(self.Ny-1):
            self.y[i+1]=self.y[i]+self.dy[i]
        self.X,self.Y=np.meshgrid(self.x,self.y)
        self.dX,self.dY=np.meshgrid(self.dx,self.dy)
        self.isMeshed=True
    
    # Define other variables for calculations after MPI
    def create_var(self, Species):
        self.eta=np.zeros_like(self.E) # extent of reaction
        self.P=np.zeros_like(self.E) # pressure
        self.T_guess=np.zeros_like(self.E)
        
        # Species
        self.rho_species={}
#        self.mu_species={}
#        self.mv_species={}
        self.Cp_species={}
        self.Cv_species={}
        self.rho_0=np.zeros_like(self.E)
        por=[self.porosity, 1-self.porosity]
        if bool(Species):
            self.species_keys=Species['keys']
            for i in range(len(self.species_keys)):
                self.rho_species[self.species_keys[i]]=np.ones_like(self.E)*Species['Specie_IC'][i]
            self.rho_0+=por[1]*self.rho_species[self.species_keys[1]]
        
    # Calculate and return the dimensions of control volumes
    def CV_dim(self):
        hx=np.zeros_like(self.E)
        hy=np.zeros_like(self.E)
        
        hx[:,1:-1]=0.5*(self.dX[:,1:-1]+self.dX[:,:-2])
        hx[:,0]=0.5*(self.dX[:,0])
        hx[:,-1]=0.5*(self.dX[:,-1])
        
        hy[1:-1,:]=0.5*(self.dY[1:-1,:]+self.dY[:-2,:])
        hy[0,:]=0.5*(self.dY[0,:])
        hy[-1,:]=0.5*(self.dY[-1,:])
        
        return hx,hy
    
    # Calculate temperature dependent properties
    def calcProp(self, T_guess=300, init=False):
        k=np.zeros_like(self.eta)
        rho=np.zeros_like(self.eta)
        rhoC=np.zeros_like(self.eta)
        Cp=np.zeros_like(self.eta)
        
        # Heat capacity (rho*Cv) when species model active
        if bool(self.rho_species):
            if (type(self.Cv) is str) and (st.find(self.Cv, 'eta')>=0):
                Cv=self.eta*self.Cv1+(1-self.eta)*(self.Cv0)
            else:
                Cv=self.Cv
            rhoC=(1-self.porosity)*self.rho_species[self.species_keys[1]]*Cv
            rhoC+=self.porosity*self.rho_species[self.species_keys[0]]*self.Cp_calc.get_Cv(T_guess, self.pore_gas)
            rho=self.rho_species[self.species_keys[1]]
            T=self.E/rhoC
            # Iteratively solve temperature (temperature dependent properties)
#            T_0=np.ones_like(self.eta)
#            T=np.ones_like(self.eta)*T_guess # Initial guess for temperature
#            i=0
#            while np.amax(np.abs(T_0-T)/T)>self.conv and i<self.max_iter:
#                T_0=T.copy()
#                rhoC=(1-self.porosity)*self.rho_species[self.species_keys[1]]*Cv
#                rhoC+=self.porosity*self.rho_species[self.species_keys[0]]*self.Cp_calc.get_Cv(T_guess, self.pore_gas)
#                T=self.E/rhoC
#                i+=1
#                if init:
#                    break            
            
        # Plain heat transfer when species model not active
        else:
            rho[:,:]=self.rho*(1-self.porosity)
            if (type(self.Cv) is str) and (st.find(self.Cv, 'eta')>=0):
                Cv=self.eta*self.Cv1+(1-self.eta)*(self.Cv0)
            else:
                Cv=self.Cv
            rhoC=rho*Cv
            T=self.E/rhoC
            # Iteratively solve temperature (temperature dependent properties)
#            T_0=np.ones_like(self.eta)
#            T=np.ones_like(self.eta)*T_guess # Initial guess for temperature
#            i=0
#            while np.amax(np.abs(T_0-T)/T)>self.conv and i<self.max_iter:
#                T_0=T.copy()
#                rhoC=rho*Cv
#                T=self.E/rhoC
#                i+=1
#                if init:
#                    break 
            
        # Specific heat (Cp) when species model active
        if bool(self.rho_species):
            # Products (only these have gas phases)
            if self.species_keys[0]=='Ar':
                # Argon as only gas specie
                Cp=self.Cp_calc.get_Cp(T, 'Ar')
            else:
                # Special mix of products, no argon or air present
#                Cv_Al2O3=self.Cp_calc.get_Cp(T,'Al2O3')
                Cv_Al2O3=self.Cp_calc.get_Cp(np.ones_like(T)*2327,'Al2O3')
#                Cv_Cu=self.Cp_calc.get_Cp(T,'Cu')
                Cv_Cu=self.Cp_calc.get_Cp(np.ones_like(T)*2843,'Cu')
                
#                Cp=self.rho_species[self.species_keys[0]]*por[0]*(0.351*Cv_Al2O3+0.649*Cv_Cu)/rho
                Cp=(0.351*Cv_Al2O3+0.649*Cv_Cu)
                
        # Thermal conductivity
        if (type(self.k) is str) and (st.find(self.k, 'eta')>=0):
#            k=(self.eta/self.k1+(1-self.eta)/self.k0)**(-1)
            k=self.eta*self.k1+(1-self.eta)*self.k0
#            ks=(self.eta/self.k1+(1-self.eta)/self.k0)**(-1)
#            kf=self.k_calc.get_k(T, self.pore_gas)
#            k[:,:]=ks*(kf/ks)**(self.porosity)
        elif type(self.k) is float:
            k[:,:]=self.k
#            k[:,:]=20
#            k[:,-10:]=80 # Right side
#            k[-10:,:]=150 # Top surface
#            kf=self.k_calc.get_k(T, self.pore_gas)
#            k[:,:]=self.k*(kf/self.k)**(self.porosity)
        
        if init:
            return rhoC
        else:
            return T, k, rho, rhoC, Cp
    