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
from MatClasses import Cp, therm_cond

class TwoDimDomain():
    def __init__(self, settings, Species, solver, rank):
        
        self.L=settings['Length']
        self.W=settings['Width']
        self.Nx=settings['Nodes_x']
        self.Ny=settings['Nodes_y']
        self.x=np.zeros(self.Nx)
        self.y=np.zeros(self.Ny)
        self.dx=np.zeros(self.Nx) # NOTE: SIZE MADE TO MATCH REST OF ARRAYS (FOR NOW)
        self.dy=np.zeros(self.Ny) # NOTE: SIZE MADE TO MATCH REST OF ARRAYS (FOR NOW)
        self.model=settings['Model']
        self.type=solver
        self.porosity_0=settings['Porosity']
        self.rank=rank
        
        # Variables for conservation equations
        self.E=np.zeros((self.Ny, self.Nx)) # Lumped energy
        self.max_iter=settings['Max_iterations']
        self.conv=settings['Convergence']
        
        # Thermal properties
        self.k=settings['k_s']
        self.k_mode=settings['k_model']
        self.rho=settings['rho_IC']
        self.Cv=settings['Cv_s']
        
        # Process the thermal properties options
        if (type(self.k) is str):
            self.k=st.split(self.k, ',')
        if (type(self.Cv) is str):
            self.Cv=st.split(self.Cv, ',')
        
        # Species model options
        if self.model=='Species':
            self.species_keys=['g','s']
            self.rho=st.split(self.rho, ',')
            self.mu=settings['Darcy_mu']
            self.part_diam=settings['Carmen_diam']
            self.kozeny=settings['Kozeny_const']
            self.R=settings['gas_constant']
            self.Cv_g=Species['Cv_g']
            self.Cp_g=Species['Cp_g']
            self.k_g=Species['k_g']
            # Process thermal properties options
            if (type(self.Cv_g) is str):
                self.Cv_g=st.split(self.Cv_g, ',')
            if (type(self.Cp_g) is str):
                self.Cp_g=st.split(self.Cp_g, ',')
            if (type(self.k_g) is str):
                self.k_g=st.split(self.k_g, ',')
        
        # Object declarations for property calculations
        self.Cp_calc=Cp()
        self.k_calc=therm_cond()
        
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
        self.T_guess=np.ones_like(self.E)
        self.porosity=np.ones_like(self.E)*self.porosity_0
        
        # Species
        self.rho_species={}
        try:
            self.rho_0=np.ones_like(self.E)*self.rho*(1-self.porosity_0)
        except:
            self.rho_0=np.zeros_like(self.E)
        por=[self.porosity, 1-self.porosity]
        if self.model=='Species':
            for i in range(len(self.species_keys)):
                self.rho_species[self.species_keys[i]]=np.ones_like(self.E)\
                    *float(self.rho[i])*por[i]
            self.rho_0=self.rho_species[self.species_keys[1]]
            self.perm=self.porosity**3*self.part_diam**2\
                /(self.kozeny*(1-self.porosity)**2)
        
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
        Cv=np.zeros_like(self.eta)
        Cp=np.zeros_like(self.eta)
        
        ##########################################################################
        # Specific heat of solid phase
        ##########################################################################
        
        if (type(self.Cv) is list) and (self.Cv[0]=='eta'):
            Cv=self.eta*float(self.Cv[2])+(1-self.eta)*float(self.Cv[1])
        # Solid phase (temperature dependent for given element)
        elif (type(self.Cv) is list) and (self.Cv[1]=='Temp'):
            # Constant temperature value
            if len(self.Cv)>2:
                Cv=self.Cp_calc.get_Cv(np.ones_like(self.E)*float(self.Cv[2]), self.Cv[0])
            # Temperature dependent
            else:
                Cv=self.Cp_calc.get_Cv(T_guess, self.Cv[0])
        # Solid phase (constant)
        else:
            Cv[:,:]=self.Cv
        
        ##########################################################################
        # Thermal conductivity of solid phase (either model)
        ##########################################################################
        
        # Solid phase (eta dependent)
        if (type(self.k) is list) and (self.k[0]=='eta'):
            k=self.eta*float(self.k[2])+(1-self.eta)*float(self.k[1])
        
        # Solid phase (temperature dependent for given element)
        elif (type(self.k) is list) and (self.k[1]=='Temp'):
            # Constant temperature value
            if len(self.k)>2:
                k=self.k_calc.get_k(np.ones_like(self.E)*float(self.k[2]), self.k[0])
            # Temperature dependent
            else:
                k=self.k_calc.get_k(T_guess, self.k[0])
        
        # Solid phase (constant)
        else:
            k[:,:]=self.k
        
        ##########################################################################
        #  When species model is active
        ##########################################################################
        
        if self.model=='Species':
            k_g=np.zeros_like(self.eta)
            # Changing porosity/permeability
#            self.porosity=self.porosity_0+\
#                (1-self.rho_species[self.species_keys[1]]/self.rho_0)*(1-self.porosity_0)
#            self.perm=self.porosity**3*self.part_diam**2\
#                /(self.kozeny*(1-self.porosity)**2)
            
            # Heat capacity of Solid phase
            rhoC=self.rho_species[self.species_keys[1]]*Cv
#            rhoC=self.rho*(1-self.porosity)*Cv # REPLICATE CASE 10 (CASE 10d,e)
            
            ##########################################################################
            #####  Heat capacity of Gas phase
            ##########################################################################
            
            if (type(self.Cv_g) is list) and (self.Cv_g[0]=='eta'):
                Cv=self.eta*float(self.Cv_g[2])+(1-self.eta)*float(self.Cv_g[1])
            
            # Gas phase (temperature dependent for given element)
            elif (type(self.Cv_g) is list) and (self.Cv_g[1]=='Temp'):
                # Constant temperature value
                if len(self.Cv_g)>2:
                    Cv=self.Cp_calc.get_Cv(np.ones_like(self.E)*float(self.Cv_g[2]), self.Cv_g[0])
                # Temperature dependent
                else:
                    T_0=np.ones_like(self.eta)
                    rhoc=rhoC.copy()
                    T=np.ones_like(self.eta)*T_guess # Initial guess for temperature
                    i=0
                    while np.amax(np.abs(T_0-T)/T)>self.conv and i<self.max_iter:
                        T_0=T.copy()
                        Cv=self.Cp_calc.get_Cv(T, self.Cv_g[0])
                        rhoc=rhoC+self.rho_species[self.species_keys[0]]*Cv
                        T=self.E/rhoc
                        i+=1
                        if init:
                            break
                    if i>=self.max_iter:
                        Cv=-10**9
                        print('***** Unable to get converging temperature')
            
            # Gas phase (constant)
            else:
                Cv[:,:]=self.Cv_g
            
            # Temperature calculation
            rhoC+=self.rho_species[self.species_keys[0]]*Cv
            T=self.E/rhoC
            self.T_guess=T
            
            ##########################################################################
            #####  Specific heat (Cp) of Gas phase
            ##########################################################################
            # eta dependent
            if (type(self.Cp_g) is list) and (self.Cp_g[0]=='eta'):
                Cp=self.eta*float(self.Cp_g[2])+(1-self.eta)*float(self.Cp_g[1])
            
            # temperature dependent for given element
            elif (type(self.Cp_g) is list) and (self.Cp_g[1]=='Temp'):
                # Constant temperature value
                if len(self.Cp_g)>2:
                    Cp=self.Cp_calc.get_Cp(np.ones_like(self.E)*float(self.Cp_g[2]), self.Cp_g[0])
                # Temperature dependent
                else:
                    Cp=self.Cp_calc.get_Cp(T_guess, self.Cp_g[0])
            
            # constant
            else:
                Cp[:,:]=self.Cp_g
            
            ##########################################################################
            ##### Thermal conductivity of gas phase
            ##########################################################################
            # eta dependent
            if (type(self.k_g) is list) and (self.k_g[0]=='eta'):
                k_g=self.eta*float(self.k_g[2])+(1-self.eta)*float(self.k_g[1])
            
            # temperature dependent for given element
            elif (type(self.k_g) is list) and (self.k_g[1]=='Temp'):
                # Constant temperature value
                if len(self.k_g)>2:
                    k_g=self.k_calc.get_k(np.ones_like(self.E)*float(self.k_g[2]), self.k_g[0])
                # Temperature dependent
                else:
                    k_g=self.k_calc.get_k(T_guess, self.k_g[0])
            
            # constant
            else:
                k_g[:,:]=self.k_g
            
            ##########################################################################
            ##### Thermal conductivity models
            ##########################################################################
            
            if self.k_mode=='Parallel':
                k=self.porosity*k_g+(1-self.porosity)*k
            elif self.k_mode=='Geometric':
                k=k*(k_g/k)**(self.porosity)
            elif self.k_mode=='Series':
                k=(self.porosity/k_g+(1-self.porosity)/k)**(-1)
        
        ##########################################################################
        #  Plain heat transfer model
        ##########################################################################
        
        else:
            rho[:,:]=self.rho*(1-self.porosity)
            
            rhoC=rho*Cv
            T=self.E/rhoC
            self.T_guess=T
        
        if init:
            return rhoC
        else:
            return T, k, rhoC, Cp
    