# -*- coding: utf-8 -*-
"""
######################################################
#             2D Heat Conduction Solver              #
#              Created by J. Mark Epps               #
#          Part of Masters Thesis at UW 2018-2020    #
######################################################

This file contains the Combustion source term classes:
    -Uniform heat source; returns energy generated at each node
    -Combustion source term from Kim paper; returns energy generated
    
Features of Source_Kim:
    -Activation energy, pre-exponential factor, enthalpy of combustion
    -Enthalpy of combustion can be density or volume based (input file)

Notes on implementing Cantera:
    -sol=ct.Solution('___.cti') -> define solution mechanisms?
    -ct.SolutionArray() -> able to define an array of states
    https://cantera.org/documentation/docs-2.4/sphinx/html/cython/importing.html#representing-multiple-states
    
"""

import numpy as np
import string as st
#import cantera as ct

class Source_terms():
    def __init__(self, Ea, A0, dH):
        self.R=8.314 # J/mol/K
        self.Ea=Ea # J/mol
        self.A0=A0
        self.dH=st.split(dH, ',')
        self.dH[1]=float(self.dH[1])
        self.n=0.2 # Temperature exponent
        
    # Uniform volumetric generation
    def Source_Uniform(self, Q, V):
        
        return Q*V
    
    # Calculate source term for combustion based on
    # K. Kim, "Computational Modeling of Combustion Wave in Nanoscale Thermite Reaction",
    # Int. J of Energy and Power engineering, vol.8, no.7, pp. 612-615, 2014.
    def Source_Comb_Kim(self, rho, T, eta, vol, dt):
        detadt=self.A0*(1-eta)*np.exp(-self.Ea/self.R/T)
        eta+=dt*detadt
        
        # Clipping to 0
#        eta[eta<10**(-10)]=0
        
        if st.find(self.dH[0], 'vol')>=0:
            return vol*self.dH[1]*detadt, detadt
        else:
            return rho*vol*self.dH[1]*detadt, detadt
    
    # Source term for combustion based on
    # Umbrajkar, S et al., "Exothermic reactions in Al-CuO nanocomposites",
    # Thermochimica Acta, vol.451, pp. 34-43, 2006.
    def Source_Comb_Umbrajkar(self, rho, T, eta, vol, dt):
        # First temp range
        A=10**(6.68)
        n=0.6
        Ea=78000
        deta1=A*n*(eta-1)*np.log((1-eta)**(1-1/n))*np.exp(-Ea/8.314/T)
        # Second temp range
        A=10**(5.15)
        n=3.9
        Ea=79000
        deta2=A*(1-eta)**n*np.exp(-Ea/8.314/T)
        # Third temp range
        A=10**(5.03)
        n=2.6
        Ea=102000
        deta3=A*(1-eta)**n*np.exp(-Ea/8.314/T)
        # Fourth temp range
        A=10**(13.3)
        n=0.75
        Ea=266000
        deta4=A*n*(eta-1)*np.log((1-eta)**(1-1/n))*np.exp(-Ea/8.314/T)
        
#        deta4=self.A0*(1-eta)*np.exp(-self.Ea/self.R/T)
        
        detadt=deta1+deta2+deta3+deta4
        eta+=dt*detadt
        
        # Clipping to 0
        eta[eta<10**(-5)]=0
        
        if st.find(self.dH[0], 'vol')>=0:
            return vol*self.dH[1]*detadt, detadt
        else:
            return rho*vol*self.dH[1]*detadt, detadt
    
    # Calculate source term for combustion (NEEDS MODIFYING)
    def Source_Comb(self, T, y_species, dx, dy):
        # Terms to include after modelling combustion:
            #resulting expression should have dimensions of W/m
            #Be sure to account for area of CV properly
        at=np.zeros_like(dx)
        q_source=np.zeros_like(dx)
        
        # CV dimensions
        at[1:-1,1:-1]=0.25*(dx[1:-1,1:-1]+dx[1:-1,:-2])*(dy[1:-1,1:-1]+dy[:-2,1:-1])
        at[0,0]      =0.25*(dx[0,0])*(dy[0,0])
        at[0,1:-1]   =0.25*(dx[0,1:-1]+dx[0,:-2])*(dy[0,1:-1])
        at[1:-1,0]   =0.25*(dx[1:-1,0])*(dy[1:-1,0]+dy[:-2,0])
        at[0,-1]     =0.25*(dx[0,-1])*(dy[0,-1])
        at[-1,0]     =0.25*(dx[-1,0])*(dy[-1,0])
        at[-1,1:-1]  =0.25*(dx[-1,1:-1]+dx[-1,:-2])*(dy[-1,1:-1])
        at[1:-1,-1]  =0.25*(dx[1:-1,-1])*(dy[1:-1,-1]+dy[:-2,-1])
        at[-1,-1]    =0.25*(dx[-1,-1])*(dy[-1,-1])
        
        # Calculate heat generated at each node
#        sol=ct.Solution('')
#        states=ct.SolutionArray(sol, numpy.shape(dx))
#        Y_species={}
#        for i in range(len(self.species)):
#            Y_species[self.species[i]]=y_species[:,:,i]
#        states.TPY=T,101325,Y_species
        
        
        q_vol=0 # Volumetric heat generation rate
        
        q_source=q_vol*at
        
        return q_source