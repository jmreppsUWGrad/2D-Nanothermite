# -*- coding: utf-8 -*-
"""
Created on Tue Nov 06 15:19:15 2018

Notes on implementing Cantera:
    -sol=ct.Solution('___.cti') -> define solution mechanisms?
    -ct.SolutionArray() -> able to define an array of states
    https://cantera.org/documentation/docs-2.4/sphinx/html/cython/importing.html#representing-multiple-states
    

@author: Joseph
"""

import numpy as np
#import cantera as ct

class Source_terms():
    def __init__(self):
        self.R=8.314 # J/mol/K
        self.Ea=30000 # J/mol
        self.n=0.2 # Temperature exponent
#        self.species=key_species
    
    # Uniform volumetric generation
    def Source_Uniform(self, Q, dx, dy):
        at=np.zeros_like(dx)
        
        at[1:-1,1:-1]=0.25*(dx[1:-1,1:-1]+dx[1:-1,:-2])*(dy[1:-1,1:-1]+dy[:-2,1:-1])
        at[0,0]      =0.25*(dx[0,0])*(dy[0,0])
        at[0,1:-1]   =0.25*(dx[0,1:-1]+dx[0,:-2])*(dy[0,1:-1])
        at[1:-1,0]   =0.25*(dx[1:-1,0])*(dy[1:-1,0]+dy[:-2,0])
        at[0,-1]     =0.25*(dx[0,-1])*(dy[0,-1])
        at[-1,0]     =0.25*(dx[-1,0])*(dy[-1,0])
        at[-1,1:-1]  =0.25*(dx[-1,1:-1]+dx[-1,:-2])*(dy[-1,1:-1])
        at[1:-1,-1]  =0.25*(dx[1:-1,-1])*(dy[1:-1,-1]+dy[:-2,-1])
        at[-1,-1]    =0.25*(dx[-1,-1])*(dy[-1,-1])
        
        return Q*at
    
    # Calculate source term for combustion based on
    # K. Kim, "Computational Modeling of Combustion Wave in Nanoscale Thermite Reaction",
    # Int. J of Energy and Power engineering, vol.8, no.7, pp. 612-615, 2014.
    def Source_Comb_Kim(self, rho, T, eta, dx, dy, dt):
        at=np.zeros_like(dx)
        Ea=40000 # [J/mol] Approx value from Kim's paper
        A0=1e8 # [1/s] Fudged value
#        dH=1200 # [J/mol] Value taken from V. Baijot et al., A multi-phase ..., Combustion and Flame, 2017.
        dH=300000 #[J/kg] approx from ...
        
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
        
        detadt=A0*(1-eta)*np.exp(-Ea/self.R/T)
        eta+=dt*detadt
        
        return rho*at*dH*detadt
    
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