# -*- coding: utf-8 -*-
"""
Created on Tue Nov 06 15:19:15 2018

Notes on implementing Cantera:
    -sol=ct.Solution('___.cti') -> define solution mechanisms?
    -ct.SolutionArray() -> able to define an array of states
    https://cantera.org/documentation/docs-2.4/sphinx/html/cython/importing.html#representing-multiple-states
    

@author: Joseph
"""

import numpy
#import cantera as ct

class Source_terms():
    def __init__(self, key_species):
        self.Ea=30000 # J/mol
        self.n=0.2 # Temperature exponent
        self.species=key_species
    
    # Uniform volumetric generation
    def Source_Uniform(self, Q, dx, dy):
        at=numpy.zeros_like(dx)
        
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
    
    # Calculate source term for combustion (NEEDS MODIFYING)
    def Source_Comb(self, T, y_species, dx, dy):
        # Terms to include after modelling combustion:
            #resulting expression should have dimensions of W/m
            #Be sure to account for area of CV properly
        at=numpy.zeros_like(dx)
        q_source=numpy.zeros_like(dx)
        
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