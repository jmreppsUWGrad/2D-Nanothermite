# -*- coding: utf-8 -*-
"""
Created on Tue Nov 06 15:19:15 2018

@author: Joseph
"""

import numpy

class Source_terms():
    def __init__(self):
        self.Ea=30000 # J/mol
        self.n=0.2 # Temperature exponent
    
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