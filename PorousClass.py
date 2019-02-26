# -*- coding: utf-8 -*-
"""
######################################################
#             2D Heat Conduction Solver              #
#              Created by J. Mark Epps               #
#          Part of Masters Thesis at UW 2018-2020    #
######################################################

This file contains the Porous medium equations:
    -
    
Features:
    -

Desired:
    -Solve Darcy's law
    -return darcy velocities
    
"""

import numpy as np
#import string as st

class Porous_Equations():
    def __init__(self, K, mu):
        self.perm=K # Permeability
        self.mu=mu # Viscosity
        
    # Porous medium equations
    def Porous(self, P, rho, dx, dy):
        u_x=np.zeros_like(rho)
        u_y=np.zeros_like(rho)
        
        u_x[:,1:-1]=self.perm*(P[:,2:]-P[:,:-2])/(dx[:,1:-1]+dx[:,:-2])
        u_y[1:-1,:]=self.perm*(P[2:,:]-P[:-2,:])/(dy[1:-1,:]+dy[:-2,:])
    