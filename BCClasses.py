# -*- coding: utf-8 -*-
"""
######################################################
#             2D Heat Conduction Solver              #
#              Created by J. Mark Epps               #
#          Part of Masters Thesis at UW 2018-2020    #
######################################################

This file contains functions to do boundary conditions:
    -
    
Features:
    -

Desired:
    -
    -
    
"""

import numpy as np
#import string as st

class BCs():
    def __init__(self, BC_dict, dx, dy, domain):
        self.BCs=BC_dict
        self.dx,self.dy=dx,dy
        self.X=np.empty_like(dx)
        self.domain=domain
        
    # Ablation considerations (basic)
    def flux_abl(self, T, E, dt, h, q):
        T_high=800.0
        E_c=E.copy()
        h_c=h.copy()
        T_c=T.copy()
        E_c[T_c>T_high]-=q*dt/h_c[T_c>T_high]*1.1
        
#        E[T>T_high]-=q*dt/h[T>T_high]*1.1
#        return E
    
    # Energy BCs
    def Energy(self, E, T_prev, dt, rhoC, hx, hy):
        # Left face
        for i in range(len(self.BCs['bc_left_E'])/3):
            st=self.BCs['bc_left_E'][2+3*i][0]
            en=self.BCs['bc_left_E'][2+3*i][1]
            if self.BCs['bc_left_E'][3*i]=='T':
                E[st:en,0]=self.BCs['bc_left_E'][1+3*i]*rhoC[st:en,0]
                
            elif self.domain!='Axisymmetric':
                if self.BCs['bc_left_E'][3*i]=='F':
                    q=self.BCs['bc_left_E'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_left_E'][1+3*i][0]*self.BCs['bc_left_E'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_left_E'][1+3*i][0]*T_prev[st:en,0] # h*Tij
                
                E[st:en,0]+=(Bi+q)*dt/hx[st:en,0]
                
                
        # Right face
        for i in range(len(self.BCs['bc_right_E'])/3):
            st=self.BCs['bc_right_E'][2+3*i][0]
            en=self.BCs['bc_right_E'][2+3*i][1]
            if self.BCs['bc_right_E'][3*i]=='T':
                E[st:en,-1]=self.BCs['bc_right_E'][1+3*i]*rhoC[st:en,-1]
                
            else:
                if self.BCs['bc_right_E'][3*i]=='F':
                    q=self.BCs['bc_right_E'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_right_E'][1+3*i][0]*self.BCs['bc_right_E'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_right_E'][1+3*i][0]*T_prev[st:en,-1] # h*Tij
                
                if self.domain=='Axisymmetric':
                    E[st:en,-1]+=(Bi+q)*dt*self.X[st:en,-1]\
                        /hx[st:en,-1]/(self.X[st:en,-1]-self.dx[st:en,-2])
                else:
                    E[st:en,-1]+=(Bi+q)*dt/hx[st:en,-1]
                
        # South face
        for i in range(len(self.BCs['bc_south_E'])/3):
            st=self.BCs['bc_south_E'][2+3*i][0]
            en=self.BCs['bc_south_E'][2+3*i][1]
            if self.BCs['bc_south_E'][3*i]=='T':
                E[0,st:en]=self.BCs['bc_south_E'][1+3*i]*rhoC[0,st:en]
            
            else:
                if self.BCs['bc_south_E'][3*i]=='F':
                    q=self.BCs['bc_south_E'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_south_E'][1+3*i][0]*self.BCs['bc_south_E'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_south_E'][1+3*i][0]*T_prev[0,st:en] # h*Tij
                
                E[0,st:en]+=(Bi+q)*dt/hy[0,st:en]
                
        # North face
        for i in range(len(self.BCs['bc_north_E'])/3):
            st=self.BCs['bc_north_E'][2+3*i][0]
            en=self.BCs['bc_north_E'][2+3*i][1]
            if self.BCs['bc_north_E'][3*i]=='T':
                E[-1,st:en]=self.BCs['bc_north_E'][1+3*i]*rhoC[-1,st:en]
            
            else:
                if self.BCs['bc_north_E'][3*i]=='F':
                    q=self.BCs['bc_north_E'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_north_E'][1+3*i][0]*self.BCs['bc_north_E'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_north_E'][1+3*i][0]*T_prev[-1,st:en] # h*Tij
                
                E[-1,st:en]+=(Bi+q)*dt/hy[-1,st:en]
#                E[T_prev>800]-=(Bi+q)*dt/hy[T_prev>800]*1.1
#                E[-1,st:en]=self.flux_abl(T_prev[-1,st:en], E[-1,st:en], dt, hy[-1,st:en], Bi+q)
#                self.flux_abl(T_prev[-1,st:en], E[-1,st:en], dt, hy[-1,st:en], Bi+q)
                
        # Apply radiation BCs
        if self.BCs['bc_left_rad']!='None' and self.domain!='Axisymmetric':
            E[:,0]+=dt/hx[:,0]*\
                self.BCs['bc_left_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_left_rad'][1]**4-T_prev[:,0]**4)
        if self.BCs['bc_right_rad']!='None':
            E[:,-1]+=dt/hx[:,-1]*\
                self.BCs['bc_right_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_right_rad'][1]**4-T_prev[:,-1]**4)
        if self.BCs['bc_south_rad']!='None':
            E[0,:]+=dt/hy[0,:]*\
                self.BCs['bc_south_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_south_rad'][1]**4-T_prev[0,:]**4)
        if self.BCs['bc_north_rad']!='None':
            E[-1,:]+=dt/hy[-1,:]*\
                self.BCs['bc_north_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_north_rad'][1]**4-T_prev[-1,:]**4)
    
    # Conservation of mass BCs
    def mass(self, m, P, Ax, Ay, vol):
        # Left face
        for i in range(len(self.BCs['bc_left_mass'])/3):
            st=self.BCs['bc_left_mass'][2+3*i][0]
            en=self.BCs['bc_left_mass'][2+3*i][1]
            # Gradient
            if self.BCs['bc_left_mass'][3*i]=='grad':
                m[st:en,0]=m[st:en,1]-self.BCs['bc_left_mass'][1+3*i]*self.dx[st:en,0]
                if len(self.BCs['bc_left_mass'])/3-i==1:
                    m[-1,0]=m[-1,1]-self.BCs['bc_left_mass'][-2]*self.dx[-1,0]
            # Pressure flux
            elif self.BCs['bc_left_mass'][3*i]=='grad_P':
                m[st:en,0]=m[st:en,1]-self.BCs['bc_left_mass'][1+3*i]*self.dx[st:en,0]
                if len(self.BCs['bc_left_mass'])/3-i==1:
                    m[-1,0]=m[-1,1]-self.BCs['bc_left_mass'][-2]*self.dx[-1,0]
            # Constant
            else:
                m[st:en,0]=self.BCs['bc_left_mass'][1+3*i]*vol[st:en,0]
                if len(self.BCs['bc_left_mass'])/3-i==1:
                    m[-1,0]=self.BCs['bc_left_mass'][-2]*vol[-1,0]
        # Right face
        for i in range(len(self.BCs['bc_right_mass'])/3):
            st=self.BCs['bc_right_mass'][2+3*i][0]
            en=self.BCs['bc_right_mass'][2+3*i][1]
            if self.BCs['bc_right_mass'][3*i]=='grad':
                m[st:en,-1]=self.BCs['bc_right_mass'][1+3*i]*self.dx[st:en,-1]+m[st:en,-2]
                if len(self.BCs['bc_right_mass'])/3-i==1:
                    m[-1,-1]=self.BCs['bc_right_mass'][-2]*self.dx[-1,-1]+m[-1,-2]
        
        # South face
        for i in range(len(self.BCs['bc_south_mass'])/3):
            st=self.BCs['bc_south_mass'][2+3*i][0]
            en=self.BCs['bc_south_mass'][2+3*i][1]
            if self.BCs['bc_south_mass'][3*i]=='grad':
                m[0,st:en]=m[1,st:en]-self.BCs['bc_south_mass'][1+3*i]*self.dy[0,st:en]
                if len(self.BCs['bc_south_mass'])/3-i==1:
                    m[0,-1]=m[1,-1]-self.BCs['bc_south_mass'][-2]*self.dy[0,-1]
                    
        # North face
        for i in range(len(self.BCs['bc_north_mass'])/3):
            st=self.BCs['bc_north_mass'][2+3*i][0]
            en=self.BCs['bc_north_mass'][2+3*i][1]
            if self.BCs['bc_north_mass'][3*i]=='grad':
                m[-1,st:en]=self.BCs['bc_north_mass'][1+3*i]*self.dy[-1,st:en]+m[-2,st:en]
                if len(self.BCs['bc_north_mass'])/3-i==1:
                    m[-1,-1]=self.BCs['bc_north_mass'][-2]*self.dy[-1,-1]+m[-2,-1]
        return 0
    
    # Pressure BCs (eventually lead to momentum)
    def P(self, P):
        # Left face
        for i in range(len(self.BCs['bc_left_P'])/3):
            st=self.BCs['bc_left_P'][2+3*i][0]
            en=self.BCs['bc_left_P'][2+3*i][1]
            if self.BCs['bc_left_P'][3*i]=='grad':
                P[st:en,0]=P[st:en,1]-self.BCs['bc_left_P'][1+3*i]*self.dx[st:en,0]
                if len(self.BCs['bc_left_P'])/3-i==1:
                    P[-1,0]=P[-1,1]-self.BCs['bc_left_P'][-2]*self.dx[-1,0]
            
        # Right face
        for i in range(len(self.BCs['bc_right_P'])/3):
            st=self.BCs['bc_right_P'][2+3*i][0]
            en=self.BCs['bc_right_P'][2+3*i][1]
            if self.BCs['bc_right_P'][3*i]=='grad':
                P[st:en,-1]=self.BCs['bc_right_P'][1+3*i]*self.dx[st:en,-1]+P[st:en,-2]
                if len(self.BCs['bc_right_P'])/3-i==1:
                    P[-1,-1]=self.BCs['bc_right_P'][-2]*self.dx[-1,-1]+P[-1,-2]
        
        # South face
        for i in range(len(self.BCs['bc_south_P'])/3):
            st=self.BCs['bc_south_P'][2+3*i][0]
            en=self.BCs['bc_south_P'][2+3*i][1]
            if self.BCs['bc_south_P'][3*i]=='grad':
                P[0,st:en]=P[1,st:en]-self.BCs['bc_south_P'][1+3*i]*self.dy[0,st:en]
                if len(self.BCs['bc_south_P'])/3-i==1:
                    P[0,-1]=P[1,-1]-self.BCs['bc_south_P'][-2]*self.dy[0,-1]
                    
        # North face
        for i in range(len(self.BCs['bc_north_P'])/3):
            st=self.BCs['bc_north_P'][2+3*i][0]
            en=self.BCs['bc_north_P'][2+3*i][1]
            if self.BCs['bc_north_P'][3*i]=='grad':
                P[-1,st:en]=self.BCs['bc_north_P'][1+3*i]*self.dy[-1,st:en]+P[-2,st:en]
                if len(self.BCs['bc_north_P'])/3-i==1:
                    P[-1,-1]=self.BCs['bc_north_P'][-2]*self.dy[-1,-1]+P[-2,-1]
            
        return 0