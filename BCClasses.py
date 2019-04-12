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
    def __init__(self, BC_dict, dx, dy):
        self.BCs=BC_dict
        self.dx,self.dy=dx,dy
        
    # Energy BCs
    def Energy(self, E, T_prev, dt, rho, Cv, vol, Ax, Ay):
        # Left face
        for i in range(len(self.BCs['bc_left_E'])/3):
            st=self.BCs['bc_left_E'][2+3*i][0]
            en=self.BCs['bc_left_E'][2+3*i][1]
            if self.BCs['bc_left_E'][3*i]=='T':
                E[st:en,0]=self.BCs['bc_left_E'][1+3*i]*rho[st:en,0]*Cv[st:en,0]*vol[st:en,0]
                if len(self.BCs['bc_left_E'])/3-i==1:
                    E[-1,0]=self.BCs['bc_left_E'][-2]*rho[-1,0]*Cv[-1,0]*vol[-1,0]
            
            else:
                if self.BCs['bc_left_E'][3*i]=='F':
                    q=self.BCs['bc_left_E'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_left_E'][1+3*i][0]*self.BCs['bc_left_E'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_left_E'][1+3*i][0]*T_prev[st:en,0] # h*Tij
                
                E[st:en,0]+=(Bi+q)*dt*Ax[st:en,0]
                if len(self.BCs['bc_left_E'])/3-i==1:
                    if self.BCs['bc_left_E'][3*i]=='C':
                        Bi=-self.BCs['bc_left_E'][1+3*i][0]*T_prev[-1,0] # h*Tij
                    E[-1,0]+=(Bi+q)*dt*Ax[-1,0]
        
        # Right face
        for i in range(len(self.BCs['bc_right_E'])/3):
            st=self.BCs['bc_right_E'][2+3*i][0]
            en=self.BCs['bc_right_E'][2+3*i][1]
            if self.BCs['bc_right_E'][3*i]=='T':
                E[st:en,-1]=self.BCs['bc_right_E'][1+3*i]*rho[st:en,-1]*Cv[st:en,-1]*vol[st:en,-1]
                if len(self.BCs['bc_right_E'])/3-i==1:
                    E[-1,-1]=self.BCs['bc_right_E'][-2]*rho[-1,-1]*Cv[-1,-1]*vol[-1,-1]
            
            else:
                if self.BCs['bc_right_E'][3*i]=='F':
                    q=self.BCs['bc_right_E'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_right_E'][1+3*i][0]*self.BCs['bc_right_E'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_right_E'][1+3*i][0]*T_prev[st:en,-1] # h*Tij
                
                E[st:en,-1]+=(Bi+q)*dt*Ax[st:en,-1]
                if len(self.BCs['bc_right_E'])/3-i==1:
                    if self.BCs['bc_right_E'][3*i]=='C':
                        Bi=-self.BCs['bc_right_E'][1+3*i][0]*T_prev[-1,-1] # h*Tij
                    E[-1,-1]+=(Bi+q)*dt*Ax[-1,-1]
        
        # South face
        for i in range(len(self.BCs['bc_south_E'])/3):
            st=self.BCs['bc_south_E'][2+3*i][0]
            en=self.BCs['bc_south_E'][2+3*i][1]
            if self.BCs['bc_south_E'][3*i]=='T':
                E[0,st:en]=self.BCs['bc_south_E'][1+3*i]*rho[0,st:en]*Cv[0,st:en]*vol[0,st:en]
                if len(self.BCs['bc_south_E'])/3-i==1:
                    E[0,-1]=self.BCs['bc_south_E'][-2]*rho[0,-1]*Cv[0,-1]*vol[0,-1]
            
            else:
                if self.BCs['bc_south_E'][3*i]=='F':
                    q=self.BCs['bc_south_E'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_south_E'][1+3*i][0]*self.BCs['bc_south_E'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_south_E'][1+3*i][0]*T_prev[0,st:en] # h*Tij
                
                E[0,st:en]+=(Bi+q)*dt*Ay[0,st:en]
                if len(self.BCs['bc_south_E'])/3-i==1:
                    if self.BCs['bc_south_E'][3*i]=='C':
                        Bi=-self.BCs['bc_south_E'][1+3*i][0]*T_prev[0,-1] # h*Tij
                    E[0,-1]+=(Bi+q)*dt*Ay[0,-1]
                    
        # North face
        for i in range(len(self.BCs['bc_north_E'])/3):
            st=self.BCs['bc_north_E'][2+3*i][0]
            en=self.BCs['bc_north_E'][2+3*i][1]
            if self.BCs['bc_north_E'][3*i]=='T':
                E[-1,st:en]=self.BCs['bc_north_E'][1+3*i]*rho[-1,st:en]*Cv[-1,st:en]*vol[-1,st:en]
                if len(self.BCs['bc_north_E'])/3-i==1:
                    E[-1,-1]=self.BCs['bc_north_E'][-2]*rho[-1,-1]*Cv[-1,-1]*vol[-1,-1]
            
            else:
                if self.BCs['bc_north_E'][3*i]=='F':
                    q=self.BCs['bc_north_E'][1+3*i]
                    Bi=0
                    
                else:
                    q=self.BCs['bc_north_E'][1+3*i][0]*self.BCs['bc_north_E'][1+3*i][1] # h*Tinf
                    Bi=-self.BCs['bc_north_E'][1+3*i][0]*T_prev[-1,st:en] # h*Tij
                
                E[-1,st:en]+=(Bi+q)*dt*Ay[-1,st:en]
                if len(self.BCs['bc_north_E'])/3-i==1:
                    if self.BCs['bc_north_E'][3*i]=='C':
                        Bi=-self.BCs['bc_north_E'][1+3*i][0]*T_prev[-1,-1] # h*Tij
                    E[-1,-1]+=(Bi+q)*dt*Ay[-1,-1]
        
        # Apply radiation BCs
        if self.BCs['bc_left_rad']!='None':
            E[:,0]+=Ax[:,0]*dt*\
                self.BCs['bc_left_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_left_rad'][1]**4-T_prev[:,0]**4)
        if self.BCs['bc_right_rad']!='None':
            E[:,-1]+=Ax[:,-1]*dt*\
                self.BCs['bc_right_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_right_rad'][1]**4-T_prev[:,-1]**4)
        if self.BCs['bc_south_rad']!='None':
            E[0,:]+=Ay[0,:]*dt*\
                self.BCs['bc_south_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_south_rad'][1]**4-T_prev[0,:]**4)
        if self.BCs['bc_north_rad']!='None':
            E[-1,:]+=Ay[-1,:]*dt*\
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