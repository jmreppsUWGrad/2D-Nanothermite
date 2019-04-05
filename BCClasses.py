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
    def Energy(self, E, T_prev, dt, rho, Cv, vol):
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
                
                E[st:en,0]+=(Bi+q)*self.dy[st:en,0]*dt
                if len(self.BCs['bc_left_E'])/3-i==1:
                    if self.BCs['bc_left_E'][3*i]=='C':
                        Bi=-self.BCs['bc_left_E'][1+3*i][0]*T_prev[-1,0] # h*Tij
                    E[-1,0]+=(Bi+q)*self.dy[-1,0]*dt
        
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
                
                E[st:en,-1]+=(Bi+q)*self.dy[st:en,-1]*dt
                if len(self.BCs['bc_right_E'])/3-i==1:
                    if self.BCs['bc_right_E'][3*i]=='C':
                        Bi=-self.BCs['bc_right_E'][1+3*i][0]*T_prev[-1,-1] # h*Tij
                    E[-1,-1]+=(Bi+q)*self.dy[-1,-1]*dt
        
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
                
                E[0,st:en]+=(Bi+q)*self.dx[0,st:en]*dt
                if len(self.BCs['bc_south_E'])/3-i==1:
                    if self.BCs['bc_south_E'][3*i]=='C':
                        Bi=-self.BCs['bc_south_E'][1+3*i][0]*T_prev[0,-1] # h*Tij
                    E[0,-1]+=(Bi+q)*self.dx[0,-1]*dt
                    
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
                
                E[-1,st:en]+=(Bi+q)*self.dx[-1,st:en]*dt
                if len(self.BCs['bc_north_E'])/3-i==1:
                    if self.BCs['bc_north_E'][3*i]=='C':
                        Bi=-self.BCs['bc_north_E'][1+3*i][0]*T_prev[-1,-1] # h*Tij
                    E[-1,-1]+=(Bi+q)*self.dx[-1,-1]*dt
        
        # Apply radiation BCs
        if self.BCs['bc_left_rad']!='None':
            E[:,0]+=self.dy[:,0]*dt*\
                self.BCs['bc_left_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_left_rad'][1]**4-T_prev[:,0]**4)
        if self.BCs['bc_right_rad']!='None':
            E[:,-1]+=self.dy[:,-1]*dt*\
                self.BCs['bc_right_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_right_rad'][1]**4-T_prev[:,-1]**4)
        if self.BCs['bc_south_rad']!='None':
            E[0,:]+=self.dx[0,:]*dt*\
                self.BCs['bc_south_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_south_rad'][1]**4-T_prev[0,:]**4)
        if self.BCs['bc_north_rad']!='None':
            E[-1,:]+=self.dx[-1,:]*dt*\
                self.BCs['bc_north_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_north_rad'][1]**4-T_prev[-1,:]**4)
    
    # Conservation of mass BCs
    def mass(self, m, Ax,Ay,vol):
        return 0