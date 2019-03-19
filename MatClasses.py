# -*- coding: utf-8 -*-
"""
Created on Sat Sep 29 12:58:39 2018

@author: Joseph

Class for calculating temperature dependent properties

"""

import numpy as np
#from CoolProp.CoolProp import PropsSI

# Class for calculating diffusion coefficients
class Diff_Coef():
    def __init__(self):
        self.D0_O2_Cu=1.16*10**(-6)
        self.D0_O2_Al2O3=9.0*10**(-5)
        self.D0_Al_Al2O3=1.0*10**(-8)
        
    def O2_Cu(self, T):
        # Baijot 2017, Brotman 2019, 
        # This also used for diffusion of O2 through CuO/Cu2O
        return self.D0_Cu*np.exp(-67300/8.314/T)
    
    def O2_Al2O3(self, T):
        # Brotman 2019, 
        return self.D0_Al2O3*np.exp(-140000/8.314/T)
    
    def Al_Al2O3(self,T):
        # Baijot 2017; experimental curve fitting
        return self.D0_Al_Al2O3*np.exp(-110000/8.314/T)
    
    # Main function to calculate diffusion coefficients
    # Dependent on specie, points to correct function
    def get_Diff(self,T,species):
        if species=='Al':
            return self.Al_Al2O3(T)
        elif species=='Cu':
            return self.O2_Cu(T)

# Class for calculating phase change rates
class Phase_Change():
    def __init__(self):
        # Enthalpies of vaporization
        self.dH_Al=294000 # kJ/mol; Baijot 2016
        self.dH_Cu=338000 # kJ/mol; Baijot 2016
        self.dH_Al2O3=1402000 # kJ/mol; Baijot 2016
        # Coefficients for calculating thermal conductivity
        self.k0
        self.k1
        self.k2
        self.k3
        # Coefficients for calcuating specific heat
        self.C0
        self.C1
        self.C2
        self.C3
        # Coefficients for calculating ???
        
    def Al(self, T):
        T_melt=np.where(T>=933,T,-1)
        T_melt=np.where(T_melt<2792,T_melt,-1)
        T_vap=np.where(T_melt>=2792, T_melt,-1)
        # Melting
        
        
        # Vaporization
        
        
        
        return rho
    
    def Al2O3(self, T):
        C # Calculate specific heat
        return C
    
    def CuO(self, T):
        k # Calcuate thermal conductivity here
        return k
    
    def Cu(self, T):
        return 2
    
    # Main function to calculate evaporation rates
    # Points to correct function based on specie
    def get_evap(self,T,species):
        if species=='Al':
            return self.Al(T)
        elif species=='Cu':
            return self.Cu(T)
    
class ArProp():
    def __init__(self):
        # Ceofficients for calculating density
        self.rho0
        self.rho1
        self.rho2
        self.rho3
        # Coefficients for calculating thermal conductivity
        self.k0
        self.k1
        self.k2
        self.k3
        # Coefficients for calcuating specific heat
        self.C0
        self.C1
        self.C2
        self.C3
        # Coefficients for calculating ???
        
    def CalcDens(self, T):
        rho # Calcuate density here
        return rho
    def CalcCond(self,T):
        k # Calcuate thermal conductivity here
        return k
    def CalcSpecHeat(self,T):
        C # Calculate specific heat
        return C