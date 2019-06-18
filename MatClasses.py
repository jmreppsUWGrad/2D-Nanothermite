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
        return self.D0_O2_Cu*np.exp(-67300/8.314/T)
    
    def O2_Al2O3(self, T):
        # Brotman 2019, 
        return self.D0_O2_Al2O3*np.exp(-140000/8.314/T)
    
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
        elif species=='g':
            return self.O2_Al2O3(T)
        else:
            return 0

# Class for calculating specific heat at constant pressure
# Returned units are in J/kg/K
class Cp():
    def __init__(self):
        # Molar masses
        self.Al_mol_mass=26.982
        self.Al2O3_mol_mass=101.96
        self.CuO_mol_mass=79.546
        self.Cu_mol_mass=63.546
        self.Ar_mol_mass=39.948
        # Density dictionary for gases
        self.rho={'Air': 1.2, 'Ar': 1.8}
        
    def Al(self, T, typ):
        molar_mass=26.982
        Cp=np.zeros_like(T)
        
        # Solid phase
        # Polynomial fit coefficients
        # Excel regression of JANAF data
        a0=11.674
        a1=0.075938
        a2=-1.5457e-4
        a3=1.5259e-7
        a4=-5.1417e-11
        
        Cp[T<933]=(a0+a1*T[T<933]+a2*T[T<933]**2+a3*T[T<933]**3+a4*T[T<933]**4)*1000/molar_mass
        
        # Gas phase; use value near 2900 K
        if typ=='Cp':
            Cp[T>2791]=771.0
        else:
            Cp[T>2791]=771.0-8.314*1000/self.Al_mol_mass
        
        # Liquid phase
        Cp[Cp==0]=1177.0
        
        return Cp
    
    def Al2O3(self, T, typ):
        molar_mass=101.96
        Cp=np.zeros_like(T)
        
        # Solid phase
        # Coefficicents for polynomial fit
        # Excel regression of JANAF data
        a0=17.34
        a1=2.8436e-1
        a2=-2.7874e-4
        a3=1.2282e-7
        a4=-1.9758e-11
        
        Cp[T<2327]=(a0+a1*T[T<2327]+a2*T[T<2327]**2+a3*T[T<2327]**3+a4*T[T<2327]**4)*1000/molar_mass
        
        # Liquid phase
        Cp[Cp==0]=1888.0
        
        return Cp
    
    def CuO(self, T, typ):
        molar_mass=79.546
        # Solid phase only
#        # Coefficicents for polynomial fit (CuO base; up to 2000 K)
#        a0=24.56
#        a1=8.3286e-2
#        a2=-9.0277e-5
#        a3=4.6464e-8
#        a4=-8.6691e-12
        
        # Coefficicents for polynomial fit (CuO (g); up to 6000 K)
        # Excel regression of JANAF data
        a0=30.821
        a1=2.6747e-2
        a2=-4.4523e-5
        a3=3.5615e-8
        a4=-1.0947e-11
        
        return (a0+a1*T+a2*T**2+a3*T**3+a4*T**4)*1000/molar_mass
        
    def Cu(self, T, typ):
        molar_mass=63.546
        Cp=np.zeros_like(T)
        # Solid phase
        # Coefficicents for polynomial fit
        # Excel regression of JANAF data
        a0=11.674
        a1=0.075938
        a2=-1.5457e-4
        a3=1.5259e-7
        a4=-5.1417e-11
        Cp[T<1358]=(a0+a1*T[T<1358]+a2*T[T<1358]**2+a3*T[T<1358]**3+a4*T[T<1358]**4)*1000/molar_mass
        
        # Gas phase; use value near 2900 K
        if typ=='Cp':
            Cp[T>2843]=388.0
        else:
            Cp[T>2843]=388.0-8.314*1000/self.Cu_mol_mass
        
        # Liquid phase
        Cp[Cp==0]=517.0
        
        return Cp
    
    def Ar(self, T, typ):
        if typ=='Cp':
            return np.ones_like(T)*520.0
        else:
            return np.ones_like(T)*520.0-8.314*1000/self.Ar_mol_mass
        
    def Air(self, T, typ):
        molar_mass=28.97
        Cp=np.zeros_like(T)
        
#        # Coefficicents for polynomial fit (273-1800 K)
#        # Taken from Cengel and Boles, Thermodynamics: An engineering approach
#        a0=28.11
#        a1=0.1967e-2
#        a2=0.4802e-5
#        a3=-1.966e-9
#        
#        Cp=(a0+a1*T+a2*T**2+a3*T**3)/molar_mass*1000
        
        # Coefficicents for polynomial fit (300-3000 K)
        # Excel regression of data taken from Bergman et al., Fundamentals of Heat and mass transfer
        a0=1.0718e3
        a1=-4.664e-1
        a2=1.0584e-3
        a3=-6.6196e-7
        a4=1.4070e-10
        
        Cp=a0+a1*T+a2*T**2+a3*T**3+a4*T**4
        
        if typ=='Cp':
            return Cp
        else:
            return Cp-8.314*1000/molar_mass
    
    # Main function to calculate specific heat at constant pressure
    # Points to correct function based on specie
    def get_Cp(self,T,species):
        if species=='Al':
            return self.Al(T,'Cp')
        elif species=='Cu':
            return self.Cu(T,'Cp')
        elif species=='Al2O3':
            return self.Al2O3(T,'Cp')
        elif species=='CuO':
            return self.CuO(T,'Cp')
        elif species=='Air':
            return self.Air(T,'Cp')
        else:
            return self.Ar(T,'Cp')
        
    # Main function to calculate specific heat at constant volume
    # Points to correct function based on specie
    def get_Cv(self,T,species):
        if species=='Al':
            return self.Al(T,'Cv')
        elif species=='Cu':
            return self.Cu(T,'Cv')
        elif species=='Al2O3':
            return self.Al2O3(T,'Cv')
        elif species=='CuO':
            return self.CuO(T,'Cv')
        elif species=='Air':
            return self.Air(T,'Cv')
        else:
            return self.Ar(T,'Cv')

# Return thermal conductivity
class therm_cond():
    def __init__(self):
        self.num='dummy'
        
    def Ar(self, T):
        k=np.zeros_like(T)
        
        # Coefficicents for quadratic fit (338-2518 K)
        # Excel regression of data taken from S.H.P Chen and S.C. Saxena, "Thermal conductivity of argon in temperature range 350 to 2500 K"
        a0=8.1079
        a1=3.9379e-2
        a2=-4.6128e-6
                
        k=(a0+a1*T+a2*T**2)/1000
        
        return k
        
    def Air(self, T):
        k=np.zeros_like(T)
        
        # Coefficicents for polynomial fit (300-3000 K)
        # Excel regression of data taken from Bergman et al., Fundamentals of Heat and mass transfer
        a0=42.467
        a1=-0.0942
        a2=2.3531e-4
        a3=-1.451e-7
        a4=3.1151e-11
        
        k=(a0+a1*T+a2*T**2+a3*T**3+a4*T**4)/1000
        
        return k
    
    # Main function to calculate specific heat at constant pressure
    # Points to correct function based on specie
    def get_k(self,T,species):
        if species=='Air':
            return self.Air(T)
        else:
            return self.Ar(T)