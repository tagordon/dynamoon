import numpy as np
from scipy.integrate import quad
from scipy.optimize import root_scalar 
from scipy.optimize import minimize
import ctypes

keplib = ctypes.CDLL("./c_src/kepler.so")

__all__ = ['keplerian_system']

class keplerian_system:
    
    days_in_year = 365.256
    earths_in_sun = 332946.08
    G = 39.478 # au ^ 3 / (yr ^2 * M_sun)
    G = G / (days_in_year ** 2) / earths_in_sun
    
    def __init__(self, primary, secondary):
        self.primary = primary
        self.secondary = secondary
        pass
    
    def set_orbital_parameters(self, t0=0, e=0, P=days_in_year, Omega=0, w=0, i=90):
        self.t0 = t0
        self.e = e
        self.P = P
        self.Omega = Omega * np.pi / 180
        self.w = w * np.pi / 180
        self.i = i * np.pi / 180
        self.n = 2 * np.pi / self.P
        self.a = (self.G * (self.primary.mass + self.secondary.mass) 
                  / (self.n ** 2)) ** (1 / 3)
        
        self.hasparams = True