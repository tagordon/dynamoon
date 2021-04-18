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
        
    def relative_coords(self, t):
                        
        if self.e == 0:
            f = (t % self.P) * 2 * np.pi / self.P
            r = self.a
        else:
            f, r = self.solve_kepler(t)            

        x = -r * np.cos(self.w + f)
        y = -r * np.sin(self.w + f)*np.cos(self.i)
        z = r * np.sin(self.w + f)*np.sin(self.i)
        
        return np.array([x, y, z])
    
    def barycentric_coords(self, t):
        if not self.hasparams:
            raise AttributeError("must set orbital parameters prior to computing coordinates")
            
        coords = self.relative_coords(t)
        mr_primary = self.secondary.mass / (self.primary.mass + self.secondary.mass)
        mr_secondary = self.primary.mass / (self.primary.mass + self.secondary.mass)
        return np.array([-mr_primary * coords, mr_secondary * coords])  
        
    def find_transit(self):
        
        if self.e == 0:
            return (np.pi / 2 - self.w) / (2 * np.pi) * self.P
        else:
            arg = (1 - (1 - self.e**2) / (1 + self.e*np.cos(np.pi/2 - self.w)))/self.e
            if arg > 1:
                arg = 1
            E = np.arccos(arg)
        M = E - self.e*np.sin(E)
        return M / self.n + self.t0
    
    def solve_kepler(self, t):
        
        r = (ctypes.c_double * len(t))(*np.zeros(len(t)))
        f = (ctypes.c_double * len(t))(*np.zeros(len(t)))
        keplib.solve_kepler_array(r, f,
                                (ctypes.c_double * len(t))(*t),
                                ctypes.c_double(self.n),
                                ctypes.c_double(self.t0),
                                ctypes.c_double(self.e),
                                ctypes.c_double(self.a),
                                ctypes.c_int(len(t)))  
        return np.array(f), np.array(r)
        
        #M = self.M(t)
        #if self.e == 0:
        #    f = 2*np.arctan(np.sqrt((1 + self.e) / (1 - self.e)) * np.tan(M / 2))
        #    r = self.a * (1 - self.e * np.cos(M))
        #    return f, r
        #else:
        #    E = self.E(t)
        #    f = 2*np.arctan(np.sqrt((1 + self.e) / (1 - self.e)) * np.tan( E / 2))
        #    r = self.a * (1 - self.e * np.cos(E))
        #    return f, r
    
    def f(self, t):
        if self.e == 0:
            return ((t - self.t0)/self.P) % self.P * 2 * np.pi
        return 2*np.arctan(np.sqrt((1 + self.e) / (1 - self.e)) * np.tan(self.E(t) / 2))
    
    def r(self, t):
        if self.e == 0:
            return self.a
        return self.a * (1 - self.e * np.cos(self.E(t)))
    
    def E(self, t, method='newton'):
        
        M = self.n * (t - self.t0)
        if self.e == 0:
            return M
        g = lambda E, t: E - self.e * np.sin(E) - self.M(t)
        dg = lambda E, t: 1 - self.e * np.cos(E)
        ddg = lambda E, t: 1 + self.e * np.sin(E)
        
        if type(t) is not np.ndarray:
            res = root_scalar(g, x0 = self.M(t), method=method, fprime=dg, args=(t)).root
        else:
            res = np.array([root_scalar(g, x0 = self.M(x), 
                                    method=method, 
                                    fprime=dg,
                                    fprime2=ddg,
                                    xtol = 1e-7,
                                    rtol = 1e-2,
                                    args=(x)).root for x in t])
        return res
    
    def M(self, t):
        return self.n * (t - self.t0)