from .kepler import keplerian_system
import numpy as np
from scipy.optimize import minimize
from astropy import constants as ac
from .phot import *
from copy import deepcopy
import ctypes
import matplotlib.pyplot as plt
from matplotlib import animation

photlib = ctypes.CDLL("./c_src/phot.so")
keplib = ctypes.CDLL("./c_src/kepler.so")

__all__ = ['system']

class system:
    
    days_in_year = 365.256
    earths_in_sun = 332946.08
    G = 39.478 # au ^ 3 / (yr ^2 * M_sun)
    G = G / (days_in_year ** 2) / earths_in_sun
    
    def __init__(self, star, planet, moon):
        self.star = deepcopy(star)
        self.star.mass = star.mass * ac.M_sun.value / ac.M_earth.value
        self.star.radius = star.radius * ac.R_sun.value / ac.R_earth.value
        self.planet = deepcopy(planet)
        self.moon = deepcopy(moon)
        
        self.starplanet = keplerian_system(self.star, self.planet)
        self.planetmoon = keplerian_system(self.planet, self.moon)
        
    def set_planet_orbit(self, t0=0, e=0, P=days_in_year, Omega=0, w=0, i=90):
        self.starplanet.set_orbital_parameters(t0, e, P, Omega, w, i)
        
    def set_moon_orbit(self, t0=0, e=0, P=days_in_year, Omega=0, w=0, i=90):
        self.planetmoon.set_orbital_parameters(t0, e, P, Omega, w, i)
        
    def coords_dep(self, t):
        
        tp = t + self.starplanet.find_transit() - self.starplanet.t0
        
        starplanet_center = self.starplanet.barycentric_coords(tp)
        planetmoon_center = self.planetmoon.barycentric_coords(tp)
        starcoords = starplanet_center[0]
        planetcoords = starplanet_center[1] + planetmoon_center[0]
        mooncoords = starplanet_center[1] + planetmoon_center[1]
        
        return starcoords, planetcoords, mooncoords
    
    def sky_projected_distance_dep(self, t):
        
        st, pl, mo = self.coords(t)
        
        stx, sty, stz = st
        plx, ply, plz = pl
        mox, moy, moz = mo
        
        starplanet_distance = np.sqrt((stx - plx)**2 + (sty + ply)**2)
        starmoon_distance = np.sqrt((stx - mox)**2 + (sty + moy)**2)
        planetmoon_distance = np.sqrt((plx - mox)**2 + (ply - moy)**2)
        
        return (starplanet_distance, starmoon_distance, planetmoon_distance)
    
    def starplanet_distance_dep(self, t):
        
        tp = t + self.starplanet.find_transit() - self.starplanet.t0
        f = self.starplanet.f(tp)
        
        return (self.starplanet.a * (1 - self.starplanet.e**2) / 
                (1 + self.starplanet.e * np.cos(f)) * 
                np.sqrt(1 - (np.sin(self.starplanet.w + f) * np.sin(self.starplanet.i))**2))
    
    def find_contacts_dep(self):
        
        mr = self.moon.radius * ac.R_earth.value / ac.au.value
        sr = self.star.radius * ac.R_earth.value / ac.au.value
        t_trans = self.starplanet.t0
        d = self.planetmoon.a + mr + sr
   
        g = lambda t: np.abs(self.starplanet_distance(t) - d)
        res = minimize(g, x0=self.starplanet.t0, method='Nelder-Mead').x[0]
        return (t_trans - np.abs(res - t_trans), t_trans + np.abs(res - t_trans))
    
    def find_dt_dep(self):
        
        f_trans = np.pi/2 - self.starplanet.w
        r_trans = (self.starplanet.a * (1 - self.starplanet.e ** 2) / 
                   (1 + self.starplanet.e * np.cos(f_trans)))
        
        strad_au = self.star.radius * ac.R_earth.value / ac.au.value
        df = np.arcsin((strad_au + (self.planetmoon.e + 1) * self.planetmoon.a) / r_trans)
        
        fa = f_trans - df
        
        e = self.starplanet.e
        dfa = (2 * np.sqrt(1 - e ** 2) * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(fa / 2)) - 
               (e * (1 - e**2) * np.sin(fa)) / (1 + e * np.cos(fa)))
        dfb = (2 * np.sqrt(1 - e ** 2) * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(f_trans / 2)) - 
               (e * (1 - e**2) * np.sin(f_trans)) / (1 + e * np.cos(f_trans)))
        
        h = np.sqrt(self.G * (self.star.mass + self.planet.mass + self.moon.mass) * 
                    self.starplanet.a * (1 - e**2))
        dt = (dfb - dfa) * (self.starplanet.a ** 2) / h
        
        return dt
            
    def reduce_t_dep(self, t):
        
        t_trans = self.starplanet.t0
        dt = self.find_dt()
        ta, tb = t_trans - dt, t_trans + dt
        nP = np.int((t[-1] - t[0]) / P)
        centers = [t_trans + i*P for i in range(nP)]
        is_transiting = (np.isclose(t % self.starplanet.P - t_trans, np.zeros_like(t), atol = contact_duration)
                         | np.isclose(t % self.starplanet.P - t_trans, np.ones_like(t)*self.starplanet.P, atol = contact_duration))

        return is_transiting
        
    def flux_dep(self, t):
        return lc(self, t, self.star.u)
    
    def flux(self, t, ld='quad'):
        
        ld_coeffs = self.star.u
        if ld == 'quad':
            l1, l2 = ld_coeffs
            c = np.array([0, l1 + 2*l2, 0, -l2])
            I_func = I_nonlinear
        elif ld == 'nonlinear':
            c = ld_coeffs
            I_func = I_nonlinear
        else:
            raise AttributeError('ld must be one of quad or nonlinear')
                                    
        pp = self.planet.radius / self.star.radius
        pm = self.moon.radius / self.star.radius
        strad_au = self.star.radius * ac.R_earth.value / ac.au.value
        c1, c2, c3, c4 = c
        
        f = (ctypes.c_double * len(t))(*np.zeros(len(t)))
        photlib.flux(f, 
                     (ctypes.c_double * len(t))(*t), 
                     ctypes.c_double(pp), 
                     ctypes.c_double(pm), 
                     ctypes.c_double(c1), 
                     ctypes.c_double(c2), 
                     ctypes.c_double(c3), 
                     ctypes.c_double(c4),
                     ctypes.c_double(self.starplanet.t0),
                     ctypes.c_double(self.planetmoon.t0),
                     ctypes.c_double(self.starplanet.P),
                     ctypes.c_double(self.starplanet.n),
                     ctypes.c_double(self.starplanet.e),
                     ctypes.c_double(self.starplanet.a),
                     ctypes.c_double(self.planet.mass),
                     ctypes.c_double(self.star.mass),
                     ctypes.c_double(self.moon.mass),
                     ctypes.c_double(self.starplanet.w),
                     ctypes.c_double(self.starplanet.Omega),
                     ctypes.c_double(self.starplanet.i),
                     ctypes.c_double(self.planetmoon.n),
                     ctypes.c_double(self.planetmoon.e),
                     ctypes.c_double(self.planetmoon.a),
                     ctypes.c_double(self.planetmoon.w),
                     ctypes.c_double(self.planetmoon.Omega),
                     ctypes.c_double(self.planetmoon.i),
                     ctypes.c_double(strad_au),
                     ctypes.c_int(len(f)))
        
        return np.array(f)
    
    def animate(self, t, stkwargs={'color':'k', 'fill':False}, 
                    plkwargs={'color':'k', 'fill':False}, 
                    mokwargs={'color':'k', 'fill':False}):
        
        fig = plt.figure(figsize=(10, 10))
        ax = plt.gca()
        plt.xlim(-1.2, 1.2)
        plt.ylim(-1.2, 1.2)

        st_patch = ax.add_patch(plt.Circle((0, 0), 0, **stkwargs))
        pl_patch = ax.add_patch(plt.Circle((0, 0), 0, **plkwargs))
        mo_patch = ax.add_patch(plt.Circle((0, 0), 0, **mokwargs))

        def init():

            pp = self.planet.radius / self.star.radius
            pm = self.moon.radius / self.star.radius
            st_patch.set_center((0, 0))
            st_patch.set_radius(1)
            pl_patch.set_center((0, 0))
            pl_patch.set_radius(pp)
            mo_patch.set_center((0, 0))
            mo_patch.set_radius(pm)
            return st_patch, pl_patch, mo_patch

        def update(t):
    
            pp = self.planet.radius / self.star.radius
            pm = self.moon.radius / self.star.radius
            strad_au = self.star.radius * ac.R_earth.value / ac.au.value
        
            st, pl, mo = self.coords(np.array([t]))
            stx, sty, stz = st / strad_au
            plx, ply, plz = pl / strad_au
            mox, moy, moz = mo / strad_au
            st_patch.set_center((stx, sty))
            pl_patch.set_center((plx, ply))
            mo_patch.set_center((mox, moy))
            return st_patch, pl_patch, mo_patch
    

        return animation.FuncAnimation(fig, update, frames=t,
                    init_func=init, blit=True, interval=20)
    
    def draw_config(self, ax, t, stkwargs={'color':'k', 'fill':False}, 
                    plkwargs={'color':'k', 'fill':False}, 
                    mokwargs={'color':'k', 'fill':False}):
        
        pp = self.planet.radius / self.star.radius
        pm = self.moon.radius / self.star.radius
        strad_au = self.star.radius * ac.R_earth.value / ac.au.value
        
        st, pl, mo = self.coords(np.array([t]))
        stx, sty, stz = st / strad_au
        plx, ply, plz = pl / strad_au
        mox, moy, moz = mo / strad_au
        
        starpatch = plt.Circle((stx, sty), 1, animated=True, **stkwargs)
        planetpatch = plt.Circle((plx, ply), pp, animated=True, **plkwargs)
        moonpatch = plt.Circle((mox, moy), pm, animated=True, **mokwargs)
        
        #fig = pl.figure(figsize=(10, 10));
        #ax = pl.gca()
        ax.add_patch(starpatch)
        ax.add_patch(planetpatch)
        ax.add_patch(moonpatch)
        
        return ax