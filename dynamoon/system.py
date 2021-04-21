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
        
    def coords(self, t):
        
        keplib.find_transit.restype = ctypes.c_double
        tt = keplib.find_transit(ctypes.c_double(self.starplanet.e), 
                                 ctypes.c_double(self.starplanet.w), 
                                 ctypes.c_double(self.starplanet.P), 
                                 ctypes.c_double(self.starplanet.t0), 
                                 ctypes.c_double(self.starplanet.n))
        
        tp = t + tt - self.starplanet.t0
        
        xyz_sp = (ctypes.c_double * 6)(*np.zeros(6))
        keplib.find_xyz(xyz_sp, 
                        ctypes.c_double(tp), 
                        ctypes.c_double(self.starplanet.n), 
                        ctypes.c_double(self.starplanet.t0), 
                        ctypes.c_double(self.starplanet.e), 
                        ctypes.c_double(self.starplanet.a), 
                        ctypes.c_double(self.star.mass), 
                        ctypes.c_double(self.planet.mass), 
                        ctypes.c_double(self.starplanet.w), 
                        ctypes.c_double(self.starplanet.Omega), 
                        ctypes.c_double(self.starplanet.i))
        
        xyz_pm = (ctypes.c_double * 6)(*np.zeros(6))
        keplib.find_xyz(xyz_pm, 
                        ctypes.c_double(tp), 
                        ctypes.c_double(self.planetmoon.n), 
                        ctypes.c_double(self.planetmoon.t0), 
                        ctypes.c_double(self.planetmoon.e), 
                        ctypes.c_double(self.planetmoon.a), 
                        ctypes.c_double(self.planet.mass), 
                        ctypes.c_double(self.moon.mass), 
                        ctypes.c_double(self.planetmoon.w), 
                        ctypes.c_double(self.planetmoon.Omega), 
                        ctypes.c_double(self.planetmoon.i))
        
        xs, ys, zs, xp_tmp, yp_tmp, zp_tmp = np.array(xyz_sp)
        xp, yp, zp, xm, ym, zm = np.array(xyz_pm)

        xp = xp_tmp + xp
        yp = yp_tmp + yp
        zp = zp_tmp + zp
        xm = xp_tmp + xm
        ym = yp_tmp + ym
        zm = zp_tmp + zm
        
        return np.array([xs, ys, zs]), np.array([xp, yp, zp]), np.array([xm, ym, zm])
    
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
        
            st, pl, mo = self.coords(t)
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
        
        st, pl, mo = self.coords(t)
        stx, sty, stz = st / strad_au
        plx, ply, plz = pl / strad_au
        mox, moy, moz = mo / strad_au
        
        starpatch = plt.Circle((stx, sty), 1, animated=True, **stkwargs)
        planetpatch = plt.Circle((plx, ply), pp, animated=True, **plkwargs)
        moonpatch = plt.Circle((mox, moy), pm, animated=True, **mokwargs)

        ax.add_patch(starpatch)
        ax.add_patch(planetpatch)
        ax.add_patch(moonpatch)
        
        return ax