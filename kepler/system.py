from kepler import keplerian_system
import numpy as np
from scipy.optimize import minimize
from astropy import constants as ac

class system:
    
    days_in_year = 365.256
    earths_in_sun = 332946.08
    G = 39.478 # au ^ 3 / (yr ^2 * M_sun)
    G = G / (days_in_year ** 2) / earths_in_sun
    
    def __init__(self, star, planet, moon):
        self.star = star 
        self.star.mass = self.star.mass * ac.M_sun.value / ac.M_earth.value
        self.star.radius = self.star.radius * ac.R_sun.value / ac.R_earth.value
        self.planet = planet
        self.moon = moon
        
        self.starplanet = keplerian_system(star, planet)
        self.planetmoon = keplerian_system(planet, moon)
        
    def set_planet_orbit(self, t0=0, e=0, P=days_in_year, Omega=0, w=0, i=90):
        self.starplanet.set_orbital_parameters(t0, e, P, Omega, w, i)
        
    def set_moon_orbit(self, t0=0, e=0, P=days_in_year, Omega=0, w=0, i=90):
        self.planetmoon.set_orbital_parameters(t0, e, P, Omega, w, i)
        
    def coords(self, t):
        
        tp = t + self.starplanet.find_transit() - self.starplanet.t0
        
        starplanet_center = self.starplanet.barycentric_coords(tp)
        planetmoon_center = self.planetmoon.barycentric_coords(tp)
        starcoords = starplanet_center[0]
        planetcoords = starplanet_center[1] + planetmoon_center[0]
        mooncoords = starplanet_center[1] + planetmoon_center[1]
        
        return starcoords, planetcoords, mooncoords
    
    def sky_projected_distance(self, t):
        
        st, pl, mo = self.coords(t)
        
        stx, sty, stz = st
        plx, ply, plz = pl
        mox, moy, moz = mo
        
        starplanet_distance = np.sqrt((stx - plx)**2 + (sty + ply)**2)
        starmoon_distance = np.sqrt((stx - mox)**2 + (sty + moy)**2)
        planetmoon_distance = np.sqrt((plx - mox)**2 + (ply - moy)**2)
        
        return (starplanet_distance, starmoon_distance, planetmoon_distance)
    
    def starplanet_distance(self, t):
        
        tp = t + self.starplanet.find_transit() - self.starplanet.t0
        f = self.starplanet.f(tp)
        
        return (self.starplanet.a * (1 - self.starplanet.e**2) / 
                (1 + self.starplanet.e * np.cos(f)) * 
                np.sqrt(1 - (np.sin(self.starplanet.w + f) * np.sin(self.starplanet.i))**2))
    
    def find_contacts(self):
        
        mr = self.moon.radius * ac.R_earth.value / ac.au.value
        sr = self.star.radius * ac.R_earth.value / ac.au.value
        t_trans = self.starplanet.t0
        d = self.planetmoon.a + mr + sr
   
        g = lambda t: np.abs(self.starplanet_distance(t) - d)
        res = minimize(g, x0=self.starplanet.t0, method='Nelder-Mead').x[0]
        return (t_trans - np.abs(res - t_trans), t_trans + np.abs(res - t_trans))
    
    def reduce_t(self, t, factor=100):
        
        t_trans = self.starplanet.t0
        contact_duration = np.diff(self.find_contacts())
        is_transiting = np.isclose(t % self.starplanet.P - t_trans, np.zeros_like(t), atol = contact_duration)

        return is_transiting
        
        
    