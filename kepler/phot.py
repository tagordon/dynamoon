import numpy as np
from scipy.integrate import quad
import astropy.constants as ac

def I_quadratic(r, l):
    mu = np.sqrt(1 - r**2)
    l1, l2 = l
    return 1 - l1*(1 - mu) - l2*(1-mu)**2

def I_nonlinear(r, c):
    mu = np.sqrt(1 - r**2)
    return 1 - np.sum(np.array([c * (1 - mu ** (n/2)) for n, c in enumerate(c)]), axis=0)

def lc_no_overlap(z, p, c, I_func, analytic=True, epsabs=1e-10, epsrel=1e-10):
    
    if z > 1 + p:
        return 1
    else:
        c1, c2, c3, c4 = c
        c0 = 1 - np.sum(c)
        c = np.concatenate([[c0], c])
        Sigma = np.sum(np.array([c*(n+4)**(-1) for n, c in enumerate(c)]))
        integrand = lambda r, c1, c2, c3, c4: I_func(r, np.array([c1, c2, c3, c4])) * 2 * r
        if (z < 1 + p) & (z > 1 - p):
            alpha = -p**2 - 2 * p * z - z**2 + 1
            beta = -p**2 + 2 * p * z - z**2 + 1
            if analytic:
                integral = (4*beta**(7/4)*c3/7 + 
                        4*beta**(5/4)*c1/5 + 
                        2*beta**(3/2)*c2/3 + 
                        beta + 
                        c1*p**2 - 
                        2*c1*p*z + 
                        c1*z**2 - 
                        c1 + 
                        c2*p**2 - 
                        2*c2*p*z + 
                        c2*z**2 - 
                        c2 + 
                        c3*p**2 - 
                        2*c3*p*z + 
                        c3*z**2 - 
                        c3 + 
                        c4*p**4/2 - 
                        2*c4*p**3*z + 
                        3*c4*p**2*z**2 - 
                        2*c4*p*z**3 + 
                        c4*z**4/2 - 
                        c4/2)
            else:
                integral = quad(integrand, z-p, 1, args=c, epsabs=epsabs, epsrel=epsrel)[0]
            Iz = integral / (1 - (z - p)**2)
            f1 = (p**2) * np.arccos((z-1)/p)
            f2 = (z-1)*np.sqrt((p**2) - (z-1)**2)
            return 1 - (Iz / (np.pi * 4 * Sigma)) * (f1 - f2)
        elif(z <= 1 - p):
            alpha = -p**2 - 2 * p * z - z**2 + 1
            beta = -p**2 + 2 * p * z - z**2 + 1
            if analytic:
                integral = (-4*alpha**(7/4)*c3/7 - 
                        4*alpha**(5/4)*c1/5 - 
                        2*alpha**(3/2)*c2/3 + 
                        4*beta**(7/4)*c3/7 + 
                        4*beta**(5/4)*c1/5 + 
                        2*beta**(3/2)*c2/3 - 
                        4*c1*p*z - 
                        4*c2*p*z - 
                        4*c3*p*z - 
                        4*c4*p**3*z - 
                        4*c4*p*z**3 + 
                        4*p*z)
            else:
                integral = quad(integrand, z-p, z+p, args=c, epsabs=epsabs, epsrel=epsrel)[0]
            Iz = integral / (4 * z * p)
            return 1 - (p**2) * Iz / (4 * Sigma)
        
def lc(z, p, ld_coeffs, ld='quad', analytic=True, epsabs=1e-10, epsrel=1e-10):
    
    if ld == 'quad':
        l1, l2 = ld_coeffs
        c = np.array([0, l1 + 2*l2, 0, -l2])
        I_func = I_nonlinear
    elif ld == 'nonlinear':
        c = ld_coeffs
        I_func = I_nonlinear
    else:
        raise AttributeError('ld must be one of quad or nonlinear')
        
    return np.array([lc_no_overlap(z, p, c, I_func, analytic=analytic, epsabs=epsabs, epsrel=epsrel) for z in z])

def flux(system, t, ld_coeffs, ld='quad', analytic=True, epsabs=1e-10, epsrel=1e-10):
    
    pp = system.planet.radius / system.star.radius
    pm = system.moon.radius / system.star.radius
    strad_au = system.star.radius * ac.R_earth.value / ac.au.value
    
    mask = system.reduce_t(t, factor=100)
    rt = t[mask]
    zp, zm = system.sky_projected_distance(rt) / strad_au
    #zp = zp[zp < 1+pp]
    #zm = zm[zm < 1+pm]
    
    lcp = lc(zp, pp, ld_coeffs, ld=ld, analytic=analytic, epsabs=epsabs, epsrel=epsrel)
    lcm = lc(zm, pm, ld_coeffs, ld=ld, analytic=analytic, epsabs=epsabs, epsrel=epsrel)
    
    full_lcp = np.ones_like(t)
    full_lcm = np.ones_like(t)
    full_lcp[mask] = lcp
    full_lcm[mask] = lcm
    
    return full_lcp, full_lcm
    
    