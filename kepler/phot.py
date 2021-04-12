import numpy as np
from scipy.integrate import quad
import astropy.constants as ac



area_dict = {1: ,
            2: ,
            3: ,
            4: ,
            5: ,
            6: ,
            7: ,
            8: ,
            9: ,
            10: ,
            11: ,
            12: ,
            13: ,
            14: ,
            15: ,
            16: ,
            17: ,
            18: ,
            19: ,
            20: ,
            21: ,
            22: ,
            23: ,
            24: ,
            25: ,
            26: ,
            27: }

def I_quadratic(r, l):
    mu = np.sqrt(1 - r**2)
    l1, l2 = l
    return 1 - l1*(1 - mu) - l2*(1-mu)**2

def I_nonlinear(r, c):
    mu = np.sqrt(1 - r**2)
    return 1 - np.sum(np.array([c * (1 - mu ** (n/2)) for n, c in enumerate(c)]), axis=0)

def flux(z, p, c, case):
    A = area(z, p, case)
    Iz = mean_intensity(z, p, c, case)
    newc = np.concatenate([[c0], c])
    Sigma = np.sum(np.array([c*(n+4)**(-1) for n, c in enumerate(newc)]))
    return -A * Iz / (4 * np.pi * Sigma)

def mean_intensity(z, p, c, case):
    pass

def area(z, p, case):
    pass

def find_area(z, p):
    zp, zm, zpm = z
    pm, pp = p
    
    if zp >= 1+pp:
        if zm >= 1+pm:
            return 0
        elif 1-pm < zm < z+pm:
            # alpha_S*
            return alpha
        else:
            return np.pi * (pm ** 2)
            
    if 1-p < zp < z+pp:
        if zm >= 1+pm:
            # alpha_P*
            return alpha
        elif 1-pm < zm < z+pm:
            if zpm >= pm+pp:
                # alpha_P* + alpha_S*
                return alpha1 + alpha2
            elif pp-pm < zpm < pp+pm:
                # determine sub-case
                return 14
            else:
                # alpha_P*
                return return alpha
        else:
            if zpm >= pm+pp:
                # alpha_P*
                return alpha + np.pi * (pm ** 2)
            elif pp-pm < zpm < pp+pm:
                # alpha_P*, alpha_PS
                return alpha1 + np.pi * (pm ** 2) - alpha2
            else:
                # alpha_P*
                return alpha
            
    else:
        if zm >= 1+pm:
            return np.pi * (pp ** 2)
        elif 1-pm < zm < z+pm:
            if zpm >= pm+pp:
                # alpha_S*
                return np.pi * (pp ** 2) + alpha
            else:
                # alpha_S*, alpha_SP
                return np.pi * (pp ** 2) + alpha1 - alpha2

        else:
            if zpm >= pm+pp:
                return np.pi * ((pp ** 2) + (pm ** 2))
            elif pp-pm < zpm < pp+pm:
                # alpha_SP
                return np.pi * (pp ** 2) + np.pi * (pm ** 2) - alpha
            else:
                return np.pi * (pp ** 2)   

def lc_no_overlap(zp, pp, zm, pm, c):
    
    flux = 0
    for z, p in zip([zp, zm], [pp, pm]):
    
        if z > 1 + p:
            pass
        else:
            c1, c2, c3, c4 = c
            c0 = 1 - np.sum(c)
            newc = np.concatenate([[c0], c])
            Sigma = np.sum(np.array([c*(n+4)**(-1) for n, c in enumerate(newc)]))
            if (z < 1 + p) & (z > 1 - p):
                alpha = -p**2 - 2 * p * z - z**2 + 1
                beta = -p**2 + 2 * p * z - z**2 + 1
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
                Iz = integral / (1 - (z - p)**2)
                f1 = (p**2) * np.arccos((z-1)/p)
                f2 = (z-1)*np.sqrt((p**2) - (z-1)**2)
                flux += - (Iz / (np.pi * 4 * Sigma)) * (f1 - f2)
            elif(z <= 1 - p):
                alpha = -p**2 - 2 * p * z - z**2 + 1
                beta = -p**2 + 2 * p * z - z**2 + 1
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
                Iz = integral / (4 * z * p)
            flux += - (p**2) * Iz / (4 * Sigma)
    return flux
        
def lc(zp, pp, zm, pm, ld_coeffs, ld='quad'):
    
    if ld == 'quad':
        l1, l2 = ld_coeffs
        c = np.array([0, l1 + 2*l2, 0, -l2])
        I_func = I_nonlinear
    elif ld == 'nonlinear':
        c = ld_coeffs
        I_func = I_nonlinear
    else:
        raise AttributeError('ld must be one of quad or nonlinear')
        
    return np.array([lc_no_overlap(zp, pp, zm, pm, c) for zp, zm in zip(zp, zm)])

def flux_dep(system, t, ld_coeffs, ld='quad'):
    
    pp = system.planet.radius / system.star.radius
    pm = system.moon.radius / system.star.radius
    strad_au = system.star.radius * ac.R_earth.value / ac.au.value
    
    mask = system.reduce_t(t, factor=100)
    rt = t[mask]
    zp, zm, zpm = system.sky_projected_distance(rt) / strad_au
    
    lcpm = lc(zp, pp, zm, pm, ld_coeffs, ld=ld)
    full_lc = np.zeros_like(t)
    full_lc[mask] = lcpm
    
    #lcp = lc(zp, pp, ld_coeffs, ld=ld)
    #lcm = lc(zm, pm, ld_coeffs, ld=ld)
    
    #full_lcp = np.ones_like(t)
    #full_lcm = np.ones_like(t)
    #full_lcp[mask] = lcp
    #full_lcm[mask] = lcm
    
    return full_lc
    
    