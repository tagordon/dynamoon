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

def mean_intensity(z, p, c):
    if p == 0:
        return 0
    if z > 1 + p:
        return 0
    else:
        c1, c2, c3, c4 = c
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
            return integral / (1 - (z - p)**2)
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
            return integral / (4 * z * p)
            
            
        #if (z < (1 + p)) & (z > (1 - p)):
        #    integral = ((4 * a - 5) * c1 * (a ** 4) / 5 + 
        #                (-42  + (42 - 28 * a ** 2) * c2 + 
        #                 6 * (7 - 4 * a **3) * c3 + 
        #                 21 * (1 + (z-p)**2)*c4) / 42)
        #    return integral / (1 - (z - p)**2) 
        #elif z < p:
        #    integral = ((4 / 5) * ((1 - b) + (z + p)**2 * (b - (5 / 4))) * c1 + 
        #               (2 / 3) * (1 - b**2 + (z + p)**2 * (b**2 - (3 / 2))) * c2 + 
        #               (1 / 14) * (14 * (z + p)**2 - 7 * (z + p)**4 * c4 + 
        #                           2 * (4 - 4 * b**3 + (z + p)**2 *(4 * b**3 - 7)) * c3))
        #    return integral / (z + p)**2
        #else:
        #    integral = ((4 / 5) * ((p**2) * (b - a) + 
        #                           (z**2 - 1) * (b - a) + 
        #                           2 * p * z * (a + b - 5 / 2)) * c1 + 
        #               (2 / 21) * (7 * (p ** 2 * (b**2 - a**2) 
        #                                + (z**2 - 1) * (b**2 - a**2) + 
        #                                2 * z * p * (a**2 + b ** 2 - 3)) * c2) - 
        #               6 * (b**3 - a**3 + 
        #                    2 * z * p * (-a**3 - b**3 + 7 / 2) + 
        #                    p**2 * (a**3 - b**3) + 
        #                    z**2 * (a**3 - b**3)) * c3 + 
        #               7 * z * p * (c4 * (z**2 + p**2) - 1))
        #    return integral / ((z + p)**2 - (z - p)**2)

def overlap(bigr, r, d):
    k0 = np.arccos((d**2 + r**2 - bigr**2) / (2*d*r))
    k1 = np.arccos((d**2 + bigr**2 - r**2) / (2*d*bigr))
    k2 = np.sqrt((2*d*bigr)**2 - (bigr**2 + d**2 - r**2)**2) / 2
    return (r**2) * k0 + (bigr**2) * k1 - k2

def area(z, p, t):
    zp, zm, zpm = z
    pp, pm = p
    
    if zp >= 1+pp:
        #print('a')
        if zm >= 1+pm:
            #print('b')
            return 0, 0
        elif 1-pm < zm < 1+pm:
            #print('c')
            alpha = overlap(1, pm, zm)
            return 0, alpha / np.pi 
        else:
            #print('d')
            return 0, pm ** 2
            
    if 1-pp < zp < 1+pp:
        #print('e')
        if zm >= 1+pm:
            #print('f')
            # alpha_P*
            alpha = overlap(1, pp, zp)
            return alpha / np.pi, 0
        elif 1-pm < zm < 1+pm:
            #print('g')
            if zpm >= pm+pp:
                #print('h')
                # alpha_P* + alpha_S*
                alpha1 = overlap(1, pp, zp)
                alpha2 = overlap(1, pm, zm)
                return alpha1 / np.pi, alpha2 / np.pi
            elif pp-pm < zpm < pp+pm:
                #print('i')
                alpha = overlap(1, pp, zp)
                # determine sub-case for case 14
                x12 = (1 - pp ** 2 + zp ** 2) / 2 / (zp ** 2)
                y12 = np.sqrt(2 * zp ** 2 * (1 + pp ** 2) - (1 - pp ** 2) ** 2 - zp ** 4) / 2 / (zp ** 2)
                cos_theta = (zp ** 2 + zm ** 2 - zpm**2) / (2 * zp * zm)
                sin_theta = np.sqrt(1 - cos_theta ** 2)
                if ((x12 - zm*cos_theta) ** 2 + (y12 - zm * sin_theta) ** 2) < pm ** 2:
                    #print('j')
                    if ((x12 - zm*cos_theta) ** 2 + (y12 + zm * sin_theta) ** 2) < pm ** 2:
                        #print('k')
                        if zp > 1:
                            #print('l')
                            return alpha / np.pi, (overlap(1, pm, zm) - overlap(1, pp, zp)) / np.pi
                        else:
                            #print('m')
                            return alpha / np.pi, pp ** 2 - (overlap(1, pm, zm) - overlap(1, pp, zp) - overlap(pp, pm, zm)) / np.pi
                    else:
                        #print('n')
                        x13_p = (1 - pm ** 2 + zm ** 2) / (2 * zm)
                        y13_p = - np.sqrt(2 * zm ** 2 * (1 + pm ** 2) - (1 - pm ** 2)**2 - zm ** 4) / (2 * zm)
                        x13 = x13_p * cos_theta - y13_p * sin_theta
                        y13 = x13_p * sin_theta + y13_p * cos_theta
                        x23_pp = (pp ** 2 - pm ** 2 + zpm ** 2) / (2 * zpm)
                        y23_pp = np.sqrt(2 * zpm ** 2 * (pp ** 2 + pm ** 2) - 
                                         (pp ** 2 - pm ** 2) ** 2 - zpm ** 4) / (2 * zpm)
                        cos_theta_pp = - (zp ** 2 + zpm ** 2 - zm ** 2) / (2 * zp * zpm)
                        sin_theta_pp = np.sqrt(1 - cos_theta_pp ** 2)
                        x23 = x23_pp * cos_theta_pp - y23_pp * sin_theta_pp + zp
                        y23 = x23_pp * sin_theta_pp + y23_pp * cos_theta_pp
                        if (zp * sin_theta) > (y13 + (zm * cos_theta - x13) * (y23 - y13) / (x23 - x13)):
                            #print('o')
                            c1 = (x12 - x13) ** 2 + (y12 - y13) ** 2
                            c2 = (x12 - x23) ** 2 + (y12 - y23) ** 2
                            c3 = (x13 - x23) ** 2 + (y13 - y23) ** 2
                            R1 = 1
                            R2 = pp
                            R3 = pm
                            c14_1a_overlap = (np.sqrt((c1 + c2 + c3) * 
                                                      (c2 + c3 - c1) * 
                                                      (c1 + c3 - c2) * 
                                                      (c1 + c2 - c3)) 
                                              / 4 + np.sum(np.array([r * np.arcsin(c / 2 / r) - 
                                                                     np.sqrt(4 * r ** 2 - c ** 2) * 
                                                                     c / 4 for (c, r) in 
                                                                     zip([c1, c2, c3], [1, pp, pm])])))
                            return alpha / np.pi, (overlap(1, pm, zm) - c14_1a_overlap) / np.pi
                        else:
                            #print('p')
                            c14_1b_overlap = (np.sqrt((c1 + c2 + c3) * 
                                                      (c2 + c3 - c1) * 
                                                      (c1 + c3 - c2) * 
                                                      (c1 + c2 - c3)) 
                                              / 4 + np.sum(np.array([r * np.arcsin(c / 2 / r) for (c, r) in 
                                                                     zip([c1, c2, c3], [1, pp, pm])]))
                                             - c1 * np.sqrt(4 - c1 ** 2) / 4 
                                             - c2 * np.sqrt(4 * pp ** 2 - c2 ** 2) / 4
                                             - c3 * np.sqrt(4 * pm ** 2 - c3 ** 2) / 4)
                            return alpha / np.pi, (overlap(1, pm, zm) - c14_1b_overlap) / np.pi
                
                else:
                    #print('q')
                    x13_p = (1 - pm ** 2 + zm ** 2) / (2 * zm)
                    y13_p = - np.sqrt(2 * zm ** 2 * (1 + pm ** 2) - (1 - pm ** 2)**2 - zm ** 4) / (2 * zm)
                    x13 = x13_p * cos_theta - y13_p * sin_theta
                    y13 = x13_p * sin_theta + y13_p * cos_theta
                    if ((x13 - zp) ** 2 + y13 ** 2) < (pp ** 2):
                        #print('r')
                        if (zm - pm) < (zp - pp):
                            #print('s')
                            return alpha / np.pi, (pm ** 2) - overlap(pp, pm, zpm) / np.pi
                        else:
                            #print('t')
                            return 0, 0
                    else:
                        #print('u')
                        x23_pp = (pp ** 2 - pm ** 2 + zpm ** 2) / (2 * zpm)
                        y23_pp = np.sqrt(2 * zpm ** 2 * (pp ** 2 + pm ** 2) - 
                                         (pp ** 2 - pm ** 2) ** 2 - zpm ** 4) / (2 * zpm)
                        cos_theta_pp = - (zp ** 2 + zpm ** 2 - zm ** 2) / (2 * zp * zpm)
                        sin_theta_pp = np.sqrt(1 - cos_theta_pp ** 2)
                        x23 = x23_pp * cos_theta_pp - y23_pp * sin_theta_pp + zp
                        y23 = x23_pp * sin_theta_pp + y23_pp * cos_theta_pp
                        if (x23 ** 2 + y23 ** 2) < 1:
                            #print('v')
                            return alpha / np.pi, (overlap(1, pm, zm) - overlap(pp, pm, zpm)) / np.pi
                        else:
                            #print('w')
                            return alpha / np.pi, overlap(1, pm, zm) / np.pi
                
            else:
                #print('x')
                # alpha_P*
                alpha = overlap(1, pp, zp)
                return alpha / np.pi, 0
        else:
            #print('y')
            if zpm >= pm+pp:
                #print('z')
                # alpha_P*
                alpha = overlap(1, pp, zp)
                return alpha / np.pi,  pm ** 2
            elif pp-pm < zpm < pp+pm:
                #print('1')
                # alpha_P*, alpha_PS
                alpha1 = overlap(1, pp, zp)
                alpha2 = overlap(pp, pm, zpm)
                return alpha1 / np.pi,  (pm ** 2) - alpha2 / np.pi
            else:
                #print('2')
                # alpha_P*
                alpha = overlap(1, pp, zp)
                return alpha / np.pi, 0
            
    else:
        #print('3')
        if zm >= 1+pm:
            #print('4')
            return pp ** 2, 0
        elif 1-pm < zm < 1+pm:
            #print('5')
            if zpm >= pm+pp:
                #print('6')
                # alpha_S*
                alpha = overlap(1, pm, zm)
                return pp ** 2,  alpha / np.pi
            else:
                #print('7')
                # alpha_S*, alpha_SP
                alpha1 = overlap(1, pm, zm)
                alpha2 = overlap(pp, pm, zpm)
                return pp ** 2,  (alpha1 - alpha2) / np.pi

        else:
            #print('8')
            if zpm >= pm+pp:
                #print('9')
                return pp ** 2, pm ** 2
            elif pp-pm < zpm < pp+pm:
                #print('0')
                # alpha_SP
                alpha = overlap(pp, pm, zpm)
                return pp ** 2, pm ** 2 - alpha / np.pi
            else:
                #print('!')
                return pp ** 2, 0 
            
def flux(z, p, c, t):
    Ap, Am = area(z, p, t)
    Izp = mean_intensity(z[0], p[0], c)
    Izm = mean_intensity(z[1], p[1], c)
    newc = np.concatenate([[1-np.sum(c)], c])
    Sigma = np.sum(np.array([c*(n+4)**(-1) for n, c in enumerate(newc)]))
    return -(Ap * Izp + Am * Izm) / (4 * Sigma)

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
        
def lc(system, t, ld_coeffs, ld='quad'):
    
    if ld == 'quad':
        l1, l2 = ld_coeffs
        c = np.array([0, l1 + 2*l2, 0, -l2])
        I_func = I_nonlinear
    elif ld == 'nonlinear':
        c = ld_coeffs
        I_func = I_nonlinear
    else:
        raise AttributeError('ld must be one of quad or nonlinear')
                                    
    pp = system.planet.radius / system.star.radius
    pm = system.moon.radius / system.star.radius
    strad_au = system.star.radius * ac.R_earth.value / ac.au.value
    p = (pp, pm)
                                    
    mask = system.reduce_t(t, factor=100)
    rt = t[mask]
    z = system.sky_projected_distance(rt) / strad_au
    f = np.zeros_like(t)
    f[mask] = np.array([flux(z, p, c, t) for (z, t) in zip(z.T, rt)])
    return f

def check_intensity_profile(system, t, ld_coeffs, ld='quad'):
    
    if ld == 'quad':
        l1, l2 = ld_coeffs
        c = np.array([0, l1 + 2*l2, 0, -l2])
        I_func = I_nonlinear
    elif ld == 'nonlinear':
        c = ld_coeffs
        I_func = I_nonlinear
    else:
        raise AttributeError('ld must be one of quad or nonlinear')
                                    
    pp = system.planet.radius / system.star.radius
    pm = system.moon.radius / system.star.radius
    strad_au = system.star.radius * ac.R_earth.value / ac.au.value
    p = (pp, pm)
                                    
    mask = system.reduce_t(t, factor=100)
    rt = t[mask]
    z = system.sky_projected_distance(rt) / strad_au
    f_pl = np.zeros_like(t)
    f_mo = np.zeros_like(t)
    f_pl[mask] = np.array([mean_intensity(z[0], p[0], c) for z in z.T])
    f_mo[mask] = np.array([mean_intensity(z[1], p[1], c) for z in z.T])
    return f_pl, f_mo

def check_occulted_area(system, t, ld_coeffs, ld='quad'):
    
    if ld == 'quad':
        l1, l2 = ld_coeffs
        c = np.array([0, l1 + 2*l2, 0, -l2])
        I_func = I_nonlinear
    elif ld == 'nonlinear':
        c = ld_coeffs
        I_func = I_nonlinear
    else:
        raise AttributeError('ld must be one of quad or nonlinear')
                                    
    pp = system.planet.radius / system.star.radius
    pm = system.moon.radius / system.star.radius
    strad_au = system.star.radius * ac.R_earth.value / ac.au.value
    p = (pp, pm)
                                    
    mask = system.reduce_t(t, factor=100)
    rt = t[mask]
    z = system.sky_projected_distance(rt) / strad_au
    f_pl = np.zeros_like(t)
    f_mo = np.zeros_like(t)
    f_pl[mask] = np.array([area(z, p, t)[0] for (z, t) in zip(z.T, rt)])
    f_mo[mask] = np.array([area(z, p, t)[1] for (z, t) in zip(z.T, rt)])
    return f_pl, f_mo

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
    
    