#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "kepler.h"

#define PI 3.14159265358979323846
#define G 8.887641315935487e-10


// compute the mean intensity in the annulus overlapped by the occulting body 
// according to the small-body approximation from Mandel & Agol 
void intensity(double * res, double z, double p, double c1, double c2, double c3, double c4){
    
    if (p == 0){
        return;
    }
    else if (z > 1 + p){
        return;
    }
    else{
        if ((z < 1 + p) & (z > 1 - p)){
            double alpha = - pow(p, 2.) - 2 * p * z - pow(z, 2.) + 1;
            double beta = alpha + 4 * p * z;
            double pmz = p - z;
            double zmp = z - p;
            double I = c1 * ((4./5.) * pow(beta, 5./4.) - beta)
                + c2 * ((2./3.) * pow(beta, 3./2.) - beta)
                + c3 * ((4./7.) * pow(beta, 7./4.) - beta)
                + c4 * (pow(z - p, 4.) - 1) / 2. + beta;
            
            double Iz = 2 * c1 * pmz * (pow(beta, 1./4.) - 1)
                + 2 * c2 * pmz * (pow(beta, 1./2.) - 1)
                + 2 * c3 * pmz * (pow(beta, 3./4.) - 1)
                + 2 * c4 * pow(zmp, 3.) + 2 * pmz;
            
            double Ip = 2 * c1 * zmp * (pow(beta, 1./4.) - 1)
                + 2 * c2 * zmp * (pow(beta, 1./2.) - 1)
                + 2 * c3 * zmp * (pow(beta, 3./4.) - 1)
                + 2 * c4 * pow(zmp, 3.) + 2 * zmp;
            
            double Ic1 = c1 * ((4./5.) * pow(beta, 5./4.) - beta);
            double Ic2 = c2 * ((2./3.) * pow(beta, 3./2.) - beta);
            double Ic3 = c3 * ((4./7.) * pow(beta, 7./4.) - beta);
            double Ic4 = (pow(zmp, 4.) - 1) / 2;
            
            double denom = (1 - pow(z - p, 2.));
            double d_denom = 2 * (z - p);
            res[0] = I / denom;
            res[1] = (denom * Iz + I * d_denom) / pow(denom, 2);
            res[2] = (denom * Ip - I * d_denom) / pow(denom, 2);
            res[3] = Ic1 / denom;
            res[4] = Ic2 / denom;
            res[5] = Ic3 / denom;
            res[6] = Ic4 / denom;
            return;
        }
        else{
            
            double alpha = - pow(p, 2.) - 2 * p * z - pow(z, 2.) + 1;
            double beta = alpha + 4 * p * z;
            double pmz = p - z;
            double ppz = p + z;
            double I = (4./5.) * c1 * (pow(beta, 5./4.) - pow(alpha, 5./4.))
                + (2./3.) * c2 * (pow(beta, 3./2.) - pow(alpha, 3./2.))
                + (4./7.) * c3 * (pow(beta, 7./4.) - pow(alpha, (7./4.))) 
                + 4 * p * z * (1 - c1 - c2 - c3 - c4 * (pow(p, 2.) + pow(z, 2.)));
            
            double Iz = 2 * c1 * (pow(beta, 1./4.) * pmz + pow(alpha, 1./4.) * ppz)
                + 2 * c2 * (pow(beta, 1./2.) * pmz + pow(alpha, 1./2.) * ppz)
                + 2 * c3 * (pow(beta, 3./4.) * pmz + pow(alpha, 3./4.) * ppz)
                + 4 * p * (1 - c1 - c2 - c3) - 4 * c4 * (pow(p, 3.) + 3 * pow(z, 2.) * p);
                
            double Ip = 2 * c1 * (-pow(beta, 1./4.) * pmz + pow(alpha, 1./4.) * ppz)
                + 2 * c2 * (-pow(beta, 1./2.) * pmz + pow(alpha, 1./2.) * ppz)
                + 2 * c3 * (-pow(beta, 3./4.) * pmz + pow(alpha, 3./4.) * ppz)
                + 4 * z * (1 - c1 - c2 - c3) - 4 * c4 * (pow(z, 3.) + 3 * pow(p, 2.) * z);
                
            double Ic1 = (4./5.) * (pow(beta, 5./4.) - pow(alpha, 5./4.)) - 4 * p * z;
            double Ic2 = (2./3.) * (pow(beta, 3./2.) - pow(alpha, 3./2.)) - 4 * p * z;
            double Ic3 = (4./7.) * (pow(beta, 7./4.) - pow(alpha, 7./4.)) - 4 * p * z;
            double Ic4 = -4 * p * z * (pow(p, 2) + pow(z, 2));
            
            double denom = 4. * z * p;
            res[0] = I / denom;
            res[1] = (denom * Iz - I * 4 * p) / pow(denom, 2);
            res[2] = (denom * Ip - I * 4 * z) / pow(denom, 2);
            res[3] = Ic1 / denom;
            res[4] = Ic2 / denom;
            res[5] = Ic3 / denom;
            res[6] = Ic4 / denom;
            return;
        }
    }
}

// overlap of two circles 
void overlap(double * res, double bigr, double r, double d, bool computebigr){
    
    double alpha = pow(d, 2.) + pow(r, 2.) - pow(bigr, 2.);
    double beta = pow(d, 2.) + pow(bigr, 2.) - pow(r, 2.);
            
    double k0 = acos(alpha / (2 * d * r));
    double k1 = acos(beta / (2 * d * bigr));
    double k2 = sqrt(pow(2. * d * bigr, 2.) - pow(beta, 2.)) / 2.;
    res[0] = pow(r, 2.) * k0 + pow(bigr, 2.) * k1 - k2;
    
    res[1] = r * k0 - pow(r, 2.) / (d * sqrt(1 - pow(alpha, 2.))) 
        + r * bigr / (d * sqrt(1 - pow(beta, 2.)))
        - beta / (2 * k2);
    
    if (computebigr) {
        res[2] = 2 * bigr * k1 - pow(bigr, 2) / (d * sqrt(1 - pow(beta, 2)))
            - bigr * k1 + bigr * r / (d * sqrt(1 - pow(alpha, 2)))
            - beta * bigr / k2;
    }
    else{
        res[2] = 0;
    }
    
    res[3] = -r / sqrt(1 - pow(alpha, 2)) - pow(r, 2) * k0 / d 
        + bigr / sqrt(1 - pow(beta, 2)) - pow(bigr, 2) * k1 / d
        + d * beta / k2;
}

// compute the area occulted by the planet + moon according to Kipping, 2012 
// //casenums are for debugging and do not correspond to the case labels from Kipping 
void area(double * a, double zp, double zm, double zpm, double pp, double pm, double t){
    
    //FILE *fp;
    //fp = fopen("./debug_phot.txt", "a");
    //FILE *fp_c;
    //fp_c = fopen("./debug_intersections.txt", "a");
    //FILE *fp_d;
    //fp_d = fopen("./debug_chords.txt", "a");
    double ap;
    double am;
    double sumfact;
    double c14_1a_overlap;
    double c14_1b_overlap;
    double x12, y12, x13, y13, x23, y23;
    double x13_p, y13_p, x23_pp, y23_pp;
    double cos_theta_pp, sin_theta_pp;
    double d1, d2, d3;
    double sqrtarg;
    double cos_theta, sin_theta;
    double alpha[4];
    double alpha1[4];
    double alpha2[4];
    
    int casenum = 0;
    
    if (zp >= 1 + pp){
        if (zm >= 1 + pm){
            ap = 0;
            am = 0;
        }
        else if ((1 - pm < zm) & (zm < 1 + pm)){
            overlap(alpha, 1, pm, zm, false);
            ap = 0;
            am = alpha[0];
        }
        else{
            ap = 0;
            am = pow(pm, 2) * PI;
        }
    }
    
    else if (((1 - pp) < zp) & (zp < (1 + pp))){
        casenum = 1;
        if (zm >= 1 + pm){
            casenum = 2;
            overlap(alpha, 1, pp, zp, false);
            ap = alpha[0];
            am = 0;
        }

        else if (((1 - pm) < zm) & (zm < (1 + pm))){
            casenum = 3;
            if (zpm >= (pm + pp)){
                casenum = 50;
                overlap(alpha1, 1, pp, zp, false);
                overlap(alpha2, 1, pm, zm, false);
                ap = alpha1[0];
                am = alpha2[0];
            }
            else if (((pp - pm) < zpm) & (zpm < (pp + pm))){
                casenum = 4;
                overlap(alpha, 1., pp, zp, false);
                x12 = (1. - pow(pp, 2.) + pow(zp, 2.)) / (2 * zp);
                y12 = sqrt(2. * pow(zp, 2.) * (1. + pow(pp, 2.)) - 
                                  pow((1. - pow(pp, 2.)), 2.) - pow(zp, 4.)) / (2. * zp);
                cos_theta = (pow(zp, 2.) + pow(zm, 2.) - pow(zpm, 2.)) / (2. * zp * zm);
                if(cos_theta > 1){
                    sin_theta = 0.0;
                }
                else{
                    sin_theta = sqrt(1. - pow(cos_theta, 2.));
                }
                if ((pow((x12 - zm * cos_theta), 2.) + pow((y12 - zm * sin_theta), 2.)) < pow(pm, 2.)){
                    if ((pow((x12 - zm*cos_theta), 2.) + pow((y12 + zm * sin_theta), 2.)) < pow(pm, 2.)){
                        if (zp > 1){
                            casenum = 5;
                            ap = alpha[0];
                            overlap(alpha1, 1, pm, zm, false);
                            overlap(alpha2, 1, pp, zp, false);
                            am = (alpha1[0] - alpha2[0]);
                        }
                        else{
                            casenum = 6;
                            overlap(alpha1, 1, pm, zm, false);
                            overlap(alpha2, pp, pm, zpm, true);
                            ap = alpha[0];
                            am = pow(pp, 2.) * PI + alpha1[0] - ap - alpha2[0];
                        }
                    }
                    else{
                        casenum = 7;
                        x13_p = (1. - pow(pm, 2.) + pow(zm, 2.)) / (2. * zm);
                        y13_p = - sqrt(2. * pow(zm, 2.) * (1. + pow(pm, 2.)) - 
                                              pow(1. - pow(pm, 2.), 2.) - pow(zm, 4.)) / (2. * zm);
                        x13 = x13_p * cos_theta - y13_p * sin_theta;
                        y13 = x13_p * sin_theta + y13_p * cos_theta;
                        x23_pp = (pow(pp, 2.) - pow(pm, 2.) + pow(zpm, 2.)) / (2. * zpm);
                        y23_pp = sqrt(2. * pow(zpm, 2.) * (pow(pp, 2.) + pow(pm, 2.)) - 
                                             pow(pow(pp, 2.) - pow(pm, 2.), 2.) - pow(zpm, 4.)) / (2. * zpm);
                        cos_theta_pp = - (pow(zp, 2.) + pow(zpm, 2.) - pow(zm, 2.)) / (2. * zp * zpm);
                        if (cos_theta_pp > 1){
                            sin_theta_pp = 0.0;
                        }
                        else{
                            sin_theta_pp = sqrt(1. - pow(cos_theta_pp, 2.));
                        }
                        x23 = x23_pp * cos_theta_pp - y23_pp * sin_theta_pp + zp;
                        y23 = x23_pp * sin_theta_pp + y23_pp * cos_theta_pp;
                        d1 = sqrt(pow((x12 - x13), 2.) + pow((y12 - y13), 2.));
                        d2 = sqrt(pow((x12 - x23), 2.) + pow((y12 - y23), 2.));
                        d3 = sqrt(pow((x13 - x23), 2.) + pow((y13 - y23), 2.));
                        if ((zm * sin_theta) > (y13 + (zm * cos_theta - x13) * (y23 - y13) / (x23 - x13))){
                            casenum = 8;
                            sumfact = 0;
                            sumfact += asin(d1 / 2.) - sqrt(4 - pow(d1, 2.)) * d1 / 4.;
                            sumfact += pow(pp, 2) * asin(d2 / (2. * pp)) - 
                                sqrt(4 * pow(pp, 2) - pow(d2, 2.)) * d2 / 4.;
                            sumfact += pow(pm, 2) * asin(d3 / (2. * pm)) - 
                                sqrt(4 * pow(pm, 2) - pow(d3, 2.)) * d3 / 4.;
                            sqrtarg = (d1 + d2 + d3) * 
                                                      (d2 + d3 - d1) * 
                                                      (d1 + d3 - d2) * 
                                                      (d1 + d2 - d3);
                            c14_1a_overlap = (sqrt(sqrtarg) / 4. + sumfact);
                            overlap(alpha1, 1, pm, zm, false);
                            ap = alpha[0];
                            am = (alpha1[0] - c14_1a_overlap);
                        }
                        else{
                            casenum = 9;
                            sumfact = asin(d1 / 2.);
                            sumfact += pow(pp, 2.) * asin(d2 / (2. * pp));
                            sumfact -= pow(pm, 2.) * asin(d3 / (2. * pm));
                            c14_1b_overlap = sqrt((d1 + d2 + d3) * 
                                                         (d2 + d3 - d1) * 
                                                         (d1 + d3 - d2) * 
                                                         (d1 + d2 - d3)) / 4. + sumfact - 
                                d1 * sqrt(4. - pow(d1, 2.)) / 4. - 
                                d2 * sqrt(4. * pow(pp, 2.) - pow(d2, 2.)) / 4. + 
                                d3 * sqrt(4. * pow(pm, 2.) - pow(d3, 2.)) / 4. + PI * pow(pm, 2.);
                            overlap(alpha1, 1, pm, zm, false);
                            ap = alpha[0];
                            am = (alpha1[0] - c14_1b_overlap);
                        }
                    }
                }
                else{
                    casenum = 10;
                    x13_p = (1 - pow(pm, 2) + pow(zm, 2)) / (2 * zm);
                    y13_p = - sqrt(2 * pow(zm, 2) * (1 + pow(pm, 2)) - 
                                          pow((1 - pow(pm, 2)), 2) - pow(zm, 4)) / (2 * zm);
                    x13 = x13_p * cos_theta - y13_p * sin_theta;
                    y13 = x13_p * sin_theta + y13_p * cos_theta;
                    if ((pow((x13 - zp), 2) + pow(y13, 2)) < (pow(pp, 2))){
                        if ((zpm > (pp - pm)) & (zm < zp)){
                            overlap(alpha1, pp, pm, zpm, true);
                            ap = alpha[0];
                            am = pow(pm, 2) * PI - alpha1[0];
                        }
                        else{
                            casenum = 11;
                            ap = alpha[0];
                            am = 0;
                        }
                    }
                    else{
                        casenum = 12;
                        x23_pp = (pow(pp, 2) - pow(pm, 2) + pow(zpm, 2)) / (2 * zpm);
                        y23_pp = sqrt(2 * pow(zpm, 2) * (pow(pp, 2) + pow(pm, 2)) - 
                                         pow((pow(pp, 2) - pow(pm, 2)), 2) - pow(zpm, 4)) / (2 * zpm);
                        cos_theta_pp = - (pow(zp, 2) + pow(zpm, 2) - pow(zm, 2)) / (2 * zp * zpm);
                        sin_theta_pp = sqrt(1 - pow(cos_theta_pp, 2));
                        x23 = x23_pp * cos_theta_pp - y23_pp * sin_theta_pp + zp;
                        y23 = x23_pp * sin_theta_pp + y23_pp * cos_theta_pp;
                        if ((pow(x23, 2) + pow(y23, 2)) < 1){
                            casenum = 13;
                            overlap(alpha1, 1, pm, zm, false);
                            overlap(alpha2, pp, pm, zpm, true);
                            ap = alpha[0];
                            am = (alpha1[0] - alpha2[0]);
                        }
                        else{
                            casenum = 14;
                            overlap(alpha1, 1, pm, zm, false);
                            ap = alpha[0];
                            am = alpha1[0];
                        }
                    }
                }
            }
            else{
                casenum = 15;
                overlap(alpha, 1, pp, zp, false);
                ap = alpha[0];
                am = 0;
            }
        }
        else{
            casenum = 16;
            if (zpm >= pm + pp){
                casenum = 17;
                overlap(alpha, 1, pp, zp, false);
                ap = alpha[0];
                am = pow(pm, 2) * PI;
            }
            else if ((pp - pm < zpm) & (zpm < pp + pm)){
                casenum = 18;
                overlap(alpha1, 1, pp, zp, false);
                overlap(alpha2, pp, pm, zpm, true);
                ap = alpha1[0];
                am = pow(pm, 2) * PI - alpha2[0];
            }
            else{
                casenum = 19;
                overlap(alpha, 1, pp, zp, false);
                ap = alpha[0];
                am = 0;
            }
        }
    }
            
    else{
        casenum = 20;
        if (zm >= (1 + pm)){
            casenum = 21;
            ap = pow(pp, 2) * PI;
            am = 0;
        }
        else if ((zm > (1 - pm)) & (zm < (1 + pm))){
            casenum = 22;
            if (zpm >= pm + pp){
                casenum = 23;
                overlap(alpha, 1, pm, zm, false);
                ap = pow(pp, 2) * PI;
                am = alpha[0];
            }
            else{
                casenum = 24;
                overlap(alpha1, 1, pm, zm, false);
                overlap(alpha1, pp, pm, zpm, true);
                ap = pow(pp, 2) * PI;
                am = (alpha1[0] - alpha2[0]);
            }
        }

        else{
            casenum = 25;
            if (zpm >= pm + pp){
                ap = pow(pp, 2) * PI;
                am = pow(pm, 2) * PI;
            }
            else if ((pp - pm < zpm) & (zpm < pp + pm)){
                casenum = 26;
                overlap(alpha, pp, pm, zpm, true);
                ap = pow(pp, 2) * PI;
                am = pow(pm, 2) * PI - alpha[0];
            }
            else{
                casenum = 27;
                ap = pow(pp, 2) * PI;
                am = 0;
            }
        }
    }
    
    a[0] = ap;
    a[1] = am;
}

// compute the time between the start of the planet + moon transit and the center of transit using 
// Kipping, 2008
double find_dt(double w, double ep, double em, double rad, double ap, double am, double ms, double mp, double mm){
    
    double f_trans;
    double r_trans;
    double df;
    double fa;
    double dfa;
    double dfb;
    double h;
    
    f_trans = PI / 2. - w;
    r_trans = ap * (1 - pow(ep, 2)) / (1 + ep * cos(f_trans));
    df = asin((rad + (em + 1) * am) / r_trans);
    fa = f_trans - df;
    dfa = 2 * sqrt(1 - pow(ep, 2)) * atan(sqrt((1 - ep) / (1 + ep)) * tan(fa / 2.)) - 
        (ep * (1 - pow(ep, 2)) * sin(fa)) / (1 + ep * cos(fa));
    dfb = 2 * sqrt(1 - pow(ep, 2)) * atan(sqrt((1 - ep) / (1 + ep)) * tan(f_trans / 2.)) - 
        (ep * (1 - pow(ep, 2)) * sin(f_trans)) / (1 + ep * cos(f_trans));
    h = sqrt(G * (ms + mp + mm) * ap * (1 - pow(ep, 2)));
        
    return (dfb - dfa) * pow(ap, 2) / h;
}

// compute the light curve
void flux(double * flux, double * t, double pp, double pm, double c1, double c2, double c3, double c4, double t0p, double t0m, double P, double np, double ep, double ap, double mp, double ms, double mm, double wp, double omegap, double ip, double nm, double em, double am, double wm, double omegam, double im, double strad_au, int m){

    //FILE *fp;
    //fp = fopen("./times.txt", "a");
    
    int j;
    double a[2];
    double Ip[7] = {};
    double Im[7] = {};
    double arp;
    double arm;
    double xyz_ps[6];
    double xyz_mp[6];
    double xp;
    double yp;
    double xm;
    double ym;
    double zp;
    double zm;
    double zpm;
    double xs;
    double ys;
    double tt;
    double tj;
    double dt;
    
    double c0 = 1. - c1 - c2 - c3 - c4;
    double Sigma = c0 / 4.0;
    Sigma += c1 / 5.0;
    Sigma += c2 / 6.0;
    Sigma += c3 / 7.0;
    Sigma += c4 / 8.0;
    
    tt = find_transit(ep, wp, P, t0p, np);
    dt = find_dt(wp, ep, em, strad_au + 2 * pm, ap, am, ms, mp, mm);
    
    for (j = 0; j < m; j++){
        tj = t[j] + tt - t0p;
        if (((fabs(fmod(tj, P) - tt) < dt) | (fabs(fmod(tj, P) - tt) > (P - dt)))){
            
            //fprintf(fp, "%f\n", t[j]);
            
            find_xyz(xyz_ps, tj, np, t0p, ep, ap, ms, mp + mm, wp, omegap, ip);
            find_xyz(xyz_mp, tj, nm, t0m, em, am, mp, mm, wm, omegam, im);
            xs = xyz_ps[0];
            ys = xyz_ps[1];
            xp = xyz_ps[3] + xyz_mp[0];
            yp = xyz_ps[4] + xyz_mp[1];
            xm = xyz_ps[3] + xyz_mp[3];
            ym = xyz_ps[4] + xyz_mp[4];
            
            zp = sqrt(pow(xp - xs, 2) + pow(yp - ys, 2)) / strad_au;
            zm = sqrt(pow(xm - xs, 2) + pow(ym - ys, 2)) / strad_au;
            zpm = sqrt(pow(xp - xm, 2) + pow(yp - ym, 2)) / strad_au;
                
            area(a, zp, zm, zpm, pp, pm, t[j]);
            arp = *(a);
            arm = *(a+1);
            intensity(Ip, zp, pp, c1, c2, c3, c4);
            intensity(Im, zm, pm, c1, c2, c3, c4);
            flux[j] = - (arp * Ip[0] + arm * Im[0]) / (4. * Sigma * PI);
        }
        else{
            flux[j] = 0;
        }
    }
    //fclose(fp);


}