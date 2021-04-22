#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "kepler.h"

#define PI 3.14159265358979323846
#define G 8.887641315935487e-10


// compute the mean intensity in the annulus overlapped by the occulting body 
// according to the small-body approximation from Mandel & Agol 
double intensity(double z, double p, double c1, double c2, double c3, double c4){
    
    if (p == 0){
        return 0.0;
    }
    else if (z > 1 + p){
        return 0.0;
    }
    else{
        if ((z < 1 + p) & (z > 1 - p)){
            double beta = - pow(p, 2.) + 2 * p * z - pow(z, 2.) + 1;
            double integral = 4. * pow(beta, 7./4.) * (c3 / 7.) +
                4 * pow(beta, 5./4.) * (c1 / 5.) + 
                2 * pow(beta, 3./2.) * (c2 / 3.) + 
                beta + c1 * (pow(p, 2.) - 
                2 * p * z + pow(z, 2.) - 1) + 
                c2 * (pow(p, 2.) - 2. * p * z + 
                pow(z, 2.) - 1) + c3 * (pow(p, 2.) - 
                2 * p * z + pow(z, 2.) - 1) + 
                c4 * (pow(p, 4.) / 2. - 2. * pow(p, 3.) * z + 
                3 * pow(p, 2.) * pow(z, 2.) - 2 * p * pow(z, 3.) + 
                pow(z, 4.) / 2. - 1. / 2.);
            return integral / (1 - pow(z - p, 2.));
        }
        else{
            double alpha = - pow(p, 2.) - 2 * p * z - pow(z, 2.) + 1;
            double beta = alpha + 4 * p * z;
            double integral = -4 * pow(alpha, 7./4.) * c3 / 7. - 
                4 * pow(alpha, 5./4.) * c1 / 5. - 
                2 * pow(alpha, 3./2.) * c2 / 3. + 
                4 * pow(beta, 7./4.) * c3 / 7. + 
                4 * pow(beta, 5./4.) * c1 / 5. + 
                2 * pow(beta, 3./2.) * c2 / 3. - 
                4 * c1 * p * z - 4 * c2 * p * z - 
                4 * c3 * p * z - 4 * c4 * pow(p, 3.) * z - 
                4 * c4 * p * pow(z, 3.) + 4. * p * z;
            return integral / (4. * z * p);
        }
    }
}

// overlap of two circles 
double overlap(double bigr, double r, double d){
    
    double k0 = acos((pow(d, 2.) + pow(r, 2.) - pow(bigr, 2.)) / (2 * d * r));
    double k1 = acos((pow(d, 2.) + pow(bigr, 2.) - pow(r, 2.)) / (2 * d * bigr));
    double k2 = sqrt(pow(2. * d * bigr, 2.) - pow(pow(bigr, 2.) + pow(d, 2.) - pow(r, 2.), 2.)) / 2.;
    return pow(r, 2.) * k0 + pow(bigr, 2.) * k1 - k2;
}

// compute the area occulted by the planet + moon according to Kipping, 2012 
// casenums are for debugging and do not correspond to the case labels from Kipping 
void area(double * a, double zp, double zm, double zpm, double pp, double pm, double t){
    
    //FILE *fp;
    //fp = fopen("./debug_phot.txt", "a");
    //FILE *fp_c;
    //fp_c = fopen("./debug_intersections.txt", "a");
    //FILE *fp_d;
    //fp_d = fopen("./debug_chords.txt", "a");
    double ap;
    double am;
    double alpha;
    double alpha1;
    double alpha2;
    double sumfact;
    double c14_1a_overlap;
    double c14_1b_overlap;
    double x12, y12, x13, y13, x23, y23;
    double x13_p, y13_p, x23_pp, y23_pp;
    double cos_theta_pp, sin_theta_pp;
    double d1, d2, d3;
    double sqrtarg;
    double cos_theta, sin_theta;
    
    int casenum = 0;
    
    if (zp >= 1 + pp){
        if (zm >= 1 + pm){
            ap = 0;
            am = 0;
        }
        else if ((1 - pm < zm) & (zm < 1 + pm)){
            alpha = overlap(1, pm, zm);
            ap = 0;
            am = alpha;
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
            alpha = overlap(1, pp, zp);
            ap = alpha;
            am = 0;
        }

        else if (((1 - pm) < zm) & (zm < (1 + pm))){
            casenum = 3;
            if (zpm >= (pm + pp)){
                casenum = 50;
                alpha1 = overlap(1, pp, zp);
                alpha2 = overlap(1, pm, zm);
                ap = alpha1;
                am = alpha2;
            }
            else if (((pp - pm) < zpm) & (zpm < (pp + pm))){
                casenum = 4;
                alpha = overlap(1., pp, zp);
                x12 = (1. - pow(pp, 2.) + pow(zp, 2.)) / (2 * zp);
                y12 = sqrt(2. * pow(zp, 2.) * (1. + pow(pp, 2.)) - 
                                  pow((1. - pow(pp, 2.)), 2.) - pow(zp, 4.)) / (2. * zp);
                cos_theta = (pow(zp, 2.) + pow(zm, 2.) - pow(zpm, 2.)) / (2. * zp * zm);
                sin_theta = sqrt(1. - pow(cos_theta, 2.));
                if ((pow((x12 - zm * cos_theta), 2.) + pow((y12 - zm * sin_theta), 2.)) < pow(pm, 2.)){
                    if ((pow((x12 - zm*cos_theta), 2.) + pow((y12 + zm * sin_theta), 2.)) < pow(pm, 2.)){
                        if (zp > 1){
                            casenum = 5;
                            ap = alpha;
                            am = (overlap(1, pm, zm) - overlap(1, pp, zp));
                        }
                        else{
                            casenum = 6;
                            ap = alpha;
                            am = pow(pp, 2.) * PI + 
                                (overlap(1, pm, zm) - 
                                 overlap(1, pp, zp) - 
                                 overlap(pp, pm, zpm));
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
                        sin_theta_pp = sqrt(1. - pow(cos_theta_pp, 2.));
                        x23 = x23_pp * cos_theta_pp - y23_pp * sin_theta_pp + zp;
                        y23 = x23_pp * sin_theta_pp + y23_pp * cos_theta_pp;
                        //fprintf(fp_c, "%f %f %f %f %f %f %f\n", t, x12, y12, x13, y13, x23, y23);
                        d1 = sqrt(pow((x12 - x13), 2.) + pow((y12 - y13), 2.));
                        d2 = sqrt(pow((x12 - x23), 2.) + pow((y12 - y23), 2.));
                        d3 = sqrt(pow((x13 - x23), 2.) + pow((y13 - y23), 2.));
                        //fprintf(fp_d, "%f %f %f\n", d1, d2, d3);
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
                            ap = alpha;
                            am = (overlap(1, pm, zm) - c14_1a_overlap);
                        }
                        else{
                            casenum = 9;
                            sumfact = asin(d1 / 2.);
                            sumfact += pow(pp, 2.) * asin(d2 / (2. * pp));
                            sumfact += pow(pm, 2.) * asin(d3 / (2. * pm));
                            c14_1b_overlap = sqrt((d1 + d2 + d3) * 
                                                         (d2 + d3 - d1) * 
                                                         (d1 + d3 - d2) * 
                                                         (d1 + d2 - d3)) / 4. + sumfact - 
                                d1 * sqrt(4. - pow(d1, 2.)) / 4. - 
                                d2 * sqrt(4. * pow(pp, 2.) - pow(d2, 2.)) / 4. + 
                                d3 * sqrt(4. * pow(pm, 2.) - pow(d3, 2.)) / 4.;
                            ap = alpha;
                            //fprintf(fp, "%f %f\n", overlap(1., pm, zm), c14_1b_overlap);
                            //fprintf(fp_d, "pp: %f, pm: %f, sumfact: %f, zp: %f, zm: %f, zpm: %f\n", pp, pm, sumfact, zp, zm, zpm);
                            am = (overlap(1., pm, zm) - c14_1b_overlap);
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
                        if ((zm - pm) < (zp - pp)){
                            ap = alpha;
                            am = pow(pm, 2) * PI - overlap(pp, pm, zpm);
                        }
                        else{
                            casenum = 11;
                            ap = alpha;
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
                            ap = alpha;
                            am = (overlap(1, pm, zm) - overlap(pp, pm, zpm));
                        }
                        else{
                            casenum = 14;
                            ap = alpha;
                            am = overlap(1, pm, zm);
                        }
                    }
                }
            }
            else{
                casenum = 15;
                alpha = overlap(1, pp, zp);
                ap = alpha;
                am = 0;
            }
        }
        else{
            casenum = 16;
            if (zpm >= pm + pp){
                casenum = 17;
                alpha = overlap(1, pp, zp);
                ap = alpha;
                am = pow(pm, 2) * PI;
            }
            else if ((pp - pm < zpm) & (zpm < pp + pm)){
                casenum = 18;
                alpha1 = overlap(1, pp, zp);
                alpha2 = overlap(pp, pm, zpm);
                ap = alpha1;
                am = pow(pm, 2) * PI - alpha2;
            }
            else{
                casenum = 19;
                alpha = overlap(1, pp, zp);
                ap = alpha;
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
                alpha = overlap(1, pm, zm);
                ap = pow(pp, 2) * PI;
                am = alpha;
            }
            else{
                casenum = 24;
                alpha1 = overlap(1, pm, zm);
                alpha2 = overlap(pp, pm, zpm);
                ap = pow(pp, 2) * PI;
                am = (alpha1 - alpha2);
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
                alpha = overlap(pp, pm, zpm);
                ap = pow(pp, 2) * PI;
                am = pow(pm, 2) * PI - alpha;
            }
            else{
                casenum = 27;
                ap = pow(pp, 2) * PI;
                am = 0;
            }
        }
    }
    //fprintf(fp, "%f %d %f %f\n", t, casenum, ap, am);
    //fclose(fp);
    //fclose(fp_c);
    //fclose(fp_d);
    
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

    FILE *fp;
    fp = fopen("./times.txt", "a");
    
    int j;
    double a[2];
    double arp;
    double arm;
    double Ip;
    double Im;
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
    dt = find_dt(wp, ep, em, strad_au, ap, am, ms, mp, mm);
    
    for (j = 0; j < m; j++){
        tj = t[j] + tt - t0p;
        if (((fabs(fmod(tj, P) - tt) < dt) | (fabs(fmod(tj, P) - tt) > (P - dt)))){
            
            fprintf(fp, "%f\n", t[j]);
            
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
            Ip = intensity(zp, pp, c1, c2, c3, c4);
            Im = intensity(zm, pm, c1, c2, c3, c4);
            flux[j] = - (arp * Ip + arm * Im) / (4. * Sigma * PI);
        }
        else{
            flux[j] = 0;
        }
    }
    fclose(fp);


}