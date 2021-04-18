#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265358979323846

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
                beta + c1 * pow(p, 2.) - 
                2 * c1 * p * z + c1 * pow(z, 2.) - c1 + 
                c2 * pow(p, 2.) - 2. * c2 * p * z + 
                c2 * pow(z, 2.) - c2 + c3 * pow(p, 2.) - 
                2 * c3 * p * z + c3 * pow(z, 2.) - c3 + 
                c4 * pow(p, 4.) / 2. - 2. * c4 * pow(p, 3.) * z + 
                3 * c4 * pow(p, 2.) * pow(z, 2.) - 2 * c4 * p * pow(z, 3.) + 
                c4 * pow(z, 4.) / 2. - c4 / 2.;
            return integral / (1 - pow(z - p, 2.));
        }
        else{
            double alpha = - pow(p, 2.) - 2 * p * z - pow(z, 2.) + 1;
            double beta = - pow(p, 2.) + 2 * p * z - pow(z, 2.) + 1;
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
    
double overlap(double bigr, double r, double d){
    
    double k0 = acos((pow(d, 2) + pow(r, 2) - pow(bigr, 2)) / (2 * d * r));
    double k1 = acos((pow(d, 2) + pow(bigr, 2) - pow(r, 2)) / (2 * d * bigr));
    double k2 = sqrt(pow(2 * d * bigr, 2) - pow(pow(bigr, 2) + pow(d, 2) - pow(r, 2), 2)) / 2.;
    return pow(r, 2) * k0 + pow(bigr, 2) * k1 - k2;
}

double * area(double zp, double zm, double zpm, double pp, double pm){
    
    double ap;
    double am;
    double alpha;
    double alpha1;
    double alpha2;
    
    //FILE *fp;
    //fp = fopen("./debug_cases.txt", "a");
    int casenum = 0;
    
    if (zp >= 1 + pp){
        if (zm >= 1 + pm){
            ap = 0;
            am = 0;
        }
        else if ((1 - pm < zm) & (zm < 1 + pm)){
            alpha = overlap(1, pm, zm);
            ap = 0;
            am = alpha / PI;
        }
        else{
            ap = 0;
            am = pow(pm, 2);
        }
    }
    
    else if ((1 - pp < zp) & (zp < 1 + pp)){
        casenum = 1;
        if (zm >= 1 + pm){
            casenum = 2;
            alpha = overlap(1, pp, zp);
            ap = alpha / PI;
            am = 0;
        }

        else if ((1 - pm < zm) & (zm < 1 + pm)){
            casenum = 3;
            if (zpm >= pm + pp){
                casenum = 50;
                alpha1 = overlap(1, pp, zp);
                alpha2 = overlap(1, pm, zm);
                ap = alpha1 / PI;
                am = alpha2 / PI;
            }
            else if ((pp - pm < zpm) & (zpm < pp+pm)){
                casenum = 4;
                alpha = overlap(1, pp, zp);
                double x12 = (1 - pow(pp, 2.) + pow(zp, 2.)) / 2. / pow(zp, 2.);
                double y12 = sqrt(2 * pow(zp, 2.) * (1 + pow(pp, 2.)) - 
                                  pow((1 - pow(pp, 2.)), 2.) - pow(zp, 4.)) / 2. / pow(zp , 2.);
                double cos_theta = (pow(zp, 2.) + pow(zm, 2.) - pow(zpm, 2.)) / (2. * zp * zm);
                double sin_theta = sqrt(1 - pow(cos_theta, 2.));
                if ((pow((x12 - zm * cos_theta), 2.) + pow((y12 - zm * sin_theta), 2.)) < pow(pm, 2.)){
                    if ((pow((x12 - zm*cos_theta), 2.) + pow((y12 + zm * sin_theta), 2.)) < pow(pm, 2.)){
                        if (zp > 1){
                            casenum = 5;
                            ap = alpha / PI;
                            am = (overlap(1, pm, zm) - overlap(1, pp, zp)) / PI;
                        }
                        else{
                            casenum = 6;
                            ap = alpha / PI;
                            am = pow(pp, 2.) + 
                                (overlap(1, pm, zm) - 
                                 overlap(1, pp, zp) - 
                                 overlap(pp, pm, zpm)) / PI;
                        }
                    }
                    else{
                        casenum = 7;
                        double x13_p = (1 - pow(pm, 2) + pow(zm, 2)) / (2 * zm);
                        double y13_p = - sqrt(2 * pow(zm, 2) * (1 + pow(pm, 2)) - 
                                              pow(1 - pow(pm, 2), 2) - pow(zm, 4)) / (2 * zm);
                        double x13 = x13_p * cos_theta - y13_p * sin_theta;
                        double y13 = x13_p * sin_theta + y13_p * cos_theta;
                        double x23_pp = (pow(pp, 2) - pow(pm, 2) + pow(zpm, 2)) / (2 * zpm);
                        double y23_pp = sqrt(2 * pow(zpm, 2) * (pow(pp, 2) + pow(pm, 2)) - 
                                             pow(pow(pp, 2) - pow(pm, 2), 2) - pow(zpm, 4)) / (2 * zpm);
                        double cos_theta_pp = - (pow(zp, 2) + pow(zpm, 2) - pow(zm, 2)) / (2 * zp * zpm);
                        double sin_theta_pp = sqrt(1 - pow(cos_theta_pp, 2));
                        double x23 = x23_pp * cos_theta_pp - y23_pp * sin_theta_pp + zp;
                        double y23 = x23_pp * sin_theta_pp + y23_pp * cos_theta_pp;
                        double d1 = sqrt(pow((x12 - x13), 2.) + pow((y12 - y13), 2.));
                        double d2 = sqrt(pow((x12 - x23), 2.) + pow((y12 - y23), 2.));
                        double d3 = sqrt(pow((x13 - x23), 2.) + pow((y13 - y23), 2.));
                        if ((zp * sin_theta) > (y13 + (zm * cos_theta - x13) * (y23 - y13) / (x23 - x13))){
                            casenum = 8;
                            double sumfact = 0;
                            sumfact += asin(d1 / 2.) - sqrt(4 - pow(d1, 2.)) * d1 / 4.;
                            sumfact += pow(pp, 2) * asin(d2 / 2. / pp) - 
                                sqrt(4 * pow(pp, 2) - pow(d2, 2.)) * d2 / 4.;
                            sumfact += pow(pm, 2) * asin(d3 / 2. / pm) - 
                                sqrt(4 * pow(pm, 2) - pow(d3, 2.)) * d3 / 4.;
                            double sqrtarg = (d1 + d2 + d3) * 
                                                      (d2 + d3 - d1) * 
                                                      (d1 + d3 - d2) * 
                                                      (d1 + d2 - d3);
                            double c14_1a_overlap = (sqrt(sqrtarg) / 4. + sumfact);
                            ap = alpha / PI;
                            am = (overlap(1, pm, zm) - c14_1a_overlap) / PI;
                        }
                        else{
                            casenum = 9;
                            double sumfact = 0;
                            sumfact += asin(d1 / 2.);
                            sumfact += pow(pp, 2) * asin(d2 / 2. / pp);
                            sumfact += pow(pm, 2) * asin(d3 / 2. / pm);
                            double c14_1b_overlap = sqrt((d1 + d2 + d3) * 
                                                         (d2 + d3 - d1) * 
                                                         (d1 + d3 - d2) * 
                                                         (d1 + d2 - d3)) / 4. + sumfact - 
                                d1 * sqrt(4 - pow(d1, 2)) / 4. - 
                                d2 * sqrt(4 * pow(pp, 2) - pow(d2, 2)) / 4. - 
                                d3 * sqrt(4 * pow(pm, 2) - pow(d3, 2)) / 4.;
                            ap = alpha / PI;
                            am = (overlap(1, pm, zm) - c14_1b_overlap) / PI;
                        }
                    }
                }
                else{
                    casenum = 10;
                    double x13_p = (1 - pow(pm, 2) + pow(zm, 2)) / (2 * zm);
                    double y13_p = - sqrt(2 * pow(zm, 2) * (1 + pow(pm, 2)) - 
                                          pow((1 - pow(pm, 2)), 2) - pow(zm, 4)) / (2 * zm);
                    double x13 = x13_p * cos_theta - y13_p * sin_theta;
                    double y13 = x13_p * sin_theta + y13_p * cos_theta;
                    if ((pow((x13 - zp), 2) + pow(y13, 2)) < (pow(pp, 2))){
                        if ((zm - pm) < (zp - pp)){
                            ap = alpha / PI;
                            am = pow(pm, 2) - overlap(pp, pm, zpm) / PI;
                        }
                        else{
                            casenum = 11;
                            ap = 0;
                            am = 0;
                        }
                    }
                    else{
                        casenum = 12;
                        double x23_pp = (pow(pp, 2) - pow(pm, 2) + pow(zpm, 2)) / (2 * zpm);
                        double y23_pp = sqrt(2 * pow(zpm, 2) * (pow(pp, 2) + pow(pm, 2)) - 
                                         pow((pow(pp, 2) - pow(pm, 2)), 2) - pow(zpm, 4)) / (2 * zpm);
                        double cos_theta_pp = - (pow(zp, 2) + pow(zpm, 2) - pow(zm, 2)) / (2 * zp * zpm);
                        double sin_theta_pp = sqrt(1 - pow(cos_theta_pp, 2));
                        double x23 = x23_pp * cos_theta_pp - y23_pp * sin_theta_pp + zp;
                        double y23 = x23_pp * sin_theta_pp + y23_pp * cos_theta_pp;
                        if ((pow(x23, 2) + pow(y23, 2)) < 1){
                            casenum = 13;
                            ap = alpha / PI;
                            am = (overlap(1, pm, zm) - overlap(pp, pm, zpm)) / PI;
                        }
                        else{
                            casenum = 14;
                            ap = alpha / PI;
                            am = overlap(1, pm, zm) / PI;
                        }
                    }
                }
            }
            else{
                casenum = 15;
                alpha = overlap(1, pp, zp);
                ap = alpha / PI;
                am = 0;
            }
        }
        else{
            casenum = 16;
            if (zpm >= pm + pp){
                casenum = 17;
                alpha = overlap(1, pp, zp);
                ap = alpha / PI;
                am = pow(pm, 2);
            }
            else if ((pp - pm < zpm) & (zpm < pp + pm)){
                casenum = 18;
                alpha1 = overlap(1, pp, zp);
                alpha2 = overlap(pp, pm, zpm);
                ap = alpha1 / PI;
                am = pow(pm, 2) - alpha2 / PI;
            }
            else{
                casenum = 19;
                alpha = overlap(1, pp, zp);
                ap = alpha / PI;
                am = 0;
            }
        }
    }
            
    else{
        casenum = 20;
        if (zm >= 1 + pm){
            casenum = 21;
            ap = pow(pp, 2);
            am = 0;
        }
        else if ((1 - pm < zm) & (zm < 1 + pm)){
            casenum = 22;
            if (zpm >= pm + pp){
                casenum = 23;
                alpha = overlap(1, pm, zm);
                ap = pow(pp, 2);
                am = alpha / PI;
            }
            else{
                casenum = 24;
                alpha1 = overlap(1, pm, zm);
                alpha2 = overlap(pp, pm, zpm);
                ap = pow(pp, 2);
                am = (alpha1 - alpha2) / PI;
            }
        }

        else{
            casenum = 25;
            if (zpm >= pm + pp){
                ap = pow(pp, 2);
                am = pow(pm, 2);
            }
            else if ((pp - pm < zpm) & (zpm < pp + pm)){
                casenum = 26;
                alpha = overlap(pp, pm, zpm);
                ap = pow(pp, 2);
                am = pow(pm, 2) - alpha / PI;
            }
            else{
                casenum = 27;
                ap = pow(pp, 2);
                am = 0;
            }
        }
    }
    
    static double a[2];
    a[0] = ap;
    a[1] = am;
    return a;
}

void flux(double * flux, double * zp, double * zm, double * zpm, double pp, double pm, double c1, double c2, double c3, double c4, int m){
    
    //FILE *fp;
    //fp = fopen("./debug_flux.txt", "a");
    
    int i;
    double * a;
    double ap;
    double am;
    double Ip;
    double Im;
    double c0 = 1. - c1 - c2 - c3 - c4;
    double Sigma = c0 / 4.0;
    Sigma += c1 / 5.0;
    Sigma += c2 / 6.0;
    Sigma += c3 / 7.0;
    Sigma += c4 / 8.0;
    
    for (i = 0; i < m; i++){
        a = area(zp[i], zm[i], zpm[i], pp, pm);
        ap = *(a);
        am = *(a+1);
        Ip = intensity(zp[i], pp, c1, c2, c3, c4);
        Im = intensity(zm[i], pm, c1, c2, c3, c4);
        flux[i] = - (ap * Ip + am * Im) / (4. * Sigma);
    }
    
    
}