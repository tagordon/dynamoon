#ifndef KEPLER_H_   
#define KEPLER_H_

void solve_kepler(double * res, double t, double n, double t0, double e, double a);

double * find_xyz(double * xy, double t, double n, double t0, double e, double a, double mprim, double msec, double w, double omega, double i);

double find_transit(double e, double w, double P, double t0, double n);

void find_xyz_array(double * xp, double * yp, double * zp, double * xs, double * ys, double * zs, double t, double n, double t0, double e, double a, double mprim, double msec, double w, double omega, double i, int m);

#endif 