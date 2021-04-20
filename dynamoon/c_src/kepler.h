#ifndef KEPLER_H_   
#define KEPLER_H_

void solve_kepler(double * res, double t, double n, double t0, double e, double a);

double * find_xy(double * xy, double t, double n, double t0, double e, double a, double mprim, double msec, double w, double omega, double i);

double find_transit(double e, double w, double P, double t0, double n);

#endif 