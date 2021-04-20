#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265358979323846

void solve_kepler(double * res, double t, double n, double t0, double e, double a){
    
    double M;
    double tol;
    double E;
    double err;
    double f;
    double r;
    
    M = n * (t - t0);
    if (e <= 1e-5){
        f = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(M / 2.));
        r = a * (1 - e * cos(M));
        res[0] = f;
        res[1] = r;
    }
    
    tol = 1e-7;
    E = M;
    err = e * sin(M);
    while (err > tol){
        err = -(E - e * sin(E) - M) / (1 - e * cos(E));
        E += err;
    }
    
    f = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2.));
    r = a * (1 - e * cos(E));
    
    res[0] = f;
    res[1] = r;
}

double find_transit(double e, double w, double P, double t0, double n){
    
    double arg;
    double E;
    double M;
    
    if (e < 1e-5){
        return P * (PI / 2. - w) / (2 * PI);
    }
    else{
        arg = (1 - (1 - pow(e, 2)) / (1 + e * cos(PI / 2. - w))) / e;
        if (arg > 1){
            arg = 1;
        }
        E = acos(arg);
    }
    M = E - e * sin(E);
    return M / n + t0;
}

void find_xy(double * xy, double t, double n, double t0, double e, double a, double mprim, double msec, double w, double omega, double i){
    
    double mrprim;
    double mrsec;
    double x;
    double y;
    double f;
    double r;
    
    double res[2];
    solve_kepler(res, t, n, t0, e, a);
    f = *(res);
    r = *(res + 1);
    
    mrprim = msec / (mprim + msec);
    mrsec = mprim / (mprim + msec);
    
    x = -r * (cos(omega) * cos(w + f) - sin(omega) * sin(w + f) * cos(i));
    y = -r * (cos(omega) * sin(w + f) * cos(i) + sin(omega) * cos(w + f));
    xy[0] = -mrprim * x;
    xy[1] = -mrprim * y;
    xy[2] = mrsec * x;
    xy[3] = mrsec * y;
}

int solve_kepler_array(double * r ,double * f, double * t, double n, double t0, double e, double a, int m){
    
    int i;
    double res[2];
    for (i=0; i < m; i++){
        solve_kepler(res, t[i], n, t0, e, a);
        f[i] = *(res);
        r[i] = *(res + 1);
    }
    return 0;
}

