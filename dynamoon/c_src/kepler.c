#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265358979323846

// Kepler solver using Newton's method
void solve_kepler(double * res, double t, double n, double t0, double e, double a){
    
    double M;
    double tol;
    double E;
    double err;
    double f;
    double r;
    
    M = n * (t - t0);
    if (e <= 1e-5){
        f = M;
        r = a;
        res[0] = f;
        res[1] = r;
    }
    
    else{
    
        tol = 1e-9;
        E = M;
        err = 1;
        while (err > tol){
            err = -(E - e * sin(E) - M) / (1 - e * cos(E));
            E += err;
        }
    
        f = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2.));
        r = a * (1 - e * cos(E));
    
        res[0] = f;
        res[1] = r;
    }
}

// find the transit time 
double find_transit(double e, double w, double P, double t0, double n){
    
    double arg;
    double E;
    double M;
    
    if (e < 1e-5){
        return P * (PI / 2. - w) / (2 * PI) + t0;
    }
    else{
        arg = (1 - (1 - pow(e, 2)) / (1 + e * cos(PI / 2. - w))) / e;
        
        // for some reason the argument ocassionally exceeds 1 (at least in the python implementation)
        // so this makes sure it doesn't. Should probably figure out what the actual issue is at some 
        // point but this seems fine for now.
        if (arg > 1){
            arg = 1;
        }
        E = acos(arg);
    }
    M = E - e * sin(E);
    return M / n + t0;
}

// find the coordinates of the primary and secondary bodies at time t
void find_xyz(double * xyz, double t, double n, double t0, double e, double a, double mprim, double msec, double w, double omega, double i){
    
    double mrprim;
    double mrsec;
    double x;
    double y;
    double z;
    double f;
    double r;
    
    double res[2];
    solve_kepler(res, t, n, t0, e, a);
    f = *(res);
    r = *(res + 1);
    
    mrprim = msec / (mprim + msec);
    mrsec = mprim / (mprim + msec);
    
    x = -r * (cos(omega) * cos(w + f) - sin(omega) * sin(w + f) * cos(i));
    y = -r * (sin(omega) * cos(w + f) + cos(omega) * sin(w + f) * cos(i));
    z = r * sin(w + f) * sin(i);
    
    xyz[0] = -mrprim * x;
    xyz[1] = -mrprim * y;
    xyz[2] = -mrprim * z;
    xyz[3] = mrsec * x;
    xyz[4] = mrsec * y;
    xyz[5] = mrsec * z;
}

// wrapping the above function to populate arrays of coordinates. I don't think I finished this, and I don't think 
// it's used. Remove if you can. 
void find_xyz_array_dep(double * xp, double * yp, double * zp, double * xs, double * ys, double * zs, double t, double n, double t0, double e, double a, double mprim, double msec, double w, double omega, double i, int m){
    
    int j;
    double xyz[6];
    
    for (j=0; j<m; j++){
        find_xyz(xyz, t, n, t0, e, a, mprim, msec, w, omega, i);
        xp[j] = xyz[0];
        yp[j] = xyz[1];
        zp[j] = xyz[2];
        xs[j] = xyz[3];
        ys[j] = xyz[4];
        zs[j] = xyz[5];
    }
}

// Wraps kepler solver to populate r and f arrays in time 
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

