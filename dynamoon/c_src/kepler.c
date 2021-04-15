#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265358979323846

double * solve_kepler(double t, double n, double t0, double e, double a){
    
    double M;
    double tol;
    double E;
    double err;
    static double ret[2];
    double f;
    double r;
    
    M = n * (t - t0);
    if (e <= 1e-5){
        f = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(M / 2));
        r = a * (1 - e * cos(M));
        ret[0] = f;
        ret[1] = r;
        return ret;
    }
    
    tol = 1e-7;
    E = M;
    err = e * sin(M);
    while (err > tol){
        err = (E - e * sin(E) - M) / (1 - e * cos(E));
        E = E - err;
        err = fabs(err);
    }
    
    f = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));
    r = a * (1 - e * cos(E));
    
    ret[0] = f;
    ret[1] = r;
    return ret;
}

double ** solve_kepler_array_dep(double* t, double n, double t0, double e, double a, int m){
    
    int i;
    double ** fr = calloc(2, sizeof(double*));
    double * r = calloc(m, sizeof(double));
    double * f = calloc(m, sizeof(double));
    
    for (i = 0; i < m; i++ ){
        double * res = solve_kepler(t[i], n, t0, e, a);
        f[i] = *(res);
        r[i] = *(res + 1);
    }
    fr[0] = f;
    fr[1] = r;
    return fr; 
}

int solve_kepler_array(double * r ,double * f, double * t, double n, double t0, double e, double a, int m){
    
    int i;
    for (i=0; i < m; i++){
        double * res = solve_kepler(t[i], n, t0, e, a);
        f[i] = *(res);
        r[i] = *(res + 1);
    }
    return 0;
}

//int main(){
//    double ** res;
//    double t[10] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
//    res = solve_kepler_array(t, 10, 0, 0.2, 1);
//    
//    int m = sizeof(t) / sizeof(double);
//    int i;
//    for (i = 0; i < m; i++ ){
//        printf("%f", res[0][i]);
//        printf(" %f\n", res[1][i]);
//    }
//   return 0;
//}

