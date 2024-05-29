//
//  3D_integration.hpp
//  3D interg
//
//  Created by John Chou on 12/11/22.
//

#ifndef _D_integration_hpp
#define _D_integration_hpp

#include <stdio.h>
#include <math.h>
#include <cmath>

class Trilinear_Trans{
public:
    Trilinear_Trans(double (*func)(double,double,double), double v000[3],double v100[3],double v010[3],double v001[3],double v110[3],double v101[3],double v011[3],double v111[3]);
    ~Trilinear_Trans();
    double* output_coorn(double p,double q,double r);
    double det(double p,double q,double r);
    double func(double p,double q,double r);
private:
    double A[3],B[3],C[3],D[3],E[3],F[3],G[3],H[3];
    double (*function)(double,double,double);
};

double simpson_1D_1step(Trilinear_Trans Trans,double x,double y,double az,double bz,double tol, int level);

double adapative_simpson_1D(Trilinear_Trans Trans,double x,double y,double az,double bz,int n_intervals,double tol,int maxLevels);

double simpson_2D_1step(Trilinear_Trans Trans,double x,double ay,double by,double az,double bz,double tol, int level,int n_intervals,int maxLevels);

double adapative_simpson_2D(Trilinear_Trans Trans,double x,double ay,double by,double az,double bz,int n_intervals,double tol,int maxLevels);

double simpson_3D_1step(Trilinear_Trans Trans,double ax,double bx,double ay,double by,double az,double bz,double tol, int level,int n_intervals,int maxLevels);

//double adapative_simpson_3D(double (*func)(double,double,double),double ax,double bx,double ay,double by,double az,double bz,int n_intervals,double tol,int maxLevels);

double AdapSimp3DTrans(int n_intervals,double tol,int maxLevels,double (*func)(double,double,double), double v000[3],double v100[3],double v010[3],double v001[3],double v110[3],double v101[3],double v011[3],double v111[3]);
#endif /* _D_integration_hpp */
