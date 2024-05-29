//
//  main.cpp
//  3D interg
//
//  Created by John Chou on 12/11/22.
//

#include <iostream>
#include <math.h>
#include <cmath>
//#include <mpi.h>
#include "3D_integration.hpp"
#define EDGE_POINT 11
using namespace std;

double func(double x,double y,double z){
    return 1;
}

double fun_value2(double r){
    if(r<=1.0/3.0){
        return (-r*sin(1.5*M_PI*r*r));
    }
    else if(fabs(r)<1.0/3.0){
        return (fabs(sin(2*M_PI*r)));
    }else if(1.0/3.0<=r){
        return (2*r-1+sin(3*M_PI*r));
    }else{
        return 0;
    }
}

double fun_value1(double x,double y,double z){
    if(x<=0.5*cos(M_PI*y)){
        return (2.0+fun_value2(x-1.0*tan(sqrt(0.5*M_PI))*y));
    }
    else if(x>0.5*cos(M_PI*y)){
        return (2.0+fun_value2(x+1.0/tan(sqrt(0.5*M_PI))*y)+cos(2*M_PI*y));
    }
    else{
        return 0;
    }
}


double fun_value(double x,double y,double z){
    return sin(sqrt(x*x+y*y+z*z))*fun_value1(x,y,z);
}

double fun_center_x(double x,double y,double z){
    return x;
}


int main(int argc, char * argv[]) {
    
    double cubed_sphere_x[EDGE_POINT][EDGE_POINT][EDGE_POINT];
    double cubed_sphere_y[EDGE_POINT][EDGE_POINT][EDGE_POINT];
    double cubed_sphere_z[EDGE_POINT][EDGE_POINT][EDGE_POINT];
    
    double r1=1;
    double r2=10;
    
    for (int k=0; k<EDGE_POINT;k++){
        double r=r1+(r2-r1)*(k/(EDGE_POINT-1.0));
        for (int i=0;i<EDGE_POINT;i++){
            float a=-M_PI/4+M_PI/2*(i/(EDGE_POINT-1.0));
            for (int j=0;j<EDGE_POINT;j++){
                double b=-M_PI/4+M_PI/2*(j/(EDGE_POINT- .0));
                double x=r/sqrt((1+tan(a)*tan(a))+(tan(b)*tan(b)));
                cubed_sphere_x[i][j][k]=x;
                cubed_sphere_y[i][j][k]=x*tan(a);
                cubed_sphere_z[i][j][k]=x*tan(b);}}
    }


    int i=0;
    int j=0;
    int k=0;
    
    double v000[3]={cubed_sphere_x[i][j][k],cubed_sphere_y[i][j][k],cubed_sphere_z[i][j][k]};
    double v100[3]={cubed_sphere_x[i+1][j][k],cubed_sphere_y[i+1][j][k],cubed_sphere_z[i+1][j][k]};
    double v010[3]={cubed_sphere_x[i][j+1][k],cubed_sphere_y[i][j+1][k],cubed_sphere_z[i][j+1][k]};
    double v001[3]={cubed_sphere_x[i][j][k+1],cubed_sphere_y[i][j][k+1],cubed_sphere_z[i][j][k+1]};
    double v110[3]={cubed_sphere_x[i+1][j+1][k],cubed_sphere_y[i+1][j+1][k],cubed_sphere_z[i+1][j+1][k]};
    double v101[3]={cubed_sphere_x[i+1][j][k+1],cubed_sphere_y[i+1][j][k+1],cubed_sphere_z[i+1][j][k+1]};
    double v011[3]={cubed_sphere_x[i][j+1][k+1],cubed_sphere_y[i][j+1][k+1],cubed_sphere_z[i][j+1][k+1]};
    double v111[3]={cubed_sphere_x[i+1][j+1][k+1],cubed_sphere_y[i+1][j+1][k+1],cubed_sphere_z[i+1][j+1][k+1]};
    

    cout<<"x: "<<cubed_sphere_x[i][j][k]<<" "<<cubed_sphere_x[i+1][j][k]<<" "<<cubed_sphere_x[i][j+1][k]<<" "<<cubed_sphere_x[i][j][k+1]<<" "<<cubed_sphere_x[i+1][j+1][k]<<" "<<cubed_sphere_x[i+1][j][k+1]<<" "<<cubed_sphere_x[i][j+1][k+1]<<" "<<cubed_sphere_x[i+1][j+1][k+1]<<endl;
    
    cout<<"y: "<<cubed_sphere_y[i][j][k]<<" "<<cubed_sphere_y[i+1][j][k]<<" "<<cubed_sphere_y[i][j+1][k]<<" "<<cubed_sphere_y[i][j][k+1]<<" "<<cubed_sphere_y[i+1][j+1][k]<<" "<<cubed_sphere_y[i+1][j][k+1]<<" "<<cubed_sphere_y[i][j+1][k+1]<<" "<<cubed_sphere_y[i+1][j+1][k+1]<<endl;
    
    cout<<"z: "<<cubed_sphere_z[i][j][k]<<" "<<cubed_sphere_z[i+1][j][k]<<" "<<cubed_sphere_z[i][j+1][k]<<" "<<cubed_sphere_z[i][j][k+1]<<" "<<cubed_sphere_z[i+1][j+1][k]<<" "<<cubed_sphere_z[i+1][j][k+1]<<" "<<cubed_sphere_z[i][j+1][k+1]<<" "<<cubed_sphere_z[i+1][j+1][k+1]<<endl;
    //double a=Trans.func(0.5, 0.5, 0.5);
    //double b=Trans.det(0.5, 0.5, 0.5);
    //double b=Trans.func(0.5, 0.5, 0.505);
    //double c=Trans.output_coorn(0.5, 0.5, 0.5)[0];
    
    double center_x=AdapSimp3DTrans(5,0.01,10,fun_center_x,v000,v100,v010,v001,v110,v101,v011,v111);
    double volumn=AdapSimp3DTrans(5,0.01,10,func,v000,v100,v010,v001,v110,v101,v011,v111);
    //double AdapSimp3DTrans(int n_intervals,double tol,int maxLevels,double (*func)(double,double,double), double v000[3],double v100[3],double v010[3],double v001[3],double v110[3],double v101[3],double v011[3],double v111[3]);

    
    cout<<"center_x: "<<center_x/volumn<<endl;
    return 0;
}
