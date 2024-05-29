//
//  3D_integration.cpp
//  3D interg
//
//  Created by John Chou on 12/11/22.
//

#include "3D_integration.hpp"


double simpson_1D_1step(Trilinear_Trans Trans,double x,double y,double az,double bz,double tol, int level){
    double sum=(bz-az)*(Trans.func(x,y,az)+Trans.func(x,y,bz)+4*Trans.func(x,y,(az+bz)/2))/6;
    double sum_c=(bz-az)*(Trans.func(x,y,az)+Trans.func(x,y,bz)+2*Trans.func(x,y,(az+bz)/2)+4*Trans.func(x,y,(3*az+bz)/4)+4*Trans.func(x,y,(az+3*bz)/4))/12;
    
    double error=abs(sum_c-sum);
    if(error>abs(tol*sum)&&level>0){
        double result_left=simpson_1D_1step(Trans,x,y,az,(az+bz)/2,tol,level-1);
        double result_right=simpson_1D_1step(Trans,x,y,(az+bz)/2,bz,tol,level-1);
        sum = result_left+result_right;
    }
    
    return sum;
}

double adapative_simpson_1D(Trilinear_Trans Trans,double x,double y,double az,double bz,int n_intervals,double tol,int maxLevels){
    
    double sum=0;
    double dz=(bz-az)/n_intervals;
    double z_curr=az;
    for(int i=0;i<n_intervals;i++){
        double sub_sum;
        sub_sum=simpson_1D_1step(Trans,x,y,z_curr,z_curr+dz,tol,maxLevels);
        sum = sum+sub_sum;
        z_curr=z_curr+dz;
    }
    return sum;
}

double simpson_2D_1step(Trilinear_Trans Trans,double x,double ay,double by,double az,double bz,double tol, int level,int n_intervals,int maxLevels){
    
    double func_s=adapative_simpson_1D(Trans,x,ay,az,bz,n_intervals,tol,maxLevels);
    double func_e=adapative_simpson_1D(Trans,x,by,az,bz,n_intervals,tol,maxLevels);
    double func_m=adapative_simpson_1D(Trans,x,(ay+by)/2,az,bz,n_intervals,tol,maxLevels);
    double func_1q=adapative_simpson_1D(Trans,x,(3*ay+by)/4,az,bz,n_intervals,tol,maxLevels);
    double func_3q=adapative_simpson_1D(Trans,x,(ay+3*by)/4,az,bz,n_intervals,tol,maxLevels);
    
    double sum=(by-ay)*(func_s+func_e+4*func_m)/6;
    double sum_c=(by-ay)*(func_s+func_e+2*func_m+4*func_1q+4*func_3q)/12;
    
    double error=abs(sum_c-sum);
    if(error>abs(tol*sum)&&level>0){
        double result_left=simpson_2D_1step(Trans,x,ay,(ay+by)/2,az,bz,tol,level-1,n_intervals,maxLevels);
        double result_right=simpson_2D_1step(Trans,x,(ay+by)/2,by,az,bz,tol,level-1,n_intervals,maxLevels);
        sum = result_left+result_right;
    }
    
    return sum;
}

double adapative_simpson_2D(Trilinear_Trans Trans,double x,double ay,double by,double az,double bz,int n_intervals,double tol,int maxLevels){
    
    double sum=0;
    double dy=(by-ay)/n_intervals;
    double y_curr=ay;
    for(int i=0;i<n_intervals;i++){
        double sub_sum;
        sub_sum=simpson_2D_1step(Trans,x,y_curr,y_curr+dy,az,bz,tol,maxLevels,n_intervals,maxLevels);
        sum = sum+sub_sum;
        y_curr=y_curr+dy;
    }
    return sum;
}

double simpson_3D_1step(Trilinear_Trans Trans,double ax,double bx,double ay,double by,double az,double bz,double tol, int level,int n_intervals,int maxLevels){
    
    double func_s=adapative_simpson_2D(Trans,ax,ay,by,az,bz,n_intervals,tol,maxLevels);
    double func_e=adapative_simpson_2D(Trans,bx,ay,by,az,bz,n_intervals,tol,maxLevels);
    double func_m=adapative_simpson_2D(Trans,(ax+bx)/2,ay,by,az,bz,n_intervals,tol,maxLevels);
    double func_1q=adapative_simpson_2D(Trans,(3*ax+bx)/4,ay,by,az,bz,n_intervals,tol,maxLevels);
    double func_3q=adapative_simpson_2D(Trans,(ax+3*bx)/4,ay,by,az,bz,n_intervals,tol,maxLevels);
    
    double sum=(bx-ax)*(func_s+func_e+4*func_m)/6;
    double sum_c=(bx-ax)*(func_s+func_e+2*func_m+4*func_1q+4*func_3q)/12;
    
    double error=abs(sum_c-sum);
    if(error>abs(tol*sum)&&level>0){
        double result_left=simpson_3D_1step(Trans,ax,(ax+bx)/2,ay,by,az,bz,tol,level-1,n_intervals,maxLevels);
        double result_right=simpson_3D_1step(Trans,(ax+bx)/2,bx,ay,by,az,bz,tol,level-1,n_intervals,maxLevels);
        sum = result_left+result_right;
    }
    
    return sum;
}

/*double adapative_simpson_3D(double (*func)(double,double,double),double ax,double bx,double ay,double by,double az,double bz,int n_intervals,double tol,int maxLevels=10){
    
    double sum=0;
    double dx=(bx-ax)/n_intervals;
    double x_curr=ax;
    for(int i=0;i<n_intervals;i++){
        double sub_sum;
        sub_sum=simpson_3D_1step(func,x_curr,x_curr+dx,ay,by,az,bz,tol,maxLevels,n_intervals,maxLevels);
        sum = sum+sub_sum;
        x_curr=x_curr+dx;
    }
    return sum;
}*/

double AdapSimp3DTrans(int n_intervals,double tol,int maxLevels,double (*func)(double,double,double), double v000[3],double v100[3],double v010[3],double v001[3],double v110[3],double v101[3],double v011[3],double v111[3]){
    
    Trilinear_Trans Trans(func,v000,v100,v010,v001,v110,v101,v011,v111);
    
    double sum=0;
    double dx=(1.0-0.0)/n_intervals;
    double x_curr=0;
    for(int i=0;i<n_intervals;i++){
        double sub_sum;
        sub_sum=simpson_3D_1step(Trans,x_curr,x_curr+dx,0,1,0,1,tol,maxLevels,n_intervals,maxLevels);
        sum = sum+sub_sum;
        x_curr=x_curr+dx;
    }
    return sum;
}

Trilinear_Trans::Trilinear_Trans(double (*funct)(double,double,double),double v000[3],double v100[3],double v010[3],double v001[3],double v110[3],double v101[3],double v011[3],double v111[3]){
    
    function=funct;
    
    A[0]=v000[0];
    A[1]=v000[1];
    A[2]=v000[2];
    
    B[0]=v100[0]-A[0];
    B[1]=v100[1]-A[1];
    B[2]=v100[2]-A[2];
    
    C[0]=v010[0]-A[0];
    C[1]=v010[1]-A[1];
    C[2]=v010[2]-A[2];
    
    D[0]=v001[0]-A[0];
    D[1]=v001[1]-A[1];
    D[2]=v001[2]-A[2];
    
    E[0]=v110[0]-A[0]-B[0]-C[0];
    E[1]=v110[1]-A[1]-B[1]-C[1];
    E[2]=v110[2]-A[2]-B[2]-C[2];
    
    F[0]=v101[0]-A[0]-B[0]-D[0];
    F[1]=v101[1]-A[1]-B[1]-D[1];
    F[2]=v101[2]-A[2]-B[2]-D[2];
    
    G[0]=v011[0]-A[0]-D[0]-C[0];
    G[1]=v011[1]-A[1]-D[1]-C[1];
    G[2]=v011[2]-A[2]-D[2]-C[2];
    
    H[0]=v111[0]-A[0]-B[0]-C[0]-D[0]-E[0]-F[0]-G[0];
    H[1]=v111[1]-A[1]-B[1]-C[1]-D[1]-E[1]-F[1]-G[1];
    H[2]=v111[2]-A[2]-B[2]-C[2]-D[2]-E[2]-F[2]-G[2];
    
}

Trilinear_Trans::~Trilinear_Trans(){}

double* Trilinear_Trans::output_coorn(double p,double q,double r){
    double *result = new double[3];
    
    for(int i=0;i<3;i++){
        result[i]=A[i]+B[i]*p+C[i]*q+D[i]*r+E[i]*q*p+F[i]*p*r+G[i]*q*r+H[i]*p*q*r;
    }
    
    return result;
}

double Trilinear_Trans::det(double p,double q,double r){
    
    double rp[3],rq[3],rr[3],result=0;
    for(int i=0;i<3;i++){
        rp[i]=B[i]+E[i]*q+F[i]*r+H[i]*q*r;
        rq[i]=C[i]+E[i]*p+G[i]*r+H[i]*p*r;
        rr[i]=D[i]+F[i]*p+G[i]*q+H[i]*p*q;
    }
    /*double Jaco[3][3];
    for(int i=0;i<3;i++){
        Jaco[i][0]=rp[i];
        Jaco[i][1]=rq[i];
        Jaco[i][2]=rr[i];
    }*/
    double rqrr[3];
    rqrr[0]=rq[1]*rr[2]-rq[2]*rr[1];
    rqrr[1]=rr[0]*rq[2]-rq[0]*rr[2];
    rqrr[2]=rq[0]*rr[1]-rr[0]*rq[1];
    
    result=rp[0]*rqrr[0]+rp[1]*rqrr[1]+rp[2]*rqrr[2];
    
    return result;
}
double Trilinear_Trans::func(double p,double q,double r){
    double *coordn = new double[3];
    coordn=output_coorn(p,q,r);
    double result=function(coordn[0],coordn[1],coordn[2]);
    return result*det(p,q,r);
}
