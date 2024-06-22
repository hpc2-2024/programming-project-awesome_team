#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include "src_mg/poisson_mat_vek.h"
#include "src_mg/utils.h"
#include "src_pcg/preconditioner.h"
#include "src_pcg/pcg_solver.h"

/*! The function f of the exercise sheet*/
double fun(double x, double y){
    return sin(y*M_PI)*sin(x*M_PI)*2.0*M_PI*M_PI;
}

double fun_solution(double x, double y){
    return sin(x*M_PI)*sin(y*M_PI);
}

void init_b(double b[],int N){
    double h = 1.0/(N+1);
    // inner points of x_0,b
    for (int i = 1;i<N+1;i++) {
        for (int j = 1;j<N+1;j++){
            b[(N+2)*i+j]=fun(i*h,j*h)*h*h; // TB: it is better to shift the h on the rhs, so you do not have a matirx that scale with h
        }
    }
}

double delta_rel(double x[], int N){
    double h = 1.0/(N+1);
    double abs_diff = 0.0;
    double abs_sol = 0.0;
    for (int i = 1;i<N+1;i++) {
        for (int j = 1;j<N+1;j++){
            double x_sol = fun_solution(i*h,j*h);
            abs_diff += fabs(x[(N+2)*i+j] - x_sol);
            abs_sol += x_sol;
        }
    }
    return (abs_diff/abs_sol);
}





int main(int argc, char** argv){

    int preconditioner=0;
    int N = 50; // N^2 is the number of inner points in our lattice
    int debug = 1; // option for printing infos of pcg (residual at each iteration)
    // Options when running code
    if (argc>1){ // preconditioner
        preconditioner=atoi(argv[1]);
    }
    if (argc>2){ // gridsize N
        N = atoi(argv[2]);
    }
    if (N>100){
        debug = 0;
    }

    // initilize variables
    double h = 1.0/(N+1);

    int N2 = (N+2)*(N+2);
    int vec_size_ghost = (N+2)*(N+2);

    double epsilon = 1e-3;

    double *x,*p,*r,*b,*m, *z,*temp;

    x=(double *)malloc(vec_size_ghost*sizeof(double));
    null_vec(x,vec_size_ghost);
    rand_vec(x, N, 2); // random start vector x_0 

    z=(double *)malloc(vec_size_ghost*sizeof(double));

    p=(double *)malloc(vec_size_ghost*sizeof(double));
    null_vec(p,vec_size_ghost);

    r=(double *)malloc(vec_size_ghost*sizeof(double));
    null_vec(r,vec_size_ghost);

    b=(double *)malloc(vec_size_ghost*sizeof(double));
    null_vec(b,vec_size_ghost);
    init_b(b,N);

    m=(double *)malloc(vec_size_ghost*sizeof(double));
    null_vec(m,vec_size_ghost);

    // Creation of matrix a
    double a[N*N][5];
    
    // PCG solve
    pcg_solve(a,N,x,r,b,z,p,m,preconditioner,epsilon,debug);

    // Compute the relative absolute difference 
    double rel_dif = delta_rel(x,N); 

    // Output
    printf("\nSetup of cg:\n");
    printf("size of grid: N = %d\n",N);
    printf("Preconditioner: %d \n",preconditioner);
    printf("break condition of cg loop: epsilon= %f\n", epsilon);
    printf("\n");

    printf("Relative absolute difference between exact and approx. solution: %f%%\n", 100.0 * rel_dif);

    // free allocated memory
    free(x);
    free(p);
    free(r);
    free(b);
    free(m);
    if (preconditioner==1){
        free(z);
    }
    return 0;
}