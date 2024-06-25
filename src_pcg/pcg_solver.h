#ifndef PCG_SOLVER
#define PCG_SOLVER

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include "../src_mg/poisson_mat_vek.h"
#include "../src_mg/utils.h"
#include "preconditioner.h"

/*! PCG method */
void pcg_solve(int N, double x[],  double b[],  int preconditioner, double epsilon, int debug){
    double N2 = (N+2)*(N+2);
    int vec_size_ghost = (N+2)*(N+2);
    
    double* Ap=(double *)malloc(vec_size_ghost*sizeof(double));
    null_vec(Ap,vec_size_ghost);

    double* r=(double *)malloc(vec_size_ghost*sizeof(double));
    null_vec(r,vec_size_ghost);

    double* z=(double *)malloc(vec_size_ghost*sizeof(double));

    double* p=(double *)malloc(vec_size_ghost*sizeof(double));
    null_vec(p,vec_size_ghost);

    double** a = (double**)malloc(N*N*sizeof(double*));
    for (int i=0; i<N*N; i++){
        a[i] = (double*)malloc(5*sizeof(double));
    }
    

    // Pre loop calculations ( calculating residuum )
    int use_stencil9 = 0;
    poisson_mat_vek(2,N,x, r, use_stencil9);
    axpy(r,-1,b,r,(N+2)*(N+2)); // r = Ax -b (together with last line)
    
    if (preconditioner==0){
        z=r;
    }
    init_preconditioner(a,r,z,N,preconditioner, 0);

    axpy(p,-1,z,0,(N+2)*(N+2)); // p = -r (first conjugated gradient direction)

    double old_r_dot, new_r_dot;
    double alpha, beta;
    double err_k;
    old_r_dot = dot(r,z,N2); // dot product for loop

    int number_of_iterations =0;
    int it_max = 5000;
    do
    {
        number_of_iterations+=1;

        poisson_mat_vek(2, N, p, Ap, use_stencil9); // Ap
        alpha=old_r_dot/dot(p,Ap,N2); // rz/pAp

        // update x,r
        axpy(x,alpha,p,x,N2); //x = x + alpha*p
        axpy(r,alpha,Ap,r,N2); // r=r + alpha*Ap

        
        // precondition r: z = M^-1*r
        apply_precon(a, r, z, N, preconditioner, 0);

        //update p
        new_r_dot = dot(r,z,N2);    // TB :need dot(r,z,N2) instaed of dot(r,r,N2)! 
        beta = new_r_dot/old_r_dot;
        axpby(p,-1,z,beta,p,(N+2)*(N+2));

        old_r_dot = new_r_dot;

        //break criteria
        if (preconditioner==0){
           err_k=sqrt(new_r_dot);
        }
        else {
            err_k=sqrt(dot(z,z,N2));
        }
        if (debug==1) {
            printf("residual it %d: %f \n",number_of_iterations,err_k); // TB: sometimes printing the residual is useful for the debug, you can monitor the convergence (can be commented)
        }

    } while (err_k >= epsilon && number_of_iterations <= it_max);   // TB: Better add a maximum number for the cg iteraions 

    //vec_print(N,x,"vector x"); // TB: I commented this print, use it only when you need it 
    printf("Number of iterations: %d\n",number_of_iterations);

    free(p);
    free(r);
    free(Ap);
    if (preconditioner==1){
        free(z);
    }
}

#endif