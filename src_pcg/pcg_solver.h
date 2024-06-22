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
void pcg_solve(double a[][5], int N, double x[], double r[], double b[], double z[], double p[], double Ap[], int preconditioner, double epsilon, int debug){
    double N2 = (N+2)*(N+2);

    // Pre loop calculations ( calculating residuum )
    mfMult(N,x,r); // r = Ax
    axpy(r,-1,b,r,(N+2)*(N+2)); // r = Ax -b (together with last line)
    
    if (preconditioner==0){
        z=r;
    }
    init_preconditioner(a,r,z,N,preconditioner);

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

        mfMult(N,p,Ap); // Ap
        alpha=old_r_dot/dot(p,Ap,N2); // rz/pAp

        // update x,r
        axpy(x,alpha,p,x,N2); //x = x + alpha*p
        axpy(r,alpha,Ap,r,N2); // r=r + alpha*Ap

        
        // precondition r: z = M^-1*r
        apply_precon(a,r,z,N,preconditioner);

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

}

#endif