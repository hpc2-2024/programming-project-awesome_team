#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include "utils.h"
#include "ilu.h"

/*! The function f of the exercise sheet*/
double fun(double x, double y){
    return sin(y*M_PI)*sin(x*M_PI)*2.0*M_PI*M_PI;
}

double fun_solution(double x, double y){
    return sin(x*M_PI)*sin(y*M_PI);
}

/*! Implementation of a matrix free multiplication with 5-star stencil*/
void mfMult(int N, double r[], double y[]){
    for (int i=1;i<N+1;i++){
        for (int j=1;j<N+1;j++){
            y[(N+2)*i+j]=4*r[(N+2)*i+j]-r[(N+2)*(i+1)+j] -r[(N+2)*i+j-1]-r[(N+2)*(i-1)+j]-r[(N+2)*i+j+1];
            //y[(N+2)*i+j]=1.0/(h*h)*y[(N+2)*i+j]; // TB: it is better to shift the h on the rhs, so you do not have a matirx that scale with h
        }
    }
}

void apply_precon(double a[][5],double r[],double temp[],double z[],int N){
    forward_solve(a,temp,r,N);
    backward_solve(a,z,temp,N);

    // inv_diag(z,r,N); TB: if you want you can try the simplest preconditioner (inverse of diagonal, can be useful in debug)
}

void pcg_solve(double a[][5], int N, double x[], double r[], double b[], double temp[], double z[], double p[], double Ap[], int preconditioner, double epsilon){
    double N2 = (N+2)*(N+2);

    lapl_matrix(a,N);
    ilu(a,N,0.001,100);
    // Pre loop calculations ( calculating residuum )
    mfMult(N,x,r); // r = Ax
    axpy(r,-1,b,r,(N+2)*(N+2)); // r = Ax -b (together with last line)
    if (preconditioner == 1) {
        //precondition z = M^-1r
        apply_precon(a,r,temp,z,N);
    }
    else {
        z = r; // no preconditioning of the residuum
    }
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

        if (preconditioner==1) {
            // precondition r: z = Mr
            apply_precon(a,r,temp,z,N);
        }
        else {      
            z = r; 
        }           

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
        printf("residual it %d: %f \n",number_of_iterations,err_k); // TB: sometimes printing the residual is useful for the debug, you can monitor the convergence (can be commented)
    } while (err_k >= epsilon && number_of_iterations <= it_max);   // TB: Better add a maximum number for the cg iteraions 

    //vec_print(N,x,"vector x"); // TB: I commented this print, use it only when you need it 
    printf("Number of iterations: %d\n",number_of_iterations);

    
}

int main(int argc, char** argv){
    int preconditioner=0;
    int number_of_iterations = 0;
    int it_Max = 1000; // TB: add a maximum for the cg iterations
    // if you want to use a precondtioner first argument 1
    if (argc>1){
        preconditioner=atoi(argv[1]);
    }
    printf("Preconditioner: %d \n",preconditioner);
    // initilize variables
    int N = 50; // N^2 is the number of inner points in our lattice
    int N2 = (N+2)*(N+2);
    double epsilon = 1e-3;
    double h = 1.0/(N+1);
    // init x0
    double alpha = 0;
    double beta=0;
    double err0, errk;
    double *x,*p,*r,*b,*m, *z,*temp;
    x=(double *)malloc((N+2)*(N+2)*sizeof(double));
    z=(double *)malloc((N+2)*(N+2)*sizeof(double));
    p=(double *)malloc((N+2)*(N+2)*sizeof(double));
    r=(double *)malloc((N+2)*(N+2)*sizeof(double));
    b=(double *)malloc((N+2)*(N+2)*sizeof(double));
    m=(double *)malloc((N+2)*(N+2)*sizeof(double));
    temp=(double *)malloc((N+2)*(N+2)*sizeof(double));

    // Creation of matrix a
    double a[N*N][5];
    

    //print_2dim(N, a,"mat ilu0");

    // fill ghost layer with zeros (and everything else also 0)
    for (int i = 0;i<N+2;i++) {
        for (int j = 0;j<N+2;j++){
            x[(N+2)*i+j]=0;
            p[(N+2)*i+j]=0;
            r[(N+2)*i+j]=0;
            b[(N+2)*i+j]=0;
            m[(N+2)*i+j]=0;
        }

    }

    int seed = 123456;
    // inner points of x_0,b
    for (int i = 1;i<N+1;i++) {
        for (int j = 1;j<N+1;j++){
            // randomly initialize x with values in (0,1)
            srand(seed + i);
            double r = (double)rand() / (double)RAND_MAX;
            x[(N+2)*i+j]=r;

            b[(N+2)*i+j]=fun(i*h,j*h)*h*h; // TB: it is better to shift the h on the rhs, so you do not have a matirx that scale with h
        }

    }

    pcg_solve(a,N,x,r,b,temp,z,p,m,preconditioner,epsilon);

    // Compute the relative absolute difference 
    double abs_diff = 0.0;
    double abs_sol = 0.0;
    for (int i = 1;i<N+1;i++) {
        for (int j = 1;j<N+1;j++){
            double x_sol = fun_solution(i*h,j*h);
            abs_diff += fabs(x[(N+2)*i+j] - x_sol);
            abs_sol += x_sol;
        }
    }
    printf("Relative absolute difference between exact and approx. solution: %f%%\n", 100.0 * abs_diff / abs_sol);

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