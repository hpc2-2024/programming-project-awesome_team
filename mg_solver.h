#include <utils.h>
#include <stdio.h>
#include <stdlib.h>

// smoothing of u
void smooth(int v, double b[], double u[], int N){

}

void mfMult_level(double u[], double r[], int N, int l){

}

void axpy_level(double sol[], double a, double x[], double y[], int arrSize){

}

void v_cycle(double u[], double b[],double r[], int N, int p){
    int vec_size = (N+2)*(N+2);
    int v = 2;//number of smoothing iterations
    
    for (int l=p;l>0;l--){
        // Smoothing
        smooth(v, b[],u[],N);

        // residual
        mfMult_level(u,r,N,l);
        axpy_level(r, -1, r, b,vec_size);

        // restriction
        // simple restriction -> nothing to do here

    }
    // exact solve
    for (int l=p;l>0;l--){
        // prolongate
        // Smoothing

    }
}

void mg_solve(double u[],double b[], double r[], double epsilon, int N, int p){
    
    // setup
    int it_max = 1000;
    int iterations = 0;

    int vec_size = (N+2)*(N+2);

    do {
        iterations += 1;

        // calculate new u with v cycle
        v_cycle(u,b,r,N,p);

        // calculate new residual
        mfMult(N,u,r);
        axpy(r,-1,r,b,vec_size);

        double err; 
        err = norm(r);
        if (iterations>it_max){
            printf("Multigrid solve stopped after %d iterations without convergence.", iterations);
            break;
        }
    } while ( err >epsilon);

    // solution:
}