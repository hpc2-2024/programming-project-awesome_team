#ifndef utils
#define utils
#endif

#include <stdio.h>
#include <stdlib.h>


// smoothing of u
void smooth(double u[], double f[], int N, int v){

}

void restriction(double r[],int N_r, double f[], int N_f){

}

void prolongation(double u_small[],int N1, double u[], int N2){

}

void exact_solve(double u[],double f[], int N){

}

void v_cycle(double** u, double **f, int N, int levels){
    int vec_size_start = (N+2)*(N+2);
    int vec_size;
    int v = 2;//number of smoothing iterations
    
    for (int l=levels-1;l>=1;l--){
        int Nlevel = N/pow(2,levels-1-l);
        vec_size = (2+Nlevel)*(Nlevel+2);

        double *r;
        r = (double*)malloc((vec_size)*sizeof(double));
        // Smoothing
        smooth(u[l],f[l],Nlevel,v);

        // residual
        mfMult(Nlevel,r,u[l]);
        axpy(r, -1, r, f[l],vec_size);

        // restriction
        int N_f = 10000000; // TODOOOOOOOOOOOO
        restriction(r,Nlevel, f[l-1],N_f);

        free(r);
    }

    // exact solve
    exact_solve(u[0],f[0],N/pow(2,levels-1));

    for (int l=1;l<levels;l++){
        int Nlevel = N/pow(2,levels-1-l);
        int Nlevel2 = N/pow(2,levels-1-(l-1));
        // prolongate
        prolongation(u[l-1],Nlevel,u[l],Nlevel2);
        // Smoothing
        smooth(u[l],f[l],N,v);

    }
}



void mg_solve(double** u, double **f, int N, int levels){
    
    // setup
    int it_max = 1000;
    int iterations = 0;

    double err;
    double epsilon=0.0001;

    int vec_size = (N+2)*(N+2);

    double* r;
    r=(double*)malloc(vec_size*sizeof(double));
    null_vec(r,vec_size);

    
    do {
        iterations += 1;

        // calculate new u with v cycle
        v_cycle(u,f,N,levels);

        // calculate new residual
        mfMult(N,u[levels-1],r);
        axpy(r,-1,r,f[levels-1],vec_size);

        err = norm(r,vec_size);
        printf("Error: %f\n",err);
        if (iterations>it_max){
            printf("Multigrid solve stopped after %d iterations without convergence.", iterations);
            break;
        }
    } while ( err >epsilon);
    
    free(r);
    // solution:
}

