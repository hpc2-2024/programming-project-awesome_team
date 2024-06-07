#ifndef utils
#define utils
#endif

#include <stdio.h>
#include <stdlib.h>

int dim_finer(int N){
    return N*2 + 1;
}

int dim_coarser(int N){
    return (N-1)/2;
}

// smoothing of u
void smooth(double u[], double f[], int N, int v){

}

/*! 
    N*N is the dimension of the fine grid
    (N+2)*(N+2) is the dimension of the vector with ghost layer 
    
    M*M is the dimension of the coarse grid
    (M+2)^2 is the dimension of the vector
*/
void restriction(double *fine_grid, int N, double *coarse_grid, int M){
    M=M+2;
    N=N+2;

    for(int i = 1; i<M-1; i++){
        for(int j = 1; j<M-1; j++){
            int k = 2*i;
            int l = 2*j;

            coarse_grid[i * M + j] = 0.125  * (fine_grid[(k+1) * N + l] + fine_grid[(k-1) * N + l] + fine_grid[k * N + (l-1)] + fine_grid[k * N + (l+1)]) + 0.5 * fine_grid[k * N + l];
        }
    }
}


void prolongation(double u_small[],int N1, double u[], int N2){

}

void exact_solve(double u[],double f[], int N){

}

void v_cycle(double** u, double **f, int N, int levels){
    int vec_size_start = (N+2)*(N+2);
    int vec_size;
    int v = 2;//number of smoothing iterations
    
    int Nlevel= N;
    for (int l=levels-1;l>=1;l--){
        // number of points in the grid (without border)
        vec_size = (2+Nlevel)*(Nlevel+2);

        double *r;
        r = (double*)malloc((vec_size)*sizeof(double));
        // Smoothing
        smooth(u[l],f[l],Nlevel,v);

        // residual
        mfMult(Nlevel,r,u[l]);
        axpy(r, -1, r, f[l],vec_size);

        // restriction
        int N_f = dim_coarser(Nlevel); 
        restriction(r,Nlevel, f[l-1],N_f);
        Nlevel=N_f;

        free(r);
    }

    // exact solve
    exact_solve(u[0],f[0],Nlevel);

    for (int l=1;l<levels;l++){
        int Nlevel_smaller = Nlevel;
        int Nlevel = dim_finer(Nlevel);

        // prolongate
        prolongation(u[l-1],Nlevel,u[l],Nlevel_smaller);
        // Smoothing
        smooth(u[l],f[l],Nlevel,v);

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

