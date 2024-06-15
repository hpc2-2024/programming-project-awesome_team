#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <stdbool.h>
#include "utils.h"
#include "mg_solver.h"

/*! The function f of the exercise sheet*/
double fun(double x, double y){
    return sin(y*M_PI)*sin(x*M_PI)*2.0*M_PI*M_PI;
}

double fun_solution(double x, double y){
    return sin(x*M_PI)*sin(y*M_PI);
}

void init_b(double b[],int N){
    double h = 1.0/(N+1);
    double h2 = h * h;
    // inner points of x_0,b
    for (int i = 1;i<N+1;i++) {
        for (int j = 1;j<N+1;j++){
            b[(N+2)*i+j]=fun(i*h,j*h)*h2; // TB: it is better to shift the h on the rhs, so you do not have a matirx that scale with h
        }
    }
}

int main (int argc, char** argv){
    // Variables init
    int N=129;
    int levels=4;
    int v=3;
    
    if (argc>3){
        N = atoi(argv[1]);
        levels = atoi(argv[2]);
        v = atoi(argv[3]);
        int k = pow(2,levels-1);
        if ((N-k-1)%k!=0){
            printf("(N - 2^(levels-1) -1)/2^(levels-1) has to be an integer (since this is the number of points in the coarsest grid)");
            exit(0);
        }
    }
    else if (argc>2){ // optional gridsize N and number of levels 
        N = atoi(argv[1]);
        levels = atoi(argv[2]);
        int k = pow(2,levels-1);
        if ((N-k-1)%k!=0){
            printf("(N - 2^(levels-1) -1)/2^(levels-1) has to be an integer (since this is the number of points in the coarsest grid)");
            exit(0);
        }
    }

    int vec_ghost = (N+2)*(N+2); //vector size of smallest grid with ghost layer

    // init u (the vector we are solving for)
    double** u;
    u = (double**)malloc(levels*sizeof(double*));
    int Nlevel = N;
    for (int i=levels-1;i>=0;i--){
        u[i]=(double*)malloc( pow( Nlevel+2 , 2) * sizeof(double) );
        null_vec(u[i], pow( Nlevel+2 , 2));
        Nlevel=dim_coarser(Nlevel);
    }

    rand_vec(u[levels-1],N);

    // init right hand side
    double **f;
    f = (double**)malloc(levels*sizeof(double*));
    Nlevel = N;
    for (int i=levels-1;i>=0;i--){
        f[i]=(double*)malloc( pow( Nlevel+2 , 2) * sizeof(double) );
        null_vec(f[i], pow( Nlevel+2 , 2));
        Nlevel=dim_coarser(Nlevel);
    }
    init_b(f[levels-1],N);

    // printf("Right-Hand Side at the Start:\n");
    // print_matrix(N+2, f[levels-1]);

    mg_solve(u,f,N,levels,v);

    //Output
    printf("\n");
    printf("Grid size, N = %d\n",N);
    printf("Number of grids, levels = %d\n",levels);
    printf("Number of smoothing iterations, v = %d\n",v);


    // Free memory
    for (int i=0;i<levels;i++){
        free(u[i]);
    }
    free(u);
    return 0;
}
