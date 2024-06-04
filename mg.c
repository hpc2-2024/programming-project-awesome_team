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
    // inner points of x_0,b
    for (int i = 1;i<N+1;i++) {
        for (int j = 1;j<N+1;j++){
            b[(N+2)*i+j]=fun(i*h,j*h)*h*h; // TB: it is better to shift the h on the rhs, so you do not have a matirx that scale with h
        }
    }
}

int main (int argc, char** argv){
    // Variables init
    int N=10;
    int levels=2;

    if (argc>2){ // optional gridsize N and number of levels 
        N = atoi(argv[1]);
        levels = atoi(argv[2]);
        int k = pow(2,levels-1);
        if (N%k!=0){
            printf("N has to be divisible by 2^(levels-1) for the program to run, please chose different N and or levels\n");
            exit(0);
        }
    }

    int vec_ghost = (N+2)*(N+2); //vector size of smallest grid with ghost layer

    // init u (the vector we are solving for)
    double** u;
    u = (double**)malloc(levels*sizeof(double*));
    for (int i=0;i<levels;i++){
        int Nlevel = N/pow(2,levels-1-i)+2;
        u[i]=(double*)malloc( pow( Nlevel , 2) * sizeof(double) );
        null_vec(u[i], pow( Nlevel , 2));
    }
    rand_vec(u[levels-1],N);

    // init right hand side
    double **f;
    f = (double**)malloc(levels*sizeof(double*));
    for (int i=0;i<levels;i++){
        int Nlevel = N/pow(2,levels-1-i);
        f[i]=(double*)malloc( pow( Nlevel+2 , 2) * sizeof(double) );
        null_vec(f[i], pow( Nlevel+2 , 2));
    }
    init_b(f[levels-1],N);

    
    mg_solve(u,f,N,levels);




    //Output
    printf("\n");
    printf("Grid size, N = %d\n",N);
    printf("Number of grids, levels = %d\n",levels);


    //Speicherfregeben
    for (int i=0;i<levels;i++){
        free(u[i]);
    }
    free(u);
}
