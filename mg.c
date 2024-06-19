/*
Compile code with: 
gcc -fopenmp ./mg.c -o mg -lm
Execute with e.g.:
./mg 2 57 3 5
*/
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <stdbool.h>
#include "src/utils.h"
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
            b[(N+2)*i+j]=fun(i*h,j*h); 
        }
    }
}

void init_b_1d(double b[], int N){
    double h = 1.0/(N+1);
    for (int i = 1;i<N+1;i++) {
        b[i]=1.0;  
    }
}

void mg_1dim(int N, int levels, int v){
    int vec_ghost = (N+2); //vector size of smallest grid with ghost layer

    //init u
    double** u = (double**)malloc(levels*sizeof(double*));
    int Nlevel = N;
    for (int i=levels-1;i>=0;i--){
        u[i]=(double*)malloc( (Nlevel+2) * sizeof(double) );
        null_vec(u[i], Nlevel+2);
        Nlevel=dim_coarser(Nlevel);
    }
    rand_vec_1d(u[levels-1],N);
    

    //init f
    double** f = (double**)malloc(levels*sizeof(double*));
    Nlevel = N;
    for (int i=levels-1;i>=0;i--){
        f[i]=(double*)malloc( (Nlevel+2) * sizeof(double) );
        null_vec(f[i], Nlevel+2);
        Nlevel=dim_coarser(Nlevel);
    }
    init_b_1d(f[levels-1],N);
    
    mg_solve(u,f,N,levels,v,1);

    //Speicherfregeben
    for (int i=0;i<levels;i++){
        free(u[i]);
    }
    free(u);

}

void mg_2dim(int N, int levels, int v){
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

    
    mg_solve(u,f,N,levels,v,2);



    //Output
    printf("\n");
    printf("Grid size, N = %d\n",N);
    printf("Number of grids, levels = %d\n",levels);
    printf("Number of smoothing iterations, v = %d\n",v);

    // Free memory
    for (int i=0;i<levels;i++){
        free(u[i]);
    }   
    for (int i=0;i<levels;i++){
        free(f[i]);
    }
    free(f);
}

int main (int argc, char** argv){
    // Variables init
    int N=23;
    int levels=2;
    int v=3;
    int dimension = 2;
    
    if (argc>4){
        dimension = atoi(argv[1]);
        N = atoi(argv[2]);
        levels = atoi(argv[3]);
        v = atoi(argv[4]);
        int k = pow(2,levels-1);
        if ((N-(k-1))%k!=0){
            printf("(N - 2^(levels-1) -1)/2^(levels-1) has to be an integer (since this is the number of points in the coarsest grid)");
            exit(0);
        }
        if (dimension!=1 && dimension!=2){
            printf("Execute file with the following params: dimension - gridsize N - levels - smoothing steps \nE.g. 2 19 2 2\n");
            printf("The dimension has to be 1 or 2\n");
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
    else {
        printf("execute file with the following params: dimension - gridsize N - levels - smoothing steps \nE.g. 2 19 2 2\n");
    }

    if (dimension == 2) {
        mg_2dim(N,levels,v);
    } 
    else {
        printf("1 dim mg\n");
        mg_1dim(N,levels,v);
    }

    return 0;

}
