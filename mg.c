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
#include "src_mg/utils.h"
#include "src_mg/mg_solver.h"

void print_usage() {
    printf("Usage: ./mg <dimension> <gridsize N> <levels> <smoothing steps>\n");
    printf("E.g.: ./mg 2 19 2 2\n");
    printf("The dimension must be 1 or 2.\n");
}

bool is_valid_input(int N, int levels) {
    int k = pow(2, levels - 1);
    return (N - (k - 1)) % k == 0;
}

/*! The function f of the exercise sheet*/
double fun(double x, double y){
    return sin(y*M_PI)*sin(x*M_PI)*2.0*M_PI*M_PI;
}

double fun_solution(double x, double y){
    return sin(x*M_PI)*sin(y*M_PI);
}

void init_b(double b[],int N, int dim){
    double h = 1.0/(N+1);
    double h2 = h * h;
    // inner points of x_0,b
    for (int i = 1;i<N+1;i++) {
        if (dim==2) {
            for (int j = 1;j<N+1;j++){
                b[(N+2)*i+j]=fun(i*h,j*h); 
            }
        }
        else {
            b[i]=1.0;  
        }
    }
}

int main (int argc, char** argv){

    if (argc != 5) {
        print_usage();
        return 1;
    }

    int dimension = atoi(argv[1]);
    int N = atoi(argv[2]);
    int levels = atoi(argv[3]);
    int v = atoi(argv[4]);

    if (dimension != 1 && dimension != 2) {
        print_usage();
        return 1;
    }

    if (!is_valid_input(N, levels)) {
        printf("Error: (N - (2^(levels-1) - 1)) must be divisible by 2^(levels-1).\n");
        printf("E.g.: N = 7, 15, 31, 63, 127, ... can have multiple levels \n");
        printf("This ensures that the number of points in the coarsest grid is an integer.\n");
        return 1;
    }

    double** u = allocate_multigrid(N, levels, dimension);
    double** f = allocate_multigrid(N, levels, dimension);

    rand_vec(u[levels-1], N, dimension);
    init_b(f[levels-1], N, dimension);

    mg_solve(u,f,N,levels,v,dimension);

    //Output
    printf("\n");
    printf("Grid size, N = %d\n",N);
    printf("Number of grids, levels = %d\n",levels);
    printf("Number of smoothing iterations, v = %d\n",v);

    free_multigrid(u,levels);
    free_multigrid(f,levels);

    return 0;

}
