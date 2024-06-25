/*
Compile code with: 
gcc -fopenmp ./mg.c -o mg -lm
Execute with e.g.:
./mg 2 57 3 5 10
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <stdbool.h>
#include "src_mg/utils.h"
#include "src_mg/mg_solver.h"

void print_usage() {
    printf("Usage: ./mg <dimension> <gridsize N> <levels> <smoothing steps> <iter>\n");
    printf("Example: ./mg 2 19 2 2 10\n");
    printf("Note: The dimension must be 1 or 2.\n");
    printf("\nOptional flags: \n");
    printf("-fopenmp: use this flag for shared memory parallelization\n");
    printf("    \"export OMP_NUM_THREADS=...\" before using, set the number of threads on your computer\n" );
    printf("-fcycle: with this flag the multigrid method uses fcycle instead of vcycle\n");
    printf("-stencil9: uses the 9-point stencil, only works if dimension=2\n\n");
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
    // optional flags
    int measure_time = 0;
    int measure_avg_time = 0;
    int fcycle = 0;
    int use_stencil9 = 0;

    int dimension = 2;
    int N = 15;
    int levels = 2;
    int v = 2; 
    int iter = 1;
    
    // Check for -time and -avg_time flags
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-fcycle") == 0) {
            fcycle = 1;
            for (int j = i; j < argc - 1; j++) {
                argv[j] = argv[j + 1];
            }
            argc--;  // reduce argument count
            i--;  // adjust index to recheck the current position
        }
        if (strcmp(argv[i], "-stencil9") == 0) {
            use_stencil9 = 1;
            for (int j = i; j < argc - 1; j++) {
                argv[j] = argv[j + 1];
            }
            argc--;  // reduce argument count
            i--;  // adjust index to recheck the current position
        }
    }

    // deactivate for debugging  
    if (argc != 6) {
        print_usage();
        return 1;
    }

    int arg_index = 1;
    dimension = atoi(argv[arg_index++]);
    N = atoi(argv[arg_index++]);
    levels = atoi(argv[arg_index++]);
    v = atoi(argv[arg_index++]);
    iter = atoi(argv[arg_index++]); // Measure average time for running mg_solve iter times

    printf("Grid size, N = %d\n",N);
    printf("Number of grids, levels = %d\n",levels);
    printf("Number of smoothing iterations, v = %d\n",v);
    printf("Number of runs = %d\n", iter);
    
    int smoother = 0;
    measure_avg_time = 1;

    if (dimension != 1 && dimension != 2) {
        print_usage();
        printf("\n Input was dimension of %d \n", dimension);
        return 1;
    }
    
    if (use_stencil9 == 1 && dimension!=2) {
        print_usage();
        printf("\n Error: Flag for using 9 point stencil was enable but dimension was not 2 \n");
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

    double total_time = 0;
    for (int i = 0; i < iter; i++) {
        clock_t start_time = clock();

        mg_solve(u, f, N, levels, v, dimension, fcycle, use_stencil9, 0, smoother);

        clock_t end_time = clock();
        total_time += (double)(end_time - start_time) / (10 * CLOCKS_PER_SEC);
        rand_vec(u[levels-1], N, dimension);
        init_b(f[levels-1], N, dimension);
    }
    double avg_time_taken = total_time / iter;
    
    //Output
    printf("\n");

    printf("Average time taken by mg_solve for Jacobi (%d runs): %f seconds\n", iter, avg_time_taken);

    rand_vec(u[levels-1], N, dimension);
    init_b(f[levels-1], N, dimension);
    
    total_time = 0;
    for (int i = 0; i < iter; i++) {
        clock_t start_time = clock();

        mg_solve(u, f, N, levels, v, dimension, fcycle, use_stencil9, 0, smoother);

        clock_t end_time = clock();
        total_time += (double)(end_time - start_time) / (10 * CLOCKS_PER_SEC);
        rand_vec(u[levels-1], N, dimension);
        init_b(f[levels-1], N, dimension);
    }
    avg_time_taken = total_time / iter;
    
    //Output
    printf("Average time taken by mg_solve for Gauss-Seidel (%d runs): %f seconds\n", iter, avg_time_taken);

    free_multigrid(u,levels);
    free_multigrid(f,levels);

    return 0;

}
