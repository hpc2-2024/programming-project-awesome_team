/*
Compile code with: 
gcc -fopenmp ./mg.c -o mg -lm
Execute with e.g.:
./mg 2 57 3 5
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <stdbool.h>
#include "src_mg/utils.h"
#include "src_mg/mg_solver.h"

void print_usage() {
    printf("Usage: ./mg <dimension> <gridsize N> <levels> <smoothing steps>\n");
    printf("Example: ./mg 2 19 2 2\n");
    printf("Note: The dimension must be 1 or 2.\n");
    printf("\nOptional flags: \n");
    printf("-fopenmp: use this flag for shared memory parallelization\n");
    printf("    \"export OMP_NUM_THREADS=...\" before using, set the number of threads on your computer\n" );
    printf("-time: measures the runtime\n");
    printf("-avg_time: runs the multigrid solver 10 times and prints the average runtime\n");
    printf("-wcycle: with this flag the multigrid method uses wcycle instead of vcycle\n");
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
    int use_wcycle = 0;
    int use_stencil9 = 0;
    
    // Check for -time and -avg_time flags
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-time") == 0) {
            measure_time = 1;
            for (int j = i; j < argc - 1; j++) {
                argv[j] = argv[j + 1];
            }
            argc--;  // reduce argument count
            i--;  // adjust index to recheck the current position
        } else if (strcmp(argv[i], "-avg_time") == 0) {
            measure_avg_time = 1;
            for (int j = i; j < argc - 1; j++) {
                argv[j] = argv[j + 1];
            }
            argc--;  // reduce argument count
            i--;  // adjust index to recheck the current position
        }
        if (strcmp(argv[i], "-wcycle") == 0) {
            use_wcycle = 1;
            for (int j = i; j < argc - 1; j++) {
                argv[j] = argv[j + 1];
            }
            argc--;  // reduce argument count
            i--;  // adjust index to recheck the current position
        }
        else if (strcmp(argv[i], "-fcycle") == 0) {
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

    if (argc != 5) {
        print_usage();
        return 1;
    }

    int arg_index = 1;
    int dimension = atoi(argv[arg_index++]);
    int N = atoi(argv[arg_index++]);
    int levels = atoi(argv[arg_index++]);
    int v = atoi(argv[arg_index++]);

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

    if (measure_avg_time) {
        // Measure average time for running mg_solve 10 times
        double total_time = 0;
        for (int i = 0; i < 10; i++) {
            clock_t start_time = clock();

            mg_solve(u, f, N, levels, v, dimension, use_wcycle, fcycle, use_stencil9, 0);

            clock_t end_time = clock();
            total_time += (double)(end_time - start_time) / (10 * CLOCKS_PER_SEC);
            rand_vec(u[levels-1], N, dimension);
            init_b(f[levels-1], N, dimension);
        }
        double avg_time_taken = total_time / 10;
        printf("Average time taken by mg_solve (10 runs): %f seconds\n", avg_time_taken);
    } 
    else if (measure_time) {
        // Measure time for a single run of mg_solve
        clock_t start_time = clock();

        mg_solve(u, f, N, levels, v, dimension, use_wcycle, fcycle, use_stencil9, 1);

        clock_t end_time = clock();
        double time_taken = (double)(end_time - start_time) / (10* CLOCKS_PER_SEC); // Clocks_per_sec should not be multiplied by 10, but for my computer it does for some reason
        printf("Time taken by mg_solve: %f seconds\n", time_taken);
    } 
    else {
        mg_solve(u, f, N, levels, v, dimension, use_wcycle, fcycle, use_stencil9, 1);
    }

    //Output
    printf("\n");
    printf("Grid size, N = %d\n",N);
    printf("Number of grids, levels = %d\n",levels);
    printf("Number of smoothing iterations, v = %d\n",v);

    free_multigrid(u,levels);
    free_multigrid(f,levels);

    return 0;

}
