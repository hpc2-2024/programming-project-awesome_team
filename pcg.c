#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include "src_mg/poisson_mat_vek.h"
#include "src_mg/utils.h"
#include "src_pcg/preconditioner.h"
#include "src_pcg/pcg_solver.h"

/*! @file */

void print_usage(){
    printf("Usage: ./pcg <gridsize N> <preconditioner>\n");
    printf("Example: ./pcg 63 4\n");
    printf("Optional preconditioner: \n");
    printf("1: Ilu(0)\n");
    printf("2: Jacobi\n");
    printf("3: Gauß-Seidel\n");
    printf("4: Multigrid\n\n");

    printf("\nOptional flags: \n");
    printf("-time: measures the runtime\n");
}

/*! 
* @brief Solve a linear system using the preconditioned Conjugate Gradients method. 
    @note Usage: ./pcg <gridsize> <preconditioner>
    @note Example: ./pcg 63 4
    @note Note: set <preconditioner> to one of the following values to select different preconditioning methods: 1 = ILU, 2 = Jacobi, 3 = Gauss-Seidel, 4 = Multi-Grid
    @note Optional flags:
    @note   -time : with this flag, perform a single PCG run and output the time.
*/
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

double delta_rel(double x[], int N){
    double h = 1.0/(N+1);
    double abs_diff = 0.0;
    double abs_sol = 0.0;
    for (int i = 1;i<N+1;i++) {
        for (int j = 1;j<N+1;j++){
            double x_sol = fun_solution(i*h,j*h);
            abs_diff += fabs(x[(N+2)*i+j] - x_sol);
            abs_sol += x_sol;
        }
    }
    return (abs_diff/abs_sol);
}

int main(int argc, char** argv){

    int preconditioner=0;
    int N = 50; // N^2 is the number of inner points in our lattice
    int debug = 1; // option for printing infos of pcg (residual at each iteration)
    int measure_time = 0;

    // Check for -time flag
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-time") == 0) {
            measure_time = 1;
            for (int j = i; j < argc - 1; j++) {
                argv[j] = argv[j + 1];
            }
            argc--;  // reduce argument count
            i--;  // adjust index to recheck the current position
        } 
    }

    // Options when running code
    if (argc>=2){ // preconditioner
        N = atoi(argv[1]);
        preconditioner=atoi(argv[2]);
    }
    else {
        print_usage();
        return 0;
    }

    if (N>100){
        debug = 0;
    }

    // initilize variables
    double h = 1.0/(N+1);

    int N2 = (N+2)*(N+2);
    int vec_size_ghost = (N+2)*(N+2);

    double epsilon = 1e-4;

    double *x,*b;

    x=(double *)malloc(vec_size_ghost*sizeof(double));
    null_vec(x,vec_size_ghost);
    rand_vec(x, N, 2); // random start vector x_0 

    b=(double *)malloc(vec_size_ghost*sizeof(double));
    null_vec(b,vec_size_ghost);
    init_b(b,N);

    // PCG solve
    if (measure_time) {
        // Measure time for a single run of mg_solve
        clock_t start_time = clock();

        pcg_solve(N,x,b,preconditioner,epsilon,debug);

        clock_t end_time = clock();
        double time_taken = (double)(end_time - start_time) / (10* CLOCKS_PER_SEC); // Clocks_per_sec should not be multiplied by 10, but for my computer it does for some reason
        printf("Time taken by pcg_solve: %f seconds\n", time_taken);
    } 
    else {
        pcg_solve(N,x,b,preconditioner,epsilon,debug);
    }

    // Compute the relative absolute difference 
    double rel_dif = delta_rel(x,N); 

    // Output
    printf("\nSetup of cg:\n");
    printf("size of grid: N = %d\n",N);
    printf("Preconditioner: %d \n",preconditioner);
    printf("break condition of cg loop: epsilon= %f\n", epsilon);
    printf("\n");

    printf("Relative absolute difference between exact and approx. solution: %f%%\n", 100.0 * rel_dif);

    // free allocated memory
    free(x);
    free(b);

    return 0;
}