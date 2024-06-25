#ifndef MG_SOLVER
#define MG_SOLVER

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vcycle.h"
#include "poisson_mat_vek.h"


/**
 * @brief Solves a linear system using the multigrid method.
 *
 * Solves the linear system by performing multiple V-cycles until convergence or maximum iterations are reached.
 *
 * @param u The solution vectors for each level of the grid hierarchy.
 * @param f The right-hand side vectors for each level of the grid hierarchy.
 * @param N The number of internal points in the finest grid.
 * @param levels The number of levels in the multigrid hierarchy.
 * @param v The number of pre- and post-smoothing steps.
 * @param dim The dimension of the problem (1 or 2).
 * @param use_stencil9  Flag for enabling 9 point stencil in the 2d case (0 or 1).
 * @param print_error Option to print the error after each cycle (0 or 1).
 */
void mg_solve(double** u, double **f, int N, int levels,int v, int dim, int fcycle,
                 int use_stencil9, int print_error, int smoother, int *num_iterations, double *final_error, int *converged){
    int iter_max = 200;
    int iter = 0;
    double err;
    double epsilon=0.0001;

    *converged = 1;

    int vec_size = get_vec_size(N,dim,1);

    double *r;
    r = malloc(vec_size * sizeof(*r));
    null_vec(r, vec_size);

    int debug = 0;

    if (debug==1){
        printf("\n Debug Mg method \n");
        printf("Grid size N: %d\n",N);
        iter_max = 1;
        printf("Max v cycle iterations: %d\n",iter_max);
        
        print_vector(u[levels-1],N,dim,1, "u after 1 iteration (finest level)");
        print_vector(f[levels-1],N,dim,1, "f after 1 iteration (finest level)");
        print_vector(r,N,dim,1, "residuum after 1 iteration: ");
    }
    
    do {
        iter += 1;

        // Perform a V-cycle to update the solution
        if (fcycle == 1) {
            f_cycle(u, f, N, levels, v, dim, use_stencil9, debug, smoother);
        }
        else {
            v_cycle(u, f, N, levels, v, dim, use_stencil9,debug, smoother);
        }


        // clean up f
        int Nlevel = dim_coarser(N);
        for (int i=levels-2;i>=0;i--){
            int vec_size_level = get_vec_size(Nlevel,dim,1);
            null_vec(f[i], vec_size_level );
            Nlevel=dim_coarser(Nlevel);
        }

        // Update residual
        poisson_mat_vek(dim,N, u[levels - 1], r, use_stencil9);                //r=Au
        axpy(r, -1, r, f[levels - 1], vec_size);    //r= f-r

        err = norm(r, vec_size);
        if (print_error==1) {
            printf("Error: %f\n",err);
        }
        if (iter>iter_max){
            *converged = 0;
            // printf("Multigrid solve stopped after %d iter without convergence.", iter);
            break;
        }
    } while (err > epsilon);
    
    *num_iterations = iter;
    *final_error = err;
    free(r);
    // solution:
}

#endif