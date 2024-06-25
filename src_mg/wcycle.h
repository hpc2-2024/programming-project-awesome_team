#ifndef WCYCLE
#define WCYCLE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "poisson_mat_vek.h"
#include "smoother.h"
#include "restriction.h"
#include "prolongation.h"
#include "exact_solver.h"


void w_cycle_rec(double **u, double **f, int N_start, int levels, int v, int dim, int use_stencil9, int debug) {
    int N = N_start;
    if (levels==1){
        exact_solve(u[0],f[0],N,dim,use_stencil9);
    }
    else {
        int vec_size = get_vec_size(N,dim,1);
        int N_coarser = dim_coarser(N); 
        int vec_size_coarser = get_vec_size(N_coarser,dim,1);

        double *r;
        r = malloc(vec_size * sizeof(*r));
        null_vec(r,vec_size);

        // Apply Smoothing
        smooth(u[levels-1], f[levels-1], N, v, dim, use_stencil9);

        // Compute Residual
        poisson_mat_vek(dim,N,u[levels-1],r, use_stencil9);
        axpy(r, -1, r, f[levels-1], vec_size);

        // Apply Restriction
        restriction_half(r, N, f[levels-2], N_coarser,dim);

        free(r);
        null_vec(u[levels-2], vec_size_coarser);

        if (levels>2) {
            w_cycle_rec(u,f, N_coarser, levels-1, v, dim, use_stencil9, debug);
            w_cycle_rec(u,f, N_coarser, levels-1, v, dim, use_stencil9, debug);
        }
        else {
            w_cycle_rec(u,f, N_coarser, levels-1, v, dim, use_stencil9, debug);
        }

        
        // init temporary vector  
        double* u_temp = malloc(vec_size * sizeof(*u_temp));
        null_vec(u_temp, vec_size);

        // Apply Prolongation
        prolongation_simple(u[levels-2], N_coarser, u_temp, N, dim);

        // correction
        axpy(u[levels-1], 1, u[levels-1], u_temp, vec_size);

        // Smoothing
        smooth(u[levels-1], f[levels-1], N, v, dim, use_stencil9);

        free(u_temp);

    }

}


#endif