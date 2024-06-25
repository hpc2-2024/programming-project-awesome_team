#ifndef VCYCLE
#define VCYCLE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "poisson_mat_vek.h"
#include "smoother.h"
#include "restriction.h"
#include "prolongation.h"
#include "exact_solver.h"


/**
 * @brief Performs a V-cycle multigrid method for solving a linear system.
 *
 * Implements a V-cycle for solving a linear system using multigrid methods.
 * The method involves smoothing, restriction, exact solve on the coarsest grid, prolongation, and correction.
 *
 * @param u Array of solution vectors for each grid level.
 * @param f Array of right-hand side vectors for each grid level.
 * @param N_start Number of internal points in the finest grid.
 * @param levels Number of multigrid levels.
 * @param v Number of smoothing iterations at each level.
 * @param dim Dimension of the problem (1 or 2).
 * @param use_stencil9  Flag for enabling 9 point stencil in the 2d case (0 or 1).
 * @param debug Flag for enabling debug output (1 to enable, 0 to disable).
 *
 * @note Arrays `u` and `f` should be pre-allocated for each level.
 */

void go_up(){
    
}

void w_cycle(double** u, double** f, int N_start, int levels, int v, int dim, int use_stencil9, int debug) {
    int vec_size, vec_size_finer, vec_size_coarser, N_finer,N_coarser;
    int N = N_start;
    double *r, *u_temp;

    int count = 0; 
    for (int l = levels-1; l>=1; l--){
        printf("Level: %d \n", l);

        vec_size = get_vec_size(N,dim,1);
        N_coarser = dim_coarser(N); 
        vec_size_coarser = get_vec_size(N_coarser,dim,1);

        r = malloc(vec_size * sizeof(*r));
        null_vec(r,vec_size);

        // Apply Smoothing
        smooth(u[l], f[l], N, v, dim, use_stencil9);

        // Compute Residual
        poisson_mat_vek(dim,N,u[l],r, use_stencil9);
        axpy(r, -1, r, f[l], vec_size);

        // Apply Restriction
        restriction_half(r, N, f[l-1], N_coarser,dim);

        null_vec(u[l-1], vec_size_coarser);

        N = N_coarser;
        free(r);
    }

    exact_solve(u[0],f[0],N,dim, use_stencil9);
    count += 1; 

    // Peaks of w cycle in a for loop
    int peaks = 2*(levels-3)-1; // this is the result when breaking down the recursion into for loops

    for (int i = 0; i<peaks; i++){
        printf("Level: %d \n", i);
        vec_size = get_vec_size(N,dim,1);
        N_finer = dim_finer(N);
        vec_size_finer = get_vec_size(N_finer,dim,1);

        //Prolongation + Correction
        u_temp = malloc(vec_size_finer * sizeof(*u_temp));
        null_vec(u_temp, vec_size_finer);

        // Apply Prolongation
        prolongation_simple(u[0], N, u_temp, N_finer,dim);

        // correction
        axpy(u[1], 1, u[1], u_temp, vec_size_finer);

        N = N_finer;
        vec_size = vec_size_finer;

        // Smoothing
        smooth(u[1], f[1], N, v, dim, use_stencil9);

        N_coarser = dim_coarser(N); 
        vec_size_coarser = get_vec_size(N_coarser,dim,1);

        r = malloc(vec_size * sizeof(*r));
        null_vec(r,vec_size);

        // Compute Residual
        poisson_mat_vek(dim,N,u[1],r, use_stencil9);
        axpy(r, -1, r, f[1], vec_size);

        // Apply Restriction
        restriction_half(r, N, f[0], N_coarser,dim);

        null_vec(u[0], vec_size_coarser);
        
        N = N_coarser;
        vec_size = vec_size_coarser;

        exact_solve(u[0],f[0],N,dim, use_stencil9);
        count += 1; 

        free(u_temp);
        free(r);

        int peak_height = (levels-2) - abs(levels-2 -2 +i); // this could be wrong

        for (int l=1;l<peak_height+1;l++){
            printf("Level: %d \n", l);

            vec_size = get_vec_size(N,dim,1);
            N_finer = dim_finer(N);
            vec_size_finer = get_vec_size(N_finer,dim,1);
            
            // init temporary vector  
            u_temp = malloc(vec_size_finer * sizeof(*u_temp));
            null_vec(u_temp, vec_size_finer);

            // Apply Prolongation
            prolongation_simple(u[l-1], N, u_temp, N_finer,dim);

            // correction
            axpy(u[l], 1, u[l], u_temp, vec_size_finer);

            // Smoothing
            smooth(u[l], f[l], N_finer, v, dim, use_stencil9);

            N = N_finer;
            free(u_temp);
        }

        for (int l = peak_height; l>=1; l--){
            printf("Level: %d \n", l);

            vec_size = get_vec_size(N,dim,1);
            N_coarser = dim_coarser(N); 
            vec_size_coarser = get_vec_size(N_coarser,dim,1);

            r = malloc(vec_size * sizeof(*r));
            null_vec(r,vec_size);

            // Apply Smoothing (but not at the peak level, since we already did smooth there)
            if (l<peak_height) {
                smooth(u[l], f[l], N, v, dim, use_stencil9);
            }

            // Compute Residual
            poisson_mat_vek(dim,N,u[l],r, use_stencil9);
            axpy(r, -1, r, f[l], vec_size);

            // Apply Restriction
            restriction_half(r, N, f[l-1], N_coarser,dim);

            null_vec(u[l-1], vec_size_coarser);

            N = N_coarser;
            free(r);
        }

        exact_solve(u[0],f[0],N,dim, use_stencil9);
        printf("Level: 0 \n");

        count += 1; 

    }

    vec_size = get_vec_size(N,dim,1);
    N_finer = dim_finer(N);
    vec_size_finer = get_vec_size(N_finer,dim,1);

    //Prolongation + Correction
    u_temp = malloc(vec_size_finer * sizeof(*u_temp));
    null_vec(u_temp, vec_size_finer);

    // Apply Prolongation
    prolongation_simple(u[0], N, u_temp, N_finer,dim);

    // correction
    axpy(u[1], 1, u[1], u_temp, vec_size_finer);

    free(u_temp);
    N = N_finer;
    vec_size = vec_size_finer;

    // Smoothing
    smooth(u[1], f[1], N, v, dim, use_stencil9);
    printf("Level: 1 \n");


    r = malloc(vec_size * sizeof(*r));
    null_vec(r,vec_size);

    N_coarser = dim_coarser(N); 
    vec_size_coarser = get_vec_size(N_coarser,dim,1);

    // Compute Residual
    poisson_mat_vek(dim,N,u[1],r, use_stencil9);
    axpy(r, -1, r, f[1], vec_size);

    // Apply Restriction
    restriction_half(r, N, f[0], N_coarser,dim);

    null_vec(u[0], vec_size_coarser);

    free(r);
    N = N_coarser;

    exact_solve(u[0],f[0],N,dim, use_stencil9);
    printf("Level: 0 \n");

    count += 1; 

    for (int l=1;l<levels;l++){
        printf("Level: %d \n", l);

        vec_size = get_vec_size(N,dim,1);
        N_finer = dim_finer(N);
        vec_size_finer = get_vec_size(N_finer,dim,1);
        
        // init temporary vector  
        u_temp = malloc(vec_size_finer * sizeof(*u_temp));
        null_vec(u_temp, vec_size_finer);

        // Apply Prolongation
        prolongation_simple(u[l-1], N, u_temp, N_finer,dim);

        // correction
        axpy(u[l], 1, u[l], u_temp, vec_size_finer);

        // Smoothing
        smooth(u[l], f[l], N_finer, v, dim, use_stencil9);

        N = N_finer;
        free(u_temp);
    }
    int i = 0;

}

#endif