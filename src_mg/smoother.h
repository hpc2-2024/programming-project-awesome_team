#ifndef SMOOTHER
#define SMOOTHER

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

/**
 * @brief Performs Jacobi smoothing on the input vector `u`.
 *
 * Applies Jacobi smoothing iterations to the input vector `u` based on the
 * right-hand side vector `f`. Supports both 1D and 2D based on `dim`.
 *
 * @param u Input vector to be smoothed.
 * @param f Right-hand side vector.
 * @param N Number of internal grid points.
 * @param v Number of smoothing iterations.
 * @param dim Dimension of the problem (1 or 2).
 *
 * @note Vectors `u` and `f` should be pre-allocated to the appropriate sizes.
 */
void smooth_jacobi(double u[], double f[], int N, int v, int dim) {
    int i, j, k;
    double h = 1.0/(N+1);
    int N_with_ghosts = N + 2;
    int vec_size = get_vec_size(N,dim,1);
    double *u_new = (double *)malloc(vec_size * sizeof(double));
    null_vec(u_new,vec_size);

    // Jacobi iteration
    for (k = 0; k < v; k++) {
        // Copy u to u_new (needed for the Jacobi update)
        memcpy(u_new, u, vec_size * sizeof(double));

        if (dim==2) {
            poisson_mat_vek(dim,N,u,u_new, 0);
            axpy(u_new,-1,u_new,f,vec_size);
            axpy(u_new,0.6/4*h*h,u_new,u,vec_size);

        } else if (dim==1) {
            poisson_mat_vek(dim,N,u,u_new, 0);
            axpy(u_new,-1,u_new,f,vec_size);
            axpy(u_new,0.6/2*h*h,u_new,u,vec_size);

        }

        // Swap u and u_new for the next iteration
        memcpy(u, u_new, vec_size * sizeof(double));
    }

    free(u_new);
}

/**
 * @brief Smooths the input vector `u` using Jacobi smoothing.
 *
 * Calls the Jacobi smoothing function with the given parameters.
 *
 * @param u Input vector to be smoothed.
 * @param f Right-hand side vector.
 * @param N Number of internal grid points.
 * @param v Number of smoothing iterations.
 * @param dim Dimension of the problem (1 or 2).
 *
 * @note Vectors `u` and `f` should be pre-allocated to the appropriate sizes.
 */
void smooth(double u[], double f[], int N, int v,int dim) {
    smooth_jacobi(u, f, N, v, dim);
}

#endif