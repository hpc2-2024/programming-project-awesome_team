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
 * @param use_stencil9  Flag for enabling 9 point stencil in the 2d case (0 or 1).
 *
 * @note Vectors `u` and `f` should be pre-allocated to the appropriate sizes.
 */
void smooth_jacobi(double u[], double f[], int N, int v, int dim, int use_stencil9) {
    int i, j, k;
    double h = 1.0/(N+1);
    int N_with_ghosts = N + 2;
    int vec_size = get_vec_size(N,dim,1);
    double *u_new = (double *)malloc(vec_size * sizeof(double));
    null_vec(u_new,vec_size);
    double omega = 0.6/4*h*h;

    // Jacobi iteration
    for (k = 0; k < v; k++) {
        // Copy u to u_new (needed for the Jacobi update)
        memcpy(u_new, u, vec_size * sizeof(double));

        if (dim==2) {
            if (use_stencil9 == 1) {
                poisson_mat_vek(dim,N,u,u_new, use_stencil9);
                axpy(u_new,-1,u_new,f,vec_size);
                axpy(u_new,0.6/3 * 6 *h*h,u_new,u,vec_size);

            }
            else {
                poisson_mat_vek(dim, N, u, u_new, 0); // u_new = A u
                axpy(u_new, -1, u_new, f, vec_size); // u_new = f - A u
                axpy(u_new, omega, u_new, u, vec_size); // u_new = u + omega * (f - A u)
            }

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
 * @brief Performs Richardson smoothing on the input vector `u`.
 *
 * Quite similar to the Jacobi method, except that we change parameter by which we weigh the update in each iteration.
 *
 * @param u Input vector to be smoothed.
 * @param f Right-hand side vector.
 * @param N Number of internal grid points.
 * @param v Number of smoothing iterations.
 * @param dim Dimension of the problem (1 or 2).
 * @param use_stencil9  Flag for enabling 9 point stencil in the 2d case (0 or 1).
 *
 * @note Vectors `u` and `f` should be pre-allocated to the appropriate sizes.
 */
void smooth_richardson(double u[], double f[], int N, int v, int dim, int use_stencil9) {
    int i, j, k;
    double h = 1.0/(N+1);
    double h2 = h * h;
    int N_with_ghosts = N + 2;
    int vec_size = get_vec_size(N,dim,1);
    double *u_new = (double *)malloc(vec_size * sizeof(double));
    null_vec(u_new,vec_size);

    double omega = 0.3 / 4 * h * h; // algorithm is pretty sensitive to the weight
    
    // Modified Richardson Iteration (very similar to Jacobi)
    // x_{k+1} = x_k + omega * (f - A x_k)
    for (k = 0; k < v; k++) {
        // Copy u to u_new (needed for the Jacobi update)
        memcpy(u_new, u, vec_size * sizeof(double));

        poisson_mat_vek(dim, N, u, u_new, 0); // u_new = A u 
        axpy(u_new, -1.0, u_new, f, vec_size); // u_new = f - A u  
        axpy(u_new, omega, u_new, u, vec_size); // u_new = u - omega * (A u + f)

        // Swap u and u_new for the next iteration
        memcpy(u, u_new, vec_size * sizeof(*u));
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
void smooth(double u[], double f[], int N, int v,int dim, int use_stencil9) {
    smooth_jacobi(u, f, N, v, dim, use_stencil9);
    // smooth_richardson(u, f, N, v, dim, use_stencil9);
}

#endif