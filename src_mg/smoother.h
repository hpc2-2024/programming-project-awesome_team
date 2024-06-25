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
                poisson_mat_vek(dim,N,u,u_new, 0);
                axpy(u_new,-1,u_new,f,vec_size);
                axpy(u_new,0.6/4*h*h,u_new,u,vec_size);
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

/*!
 * @brief Performs Gauss-Seidel smoothing on the input vector `u`.
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
void smooth_gauss_seidel(double *X, double *B, int N, int v, int dim, int use_stencil9){
    int N_pad = N + 2; 
    int X_size = get_vec_size(N, dim, 1); // (N+2) * (N+2)
    
    double h = 1.0/(N+1);
    double h2 = h*h;

    // double *X_new;
    // X_new = malloc(X_size * sizeof(*X_new));
    int iter;
    if(use_stencil9 == 0){
        for (iter = 0; iter < v; iter++) {

            for (int i = 1; i < N_pad-1; i++){
                for(int j = 1; j < N_pad; j++){
                    double sum = 0.0;
                    
                    sum = - X[i * N_pad + (j+1)] 
                        - X[i * N_pad + (j-1)] 
                        - X[(i+1) * N_pad + j] 
                        - X[(i-1) * N_pad + j];
                    
                    sum = sum / h2;
                    X[i * N_pad + j] = (B[i * N_pad + j] - sum) / (4 / h2);
                }
            }
        }
    }
    else{
        for (iter = 0; iter < v; iter++) {

            for (int i = 1; i < N_pad-1; i++){
                for(int j = 1; j < N_pad; j++){
                    double sum = 0.0;
                    
                    sum =   -0.5*X[N_pad*i+j-1] 
                            -0.5*X[N_pad*i+j+1] 
                            -0.5*X[N_pad*(i+1)+j]
                            -0.5*X[(N+2)*(i-1)+j]
                            -0.25*X[N_pad*(i+1)+j+1] 
                            -0.25*X[N_pad*(i+1)+j-1]  
                            -0.25*X[N_pad*(i-1)+j+1] 
                            -0.25*X[N_pad*(i-1)+j-1];
                    
                    sum = sum / (h2 * 6.0);

                    X[i * N_pad + j] = (B[i * N_pad + j] - sum) / (3 / (h2 * 6.0));
                }
            }
        }
    }
}

/*!
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
void smooth(double u[], double f[], int N, int v,int dim, int use_stencil9, int smoother) {
    
    if(smoother == 0){
        smooth_jacobi(u, f, N, v, dim, use_stencil9);
    }
    else if(smoother == 1){
        smooth_gauss_seidel(u, f, N, v, dim, use_stencil9);
    }
}

#endif