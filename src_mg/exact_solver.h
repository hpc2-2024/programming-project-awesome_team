#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "math.h"
#include "utils.h"

/**
 * @brief Performs Gaussian elimination to solve a dense linear system.
 *
 * Solves the system of linear equations Ax = b using Gaussian elimination with partial pivoting.
 *
 * @param A The coefficient matrix of size n x n (no ghostlayer).
 * @param b The right-hand side vector of size n (no ghostlayer).
 * @param x The solution vector of size n.
 * @param n The number of equations (size of the matrix and vectors).
 */
void gaussian_elimination(double A[], double b[], double x[], int n) {
    int i, j, k;
    for (i = 0; i < n; i++) {
        // Pivoting
        for (k = i + 1; k < n; k++) {
            if (fabs(A[k * n + i]) > fabs(A[i * n + i])) {
                for (j = 0; j < n; j++) {
                    double temp = A[i * n + j];
                    A[i * n + j] = A[k * n + j];
                    A[k * n + j] = temp;
                }
                double temp = b[i];
                b[i] = b[k];
                b[k] = temp;
            }
        }

        // Elimination
        for (k = i + 1; k < n; k++) {
            double t = A[k * n + i] / A[i * n + i];
            for (j = 0; j < n; j++) {
                A[k * n + j] -= t * A[i * n + j];
            }
            b[k] -= t * b[i];
        }
    }

    // Back substitution
    for (i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (j = i + 1; j < n; j++) {
            x[i] -= A[i * n + j] * x[j];
        }
        x[i] /= A[i * n + i];
    }
}

/**
 * @brief Solves the 2D Poisson equation using Gaussian elimination.
 *
 * Solves the 2D Poisson equation on a grid with ghost layers by constructing the linear system
 * and solving it using Gaussian elimination.
 *
 * @param u The solution vector with ghost layers.
 * @param f The right-hand side vector with ghost layers.
 * @param N The number of internal points in each dimension.
 * @param use_stencil9  Flag for enabling 9 point stencil in the 2d case (0 or 1).
 */
void exact_solve_poisson_2D(double u[], double f[], int N, int use_stencil9) {
    int i, j, k, idx;
    double h = 1.0/(N+1);
    double h2 = h*h;
    int NN = N * N;

    // Create matrix A (NN x NN) and vector b (NN)
    double *A = (double *)malloc(NN * NN * sizeof(double));
    double *b = (double *)malloc(NN * sizeof(double));
    double *temp_u = (double *)malloc(NN * sizeof(double)); // Temporary array for solution
    
    // Initialize A to be the Laplacian matrix and b to be the right-hand side vector
    for (i = 0; i < NN; i++) {
        b[i] = 0.0;
        for (j = 0; j < NN; j++) {
            A[i * NN + j] = 0.0;
        }
    }

    if (use_stencil9==1) {

        double h_sq = (1.0/6.0) * pow(1.0/h, 2);
        // Fill matrix A with the 9-point stencil
        for (i = 1; i <= N; i++) {
            for (j = 1; j <= N; j++) {
                idx = (i - 1) * N + (j - 1);
                b[idx] = f[i * (N+2) + j];

                // Diagonal element
                A[idx * NN + idx] = 3 * h_sq;

                // 9-point stencil elements
                if (j > 1) {
                    A[idx * NN + (idx - 1)] = -0.5 * h_sq; // left
                }
                if (j < N) {
                    A[idx * NN + (idx + 1)] = -0.5 * h_sq; // right
                }
                if (i > 1) {
                    A[idx * NN + (idx - N)] = -0.5 * h_sq; // up
                }
                if (i < N) {
                    A[idx * NN + (idx + N)] = -0.5 * h_sq; // down
                }
                if (i > 1 && j > 1) {
                    A[idx * NN + (idx - N - 1)] = - 0.25 * h_sq; // up-left
                }
                if (i > 1 && j < N) {
                    A[idx * NN + (idx - N + 1)] = -0.25 * h_sq; // up-right
                }
                if (i < N && j > 1) {
                    A[idx * NN + (idx + N - 1)] = -0.25 * h_sq; // down-left
                }
                if (i < N && j < N) {
                    A[idx * NN + (idx + N + 1)] = -0.25 * h_sq; // down-right
                }
            }
        }

    }
    else {

        // Fill matrix A with the discretized Laplacian operator
        for (i = 1; i <= N; i++) {
            for (j = 1; j <= N; j++) {
                int idx = (i-1) * N + (j-1);
                b[idx] = f[i * (N + 2) + j];
                A[idx * NN + idx] =  4.0/h2;
                if (i > 1) A[idx * NN + (idx - N)] = -1.0/h2; // up
                if (i < N) A[idx * NN + (idx + N)] = -1.0/h2; // down
                if (j > 1) A[idx * NN + (idx - 1)] = -1.0/h2; // left
                if (j < N) A[idx * NN + (idx + 1)] = -1.0/h2; // right
            }
        }

    }

    gaussian_elimination(A, b, temp_u, NN);

    // Map the solution back to the u array with ghost layers
    for (i = 1; i <= N; i++) {
        for (j = 1; j <= N; j++) {
            int idx = (i-1) * N + (j-1);
            u[i * (N + 2) + j] = temp_u[idx];
        }
    }

    // Free allocated memory
    free(A);
    free(b);
    free(temp_u);
}

/**
 * @brief Solves the 1D Poisson equation using Gaussian elimination.
 *
 * Solves the 1D Poisson equation on a grid with ghost layers by constructing the linear system
 * and solving it using Gaussian elimination.
 *
 * @param u The solution vector with ghost layers.
 * @param f The right-hand side vector with ghost layers.
 * @param N The number of internal points.
 */
void exact_solve_poisson_1D(double u[], double f[], int N){
    int NN = N * N;
    double h = 1.0/(N+1);
    double h2 = h*h;
    double *A = (double *)malloc(NN * sizeof(double));
    double *b = (double *)malloc(N * sizeof(double));
    double *temp_u = (double *)malloc(N * sizeof(double)); // Temporary array for solution
    null_vec(A,NN);

    //Fill matrix A with the discretized Laplacian operator
    for (int i=0;i<N;i++){
        A[i*N + i]=2/h2;
        b[i]=f[i+1];
    }
    for (int i=0;i<N-1;i++){
        A[i*N +i+1]=-1/h2;
        A[(i+1)*N + i]=-1/h2;
    }

    gaussian_elimination(A, b, temp_u, N);

    // Map the solution back to the array u with ghost layer
    for (int i=1;i<=N;i++) {
        u[i]=temp_u[i-1];
    }

    free(A);
    free(b);
    free(temp_u);

}

/**
 * @brief Solves the Poisson equation (1D or 2D) using Gaussian elimination.
 *
 * Determines the dimension of the problem and calls the appropriate function to solve
 * the Poisson equation using Gaussian elimination.
 *
 * @param u The solution vector with ghost layers.
 * @param f The right-hand side vector with ghost layers.
 * @param N The number of internal points.
 * @param dim Dimension of the problem (1 or 2).
 * @param use_stencil9  Flag for enabling 9 point stencil in the 2d case (0 or 1).
 */
void exact_solve(double u[], double f[], int N, int dim, int use_stencil9){
    if (dim==2){
        exact_solve_poisson_2D(u,f,N, use_stencil9);
    }
    else if (dim==1){
        exact_solve_poisson_1D(u,f,N);
    }
}