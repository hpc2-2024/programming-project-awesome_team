#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "math.h"
#include "utils.h"

// Simple Gaussian elimination solver for dense systems
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

// Function to solve the 2D Poisson problem with ghost layers
void exact_solve_poisson_2D(double u[], double f[], int N) {
    int i, j, k;
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

void exact_solve(double u[], double f[], int N, int dim){
    if (dim==2){
        exact_solve_poisson_2D(u,f,N);
    }
    else if (dim==1){
        exact_solve_poisson_1D(u,f,N);
    }
}