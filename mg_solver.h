#ifndef MG_SOLVER
#define MG_SOLVER

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "src/poisson_mat_vek.h"
#include "src/smoother.h"
#include "src/restriction.h"
#include "src/prolongation.h"


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



void v_cycle(double** u, double **f, int N_start, int levels, int v, int dim, int debug){
    int vec_size;
    int N = N_start;

    for (int l = levels-1; l>=1; l--){
        vec_size = get_vec_size(N,dim,1);
        int N_coarser = dim_coarser(N); 
        int vec_size_coarser = get_vec_size(N_coarser,dim,1);

        double *r;
        r = malloc(vec_size * sizeof(*r));
        null_vec(r,vec_size);

        // Apply Smoothing
        smooth(u[l], f[l], N, v, dim);
        if (debug==1){
            printf("u_%d after smoothing:\n",l);
            vec_print(N,u[l],"u");
        }

        // Compute Residual
        poisson_mat_vek(dim,N,u[l],r);
        axpy(r, -1, r, f[l], vec_size);

        // Apply Restriction
        restriction_half(r, N, f[l-1], N_coarser,dim);
        if (debug==1){
            printf("f_%d after smoothing:\n",l-1);
            vec_print(N_coarser,f[l-1],"f");
        }

        // print_matrix(N_coarser + 2, f[l-1]);

        null_vec(u[l-1], vec_size_coarser);

        N = N_coarser;
        free(r);
    }

    exact_solve(u[0],f[0],N,dim);

    for (int l=1;l<levels;l++){
        vec_size = get_vec_size(N,dim,1);
        int N_finer = dim_finer(N);
        int vec_size_finer = get_vec_size(N_finer,dim,1);
        
        // init temporary vector  
        double *u_temp;
        u_temp = malloc(vec_size_finer * sizeof(*u_temp));
        null_vec(u_temp, vec_size_finer);

        // print_matrix(N+2, u[l-1]);

        // Apply Prolongation
        prolongation_simple(u[l-1], N, u_temp, N_finer,dim);

        // print_matrix(N_finer + 2, u_temp);

        // correction
        axpy(u[l], 1, u[l], u_temp, vec_size_finer);

        // Smoothing
        smooth(u[l], f[l], N_finer, v, dim);

        N = N_finer;
        free(u_temp);
    }
}

void mg_solve(double** u, double **f, int N, int levels,int v, int dim){
    int iter_max = 200;
    int iter = 0;
    double err;
    double epsilon=0.0001;

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
        
        vec_print(N,u[levels-1],"u_0");
        vec_print(N,f[levels-1],"f_0");
        vec_print(N,r,"r_0");
    }
    
    do {
        iter += 1;

        // Perform a V-cycle to update the solution
        v_cycle(u, f, N, levels, v, dim, debug);


        // clean up f
        int Nlevel = dim_coarser(N);
        for (int i=levels-2;i>=0;i--){
            int vec_size_level = get_vec_size(Nlevel,dim,1);
            null_vec(f[i], vec_size_level );
            Nlevel=dim_coarser(Nlevel);
        }

        // Update residual
        poisson_mat_vek(dim,N, u[levels - 1], r);                //r=Au
        axpy(r, -1, r, f[levels - 1], vec_size);    //r= f-r

        err = norm(r, vec_size);
        printf("Error: %f\n",err);
        if (iter>iter_max){
            printf("Multigrid solve stopped after %d iter without convergence.", iter);
            break;
        }
    } while (err > epsilon);
    
    free(r);
    // solution:
}

#endif