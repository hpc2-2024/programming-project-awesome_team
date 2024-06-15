#ifndef utils
#define utils
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int dim_finer(int N){
    return N*2 + 1;
}

int dim_coarser(int N){
    return (N-1)/2;
}

// Function to perform restriction (fine grid to coarse grid)
void restriction(double* fine, int N_fine, double* coarse, int N_coarse) {
    // N_fine = N, N_coarse = (N-1)/2 (because we assume N_fine = 2 * N_coarse + 1)
    for (int i = 1; i <= N_coarse; i++) {
        for (int j = 1; j <= N_coarse; j++) {
            coarse[i * (N_coarse + 2) + j] = 0.25 * (
                fine[2 * i * (N_fine + 2) + 2 * j] +
                0.5 * (fine[2 * i * (N_fine + 2) + 2 * j - 1] + fine[2 * i * (N_fine + 2) + 2 * j + 1]) +
                0.5 * (fine[(2 * i - 1) * (N_fine + 2) + 2 * j] + fine[(2 * i + 1) * (N_fine + 2) + 2 * j]) +
                0.25 * (fine[(2 * i - 1) * (N_fine + 2) + 2 * j - 1] + fine[(2 * i - 1) * (N_fine + 2) + 2 * j + 1] +
                        fine[(2 * i + 1) * (N_fine + 2) + 2 * j - 1] + fine[(2 * i + 1) * (N_fine + 2) + 2 * j + 1])
            );
        }
    }
}

// Function to perform prolongation (coarse grid to fine grid)
void prolongation(double* coarse, int N_coarse, double* fine, int N_fine) {
    // N_fine = N, N_coarse = (N-1)/2 (because we assume N_fine = 2 * N_coarse + 1)
    for (int i = 1; i <= N_coarse; i++) {
        for (int j = 1; j <= N_coarse; j++) {
            fine[2 * i * (N_fine + 2) + 2 * j] = coarse[i * (N_coarse + 2) + j];
            fine[2 * i * (N_fine + 2) + 2 * j - 1] = 0.5 * (coarse[i * (N_coarse + 2) + j] + coarse[i * (N_coarse + 2) + j - 1]);
            fine[2 * i * (N_fine + 2) + 2 * j + 1] = 0.5 * (coarse[i * (N_coarse + 2) + j] + coarse[i * (N_coarse + 2) + j + 1]);
            fine[(2 * i - 1) * (N_fine + 2) + 2 * j] = 0.5 * (coarse[i * (N_coarse + 2) + j] + coarse[(i - 1) * (N_coarse + 2) + j]);
            fine[(2 * i + 1) * (N_fine + 2) + 2 * j] = 0.5 * (coarse[i * (N_coarse + 2) + j] + coarse[(i + 1) * (N_coarse + 2) + j]);
            fine[(2 * i - 1) * (N_fine + 2) + 2 * j - 1] = 0.25 * (
                coarse[i * (N_coarse + 2) + j] + coarse[i * (N_coarse + 2) + j - 1] +
                coarse[(i - 1) * (N_coarse + 2) + j] + coarse[(i - 1) * (N_coarse + 2) + j - 1]
            );
            fine[(2 * i - 1) * (N_fine + 2) + 2 * j + 1] = 0.25 * (
                coarse[i * (N_coarse + 2) + j] + coarse[i * (N_coarse + 2) + j + 1] +
                coarse[(i - 1) * (N_coarse + 2) + j] + coarse[(i - 1) * (N_coarse + 2) + j + 1]
            );
            fine[(2 * i + 1) * (N_fine + 2) + 2 * j - 1] = 0.25 * (
                coarse[i * (N_coarse + 2) + j] + coarse[i * (N_coarse + 2) + j - 1] +
                coarse[(i + 1) * (N_coarse + 2) + j] + coarse[(i + 1) * (N_coarse + 2) + j - 1]
            );
            fine[(2 * i + 1) * (N_fine + 2) + 2 * j + 1] = 0.25 * (
                coarse[i * (N_coarse + 2) + j] + coarse[i * (N_coarse + 2) + j + 1] +
                coarse[(i + 1) * (N_coarse + 2) + j] + coarse[(i + 1) * (N_coarse + 2) + j + 1]
            );
        }
    }
}

// Function to perform Jacobi smoothing
void smooth_jacobi(double u[], double f[], int N, int v) {
    int i, j, k;
    int N_with_ghosts = N + 2;
    double *u_new = (double *)malloc(N_with_ghosts * N_with_ghosts * sizeof(double));

    // Jacobi iteration
    for (k = 0; k < v; k++) {
        // Copy u to u_new (needed for the Jacobi update)
        memcpy(u_new, u, N_with_ghosts * N_with_ghosts * sizeof(double));

        for (i = 1; i <= N; i++) {
            for (j = 1; j <= N; j++) {
                u_new[i * N_with_ghosts + j] = 0.25 * (
                    u[(i-1) * N_with_ghosts + j] + // up
                    u[(i+1) * N_with_ghosts + j] + // down
                    u[i * N_with_ghosts + (j-1)] + // left
                    u[i * N_with_ghosts + (j+1)]   // right
                    + f[i * N_with_ghosts + j]     // source term
                );
            }
        }

        // Swap u and u_new for the next iteration
        memcpy(u, u_new, N_with_ghosts * N_with_ghosts * sizeof(double));
    }

    free(u_new);
}

void smooth(double u[], double f[], int N, int v) {
    smooth_jacobi(u, f, N, v);
}

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
void exact_solve(double u[], double f[], int N) {
    int i, j, k;
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
            A[idx * NN + idx] = 4.0;
            if (i > 1) A[idx * NN + (idx - N)] = -1.0; // up
            if (i < N) A[idx * NN + (idx + N)] = -1.0; // down
            if (j > 1) A[idx * NN + (idx - 1)] = -1.0; // left
            if (j < N) A[idx * NN + (idx + 1)] = -1.0; // right
        }
    }

    // Use a simple Gaussian elimination method (or call an existing library function)
    // to solve A * temp_u = b

    // This is a placeholder for a solver function
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

void v_cycle(double** u, double **f, int N, int levels,int v,int debug){
    int vec_size_start = (N+2)*(N+2);
    int vec_size;
    
    int Nlevel= N;
    for (int l=levels-1;l>=1;l--){
        // number of points in the grid (without border)
        vec_size = (2+Nlevel)*(Nlevel+2);

        double *r;
        r = (double*)malloc((vec_size)*sizeof(double));
        null_vec(r,vec_size);

        // Smoothing
        smooth(u[l],f[l],Nlevel,v);
        if (debug==1){
            printf("u_%d after smoothing:\n",l);
            vec_print(Nlevel,u[l],"u");
        }

        // residual
        mfMult(Nlevel,u[l],r);              //r=Au_l
        axpy(r, -1, r, f[l],vec_size);      //r=f_l-Au_l

        // restriction
        int N_f = dim_coarser(Nlevel); 
        restriction(r,Nlevel, f[l-1],N_f);
        Nlevel=N_f;
        if (debug==1){
            printf("f_%d after smoothing:\n",l-1);
            vec_print(Nlevel,f[l-1],"f");
        }

        null_vec(u[l-1],(Nlevel + 2) * (Nlevel + 2));

        free(r);
    }

    printf("Nlevel of coarsest grid: %d \n",Nlevel);
    // exact solve
    exact_solve(u[0],f[0],Nlevel);

    for (int l=1;l<levels;l++){
        int Nlevel_smaller = Nlevel;
        Nlevel = dim_finer(Nlevel);
        double *temp_u = (double *)malloc(pow(Nlevel+2,2) * sizeof(double));
        null_vec(temp_u,pow(Nlevel+2,2));

        // prolongate
        prolongation(u[l-1],Nlevel_smaller,temp_u,Nlevel);
        // correction
        axpy(u[l],1,u[l],temp_u,pow(Nlevel+2,2));
        // Smoothing
        smooth(u[l],f[l],Nlevel,v);

        free(temp_u);
    }
}

int get_vec_size(int N, int dim, int ghostlayer){
    if (ghostlayer==1){
        N = N+2;
    }
    int vec_size = pow(N,dim);
    return vec_size;
}

void mg_solve(double** u, double **f, int N, int levels,int v, int dim){
    
    // setup
    int it_max = 99;
    int iterations = 0;

    double err;
    double epsilon=0.0001;

    int vec_size = get_vec_size(N,dim,1);

    double* r =(double*)malloc(vec_size*sizeof(double));
    null_vec(r,vec_size);

    int debug = 1;

    if (debug==1){
        printf("\n Debug Mg method \n");
        printf("Grid size N: %d\n",N);
        it_max = 1;
        printf("Max v cycle iterations: %d\n",it_max);
        
        vec_print(N,u[levels-1],"u_0");
        vec_print(N,f[levels-1],"f_0");
        vec_print(N,r,"r_0");
    }
    
    do {
        iterations += 1;

        // Perform a V-cycle to update the solution
        v_cycle(u, f, N, levels, v,debug);


        // clean up f
        int Nlevel = dim_coarser(N);
        for (int i=levels-2;i>=0;i--){
            int vec_size_level = get_vec_size(N,dim,1);
            null_vec(f[i], vec_size_level );
            Nlevel=dim_coarser(Nlevel);
        }

        // calculate new residual
        mfMult(N, u[levels - 1], r);                //Au
        axpy(r, -1, r, f[levels - 1], vec_size);    //r= f-Au

        err = norm(r,vec_size);
        printf("Error: %f\n",err);
        if (iterations>it_max){
            printf("Multigrid solve stopped after %d iterations without convergence.", iterations);
            break;
        }
    } while ( err >epsilon);
    
    free(r);
    // solution:
}

