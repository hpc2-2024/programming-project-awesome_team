#ifndef MG_SOLVER
#define MG_SOLVER

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int dim_finer(int N){
    return N*2 - 1;
}

int dim_coarser(int M){
    return (M+1)/2;
}

void restriction_simple(double *fine_grid, int M, double *coarse_grid, int N){
    int M_pad = M+2;
    int N_pad = N+2;

    for(int i = 0; i<N_pad; i++){
        for(int j = 0; j<N_pad; j++){
            int k = 2*i-1;
            int l = 2*j-1;

            coarse_grid[i * N_pad + j] = fine_grid[k * M_pad + l];

        }
    }
}

void restriction_half(double *fine_grid, int M, double *coarse_grid, int N){
    int M_pad = M+2;
    int N_pad = N+2;

    for(int i = 1; i<N_pad-1; i++){
        for(int j = 1; j<N_pad-1; j++){
            int k = 2*i-1;
            int l = 2*j-1;

            coarse_grid[i * N_pad + j] = 0.125  * (
                fine_grid[(k+1) * M_pad + l]
                + fine_grid[(k-1) * M_pad + l]
                + fine_grid[k * M_pad + (l-1)] 
                + fine_grid[k * M_pad + (l+1)]) 
                + 0.5 * fine_grid[k * M_pad + l];
        }
    }
}

// expects a padded coarse grid of size (N+2) x (N+2) and an empty, padded fine grid of size (M+2) x (M+2)
void prolongation_simple(double *coarse_grid, int N, double* fine_grid, int M){
    int N_pad = N + 2;
    int M_pad = M + 2;

    // only iterate over the inner points and interpolate them from the coarse grid
    for(int i = 1; i<M_pad-1; i++){
        for(int j = 1; j<M_pad-1; j++){
            int k = (int) (i+1)/2;
            int l = (int) (j+1)/2;

            if(i%2 == 1 && j%2==1){
                fine_grid[i * M_pad + j] = coarse_grid[k * N_pad + l];
            }
            else if(i%2 == 0 && j%2==1){
                fine_grid[i * M_pad + j] = 0.5 * coarse_grid[k * N_pad + l] + 0.5 * coarse_grid[(k+1) * N_pad + l];

            }
            else if(i%2 == 1 && j%2==0){
                fine_grid[i * M_pad + j] = 0.5 * coarse_grid[k * N_pad + l] + 0.5 * coarse_grid[k * N_pad + (l+1)];

            }
            else if(i%2 == 0 && j%2==0){
                fine_grid[i * M_pad + j]=  0.25 * coarse_grid[k * N_pad + l] 
                                        + 0.25 * coarse_grid[(k+1) * N_pad + l]
                                        + 0.25 * coarse_grid[k * N_pad + (l+1)]
                                        + 0.25 * coarse_grid[(k+1) * N_pad + (l+1)];
            }
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

void v_cycle(double** u, double **f, int N_start, int levels, int v, int debug){
    int vec_size;
    int N = N_start;

    for (int l = levels-1; l>=1; l--){
        vec_size = (N + 2) * (N + 2);
        int N_coarser = dim_coarser(N); 
        int vec_size_coarser = (N_coarser + 2) * (N_coarser + 2);

        double *r;
        r = malloc(vec_size * sizeof(*r));
        null_vec(r,vec_size);

        // Apply Smoothing
        smooth(u[l], f[l], N, v);
        if (debug==1){
            printf("u_%d after smoothing:\n",l);
            vec_print(N,u[l],"u");
        }

        // Compute Residual
        mfMult(N,u[l],r);
        axpy(r, -1, r, f[l], vec_size);

        // Apply Restriction
        restriction_half(r, N, f[l-1], N_coarser);
        if (debug==1){
            printf("f_%d after smoothing:\n",l-1);
            vec_print(N_coarser,f[l-1],"f");
        }

        // print_matrix(N_coarser + 2, f[l-1]);

        null_vec(u[l-1], vec_size_coarser);

        N = N_coarser;
        free(r);
    }

    exact_solve(u[0],f[0],N);

    for (int l=1;l<levels;l++){
        vec_size = (N + 2) * (N + 2);
        int N_finer = dim_finer(N);
        int vec_size_finer = (N_finer + 2) * (N_finer + 2);
        
        // init temporary vector  
        double *u_temp;
        u_temp = malloc(vec_size_finer * sizeof(*u_temp));
        null_vec(u_temp, vec_size_finer);

        // print_matrix(N+2, u[l-1]);

        // Apply Prolongation
        prolongation_simple(u[l-1], N, u_temp, N_finer);

        // print_matrix(N_finer + 2, u_temp);

        // correction
        axpy(u[l], 1, u[l], u_temp, vec_size_finer);

        // Smoothing
        smooth(u[l], f[l], N_finer, v);

        N = N_finer;
        free(u_temp);
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
        v_cycle(u, f, N, levels, v,debug);


        // clean up f
        int Nlevel = dim_coarser(N);
        for (int i=levels-2;i>=0;i--){
            int vec_size_level = get_vec_size(Nlevel,dim,1);
            null_vec(f[i], vec_size_level );
            Nlevel=dim_coarser(Nlevel);
        }

        // Update residual
        poisson_mat_vek(dim,N, u[levels - 1], r);                //Au
        axpy(r, -1, r, f[levels - 1], vec_size);    //r= f-Au

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
