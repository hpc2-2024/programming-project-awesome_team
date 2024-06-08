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

/*! 
    N*N is the dimension of the fine grid
    (N+2)*(N+2) is the dimension of the vector with ghost layer 
    
    M*M is the dimension of the coarse grid
    (M+2)^2 is the dimension of the vector
*/
void restriction(double *fine_grid, int N, double *coarse_grid, int M){
    M=M+2;
    N=N+2;

    for(int i = 1; i<M-1; i++){
        for(int j = 1; j<M-1; j++){
            int k = 2*i;
            int l = 2*j;

            coarse_grid[i * M + j] = 0.125  * (fine_grid[(k+1) * N + l] + fine_grid[(k-1) * N + l] + fine_grid[k * N + (l-1)] + fine_grid[k * N + (l+1)]) + 0.5 * fine_grid[k * N + l];
        }
    }
}

// N x N is the dimension of the coarse grid
void prolongation(double *coarse_grid, int N, double* fine_grid, int M){
    N+=2;
    M+=2;
    // only iterate over the inner points and interpolate them from the coarse grid
    for(int i = 1; i<M-1; i++){
        for(int j=1; j<M-1; j++){
            int k = i/2;
            int l = j/2;
            if(i%2 == 0 && j%2==0){
                fine_grid[i * M + j] = coarse_grid[k * N + l];
            }
            else if(i%2 == 1 && j%2==0){
                fine_grid[i * M + j] = 0.5 * coarse_grid[k * N + l] + 0.5 * coarse_grid[(k+1) * N + l];

            }
            else if(i%2 == 0 && j%2==1){
                fine_grid[i * M + j] = 0.5 * coarse_grid[k * N + l] + 0.5 * coarse_grid[k * N + (l+1)];

            }
            else if(i%2 == 1 && j%2==1){
                fine_grid[i * M + j] =    0.25 * coarse_grid[k * N + l] 
                                        + 0.25 * coarse_grid[(k+1) * N + l]
                                        + 0.25 * coarse_grid[k * N + (l+1)]
                                        + 0.25 * coarse_grid[(k+1) * N + (l+1)];

            }
        }
    }
}

// Function to perform Jacobi smoothing
void smooth(double u[], double f[], int N, int v) {
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

void v_cycle(double** u, double **f, int N, int levels,int v){
    int vec_size_start = (N+2)*(N+2);
    int vec_size;
    
    int Nlevel= N;
    for (int l=levels-1;l>=1;l--){
        // number of points in the grid (without border)
        vec_size = (2+Nlevel)*(Nlevel+2);

        double *r;
        r = (double*)malloc((vec_size)*sizeof(double));
        // Smoothing
        smooth(u[l],f[l],Nlevel,v);

        // residual
        mfMult(Nlevel,u[l],r);
        axpy(r, -1, r, f[l],vec_size);

        // restriction
        int N_f = dim_coarser(Nlevel); 
        restriction(r,Nlevel, f[l-1],N_f);
        Nlevel=N_f;

        null_vec(u[l-1],Nlevel);

        free(r);
    }

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



void mg_solve(double** u, double **f, int N, int levels,int v){
    
    // setup
    int it_max = 49;
    int iterations = 0;

    double err;
    double epsilon=0.0001;

    int vec_size = (N+2)*(N+2);

    double* r;
    r=(double*)malloc(vec_size*sizeof(double));
    null_vec(r,vec_size);

    
    do {
        iterations += 1;

        // calculate new u with v cycle
        v_cycle(u,f,N,levels,v);

        // calculate new residual
        mfMult(N,u[levels-1],r);
        axpy(r,-1,r,f[levels-1],vec_size);

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

