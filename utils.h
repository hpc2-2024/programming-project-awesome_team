#ifndef UTILS
#define UTILS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int dim_finer(int N){
    return (N-2)*2 + 2;
}

int dim_coarser(int N){
    return (N-2)/2 + 2;
}

void print_matrix(int N, double *matrix){
    for(int i = 0; i<N; i++){
        for(int j=0; j<N; j++){
            printf("%f    ", matrix[i * N + j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_matrix_slice(int N, int slice, double *matrix){
    for(int i = 0; i<slice; i++){
        for(int j=0; j<slice; j++){
            printf("%f    ", matrix[i * N + j]);
        }
        printf("\n");
    }
    printf("\n");
}


// N x N is the dimension of the coarse grid
void prolongation(int N, double *fine_grid, double *coarse_grid){
    int M = (N-2)*2 + 2;
    double boundary_val = 0;

    // explicitly set the boundary points
    for(int i = 0; i<M; i++){
        fine_grid[i * M + 0] = boundary_val;
        fine_grid[0 * M + i] = boundary_val;
        fine_grid[i * M + (M-1)] = boundary_val;
        fine_grid[(M-1) * M + i] = boundary_val;
    }

    // only iterate over the inner points and interpolate them from the coarse grid
    for(int i = 1; i<M-1; i++){
        for(int j=1; j<M-1; j++){
            int k = i/2;
            int l = j/2;
            if(i%2 == 0 && j%2==0){
                fine_grid[i * M + j] = coarse_grid[k * N + l];
                int p = 0;
            }
            else if(i%2 == 1 && j%2==0){
                fine_grid[i * M + j] = 0.5 * coarse_grid[k * N + l] + 0.5 * coarse_grid[(k+1) * N + l];
                int p = 0;

            }
            else if(i%2 == 0 && j%2==1){
                fine_grid[i * M + j] = 0.5 * coarse_grid[k * N + l] + 0.5 * coarse_grid[k * N + (l+1)];
                int p = 0;

            }
            else if(i%2 == 1 && j%2==1){
                fine_grid[i * M + j] =    0.25 * coarse_grid[k * N + l] 
                                        + 0.25 * coarse_grid[(k+1) * N + l]
                                        + 0.25 * coarse_grid[k * N + (l+1)]
                                        + 0.25 * coarse_grid[(k+1) * N + (l+1)];
                int p = 0;

            }
        }
    }
}

// N x N is the dimension of the fine grid
void restriction(int N, double *fine_grid, double *coarse_grid){
    int M = (N-2)/2 + 2;
    double boundary_val = 0;

    // explicitly set the boundary points
    for(int i = 0; i<M; i++){
        coarse_grid[i * M + 0] = boundary_val;
        coarse_grid[0 * M + i] = boundary_val;
        coarse_grid[i * M + (M-1)] = boundary_val;
        coarse_grid[(M-1) * M + i] = boundary_val;
    }

    for(int i = 1; i<M-1; i++){
        for(int j = 1; j<M-1; j++){
            int k = 2*i;
            int l = 2*j;

            coarse_grid[i * M + j] = 0.125  * (fine_grid[(k+1) * N + l] + fine_grid[(k-1) * N + l] + fine_grid[k * N + (l-1)] + fine_grid[k * N + (l+1)]) + 0.5 * fine_grid[k * N + l];
            int p = 0;

            // if (i == 1 && j == 1){
            //     coarse_grid[i * M + j] = 0.25 * (fine_grid[(k+1) * N + l] + fine_grid[k * N + (l+1)]) + (1/2) + fine_grid[k * N + l]; 
            //     int p = 0;

            // }
            // else if (i == 1 && j == M-1){
            //     coarse_grid[i * M + j] =  0.25 * (fine_grid[(k+1) * N + l] + fine_grid[k * N + (l-1)]) + (1/2) * fine_grid[k * N + l];

            // }
            // else if (i == M-1 && j == 1){
            //     coarse_grid[i * M + j] =  0.25 * (fine_grid[(k-1) * N + l] + fine_grid[k * N + (l+1)]) + (1/2) * fine_grid[k * N + l];

            // }
            // else if (i == M-1 && j == M-1){
            //     coarse_grid[i * M + j] =  0.25 * (fine_grid[(k-1) * N + l] + fine_grid[k * N + (l-1)]) + (1/2) * fine_grid[k * N + l];

            // }
            // else if (i == 1){
            //     coarse_grid[i * M + j] = 0.167  * (fine_grid[(k+1) * N + l] + fine_grid[k * N + (l-1)] + fine_grid[k * N + (l+1)]) + (1/2) * fine_grid[k * N + l];
                
            // }
            // else if (j == 1){
            //     coarse_grid[i * M + j] = 0.167  * (fine_grid[(k+1) * N + l] + fine_grid[(k-1) * N + l] + fine_grid[k * N + (l+1)]) + (1/2) * fine_grid[k * N + l];
                
            // }
            // else if (i == M-1){
            //     coarse_grid[i * M + j] = 0.167  * (fine_grid[(k-1) * N + l] + fine_grid[k * N + (l-1)] + fine_grid[k * N + (l+1)]) + (1/2) * fine_grid[k * N + l];
                
            // }
            // else if (j == M-1){
            //     coarse_grid[i * M + j] = 0.167  * (fine_grid[(k+1) * N + l] + fine_grid[(k-1) * N + l] + fine_grid[k * N + (l-1)]) + (1/2) * fine_grid[k * N + l];
                
            // }
            // else{
            //     coarse_grid[i * M + j] = 0.125  * (fine_grid[(k+1) * N + l] + fine_grid[(k-1) * N + l] + fine_grid[k * N + (l-1)] + fine_grid[k * N + (l+1)]) + 0.5 * fine_grid[k * N + l];
            //     int p = 0;

            // }
            
        }
    }
}

double abs_diff(int N, double *a, double *b){
    double diff = 0;
    for(int i = 0; i<N; i++){
            diff += fabs(a[i] - b[i]);
    }
    return diff / N;
}

double abs_norm(int N, double *a){
    double sum = 0;
    for(int i = 0; i<N; i++){
            sum += fabs(a[i]);
    }
    return sum / N;
}

/*! Implementation of a matrix free multiplication with 5-star stencil*/
void mfMult(int N, double r[], double y[], double h){
    for (int i=1;i<N+1;i++){
        for (int j=1;j<N+1;j++){
            y[(N+2)*i+j]=4*r[(N+2)*i+j]-r[(N+2)*(i+1)+j] -r[(N+2)*i+j-1]-r[(N+2)*(i-1)+j]-r[(N+2)*i+j+1];
            y[(N+2)*i+j]=h*h*y[(N+2)*i+j];
        }
    }
}


#endif