#ifndef RESTRICTION
#define RESTRICTION

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief Restricts a finer grid to a coarser grid using simple injection.
 *
 * Copies values from the fine grid to the coarse grid by taking every other point.
 *
 * @param fine_grid Input fine grid array.
 * @param M Number of internal points in the fine grid.
 * @param coarse_grid Output coarse grid array.
 * @param N Number of internal points in the coarse grid.
 * @param dim Dimension of the grid (1 or 2).
 *
 * @note The fine grid should be of size (M+2)^dim and the coarse grid should be of size (N+2)^dim.
 */
void restriction_simple(double *fine_grid, int M, double *coarse_grid, int N, int dim){
    int M_pad = M+2;
    int N_pad = N+2;

    #pragma omp parallel for
    for(int i = 0; i<N_pad; i++){
        for(int j = 0; j<N_pad; j++){
            int k = 2*i;
            int l = 2*j;

            coarse_grid[i * N_pad + j] = fine_grid[k * M_pad + l];

        }
    }
}

/**
 * @brief Restricts a finer grid to a coarser grid using half-weighted averaging.
 *
 * Computes the coarse grid values by weighted averaging of the fine grid values.
 *
 * @param fine_grid Input fine grid array.
 * @param M Number of internal points in the fine grid.
 * @param coarse_grid Output coarse grid array.
 * @param N Number of internal points in the coarse grid.
 * @param dim Dimension of the grid (1 or 2).
 *
 * @note The fine grid should be of size (M+2)^dim and the coarse grid should be of size (N+2)^dim.
 */
void restriction_half(double *fine_grid, int M, double *coarse_grid, int N,int dim){
    int M_pad = M+2;
    int N_pad = N+2;
    if (dim==2){
        #pragma omp parallel for
        for(int i = 1; i<N_pad-1; i++){
            for(int j = 1; j<N_pad-1; j++){
                int k = 2*i;
                int l = 2*j;

                coarse_grid[i * N_pad + j] = 0.125  * (
                    fine_grid[(k+1) * M_pad + l]
                    + fine_grid[(k-1) * M_pad + l]
                    + fine_grid[k * M_pad + (l-1)] 
                    + fine_grid[k * M_pad + (l+1)]) 
                    + 0.5 * fine_grid[k * M_pad + l];
            }
        }
    }
    else if (dim==1){
        #pragma omp parallel for
        for (int i=1;i<N_pad-1;i++){
            coarse_grid[i] = 0.25  * (fine_grid[2*i-1]+2*fine_grid[2*i]+fine_grid[2*i+1]);
        }
    }
}

int conv_3d(int N, int x, int y, int z){
    int N_pad = N+2;
    int N_pad2 = N_pad*N_pad;
    return N_pad2*x+N_pad*y+z;
}

void restriction_full(double *fine_grid, int M, double *coarse_grid, int N,int dim){
    int M_pad = M+2;
    int N_pad = N+2;

    if (dim==3){
        #pragma omp parallel for
        for(int i = 1; i<N_pad-1; i++){
            for(int j = 1; j<N_pad-1; j++){
                for (int k = 1; k<N_pad-1; k++){

                int a = 2*i;
                int b = 2*j;
                int c = 2*k;

                coarse_grid[i * N_pad + j] = 
                    ( 
                    2 * fine_grid[conv_3d(M, i-1,j-1,k-1)]
                    + 2 * fine_grid[conv_3d(M,i,j-1,k-1)]
                    + 2 * fine_grid[conv_3d(M,i+1,j-1,k-1)] 
                    + 2 * fine_grid[conv_3d(M,i-1,j,k-1)] 
                    + 4 * fine_grid[conv_3d(M,i,j,k-1)]
                    + fine_grid[conv_3d(M,i+1,j,k-1)]
                    + fine_grid[conv_3d(M,i-1,j+1,k-1)]
                    + fine_grid[conv_3d(M,i,j+1,k-1)]
                    + fine_grid[conv_3d(M,i+1,j+1,k-1)]

                    + 2 * fine_grid[conv_3d(M,i-1,j-1,k)]
                    + 2 * fine_grid[conv_3d(M,i,j-1,k)]
                    + 2 * fine_grid[conv_3d(M,i+1,j-1,k)] 
                    + 2 * fine_grid[conv_3d(M,i-1,j,k)] 
                    + 4 * fine_grid[conv_3d(M,i,j,k)]
                    + fine_grid[conv_3d(M,i+1,j,k)]
                    + fine_grid[conv_3d(M,i-1,j+1,k)]
                    + fine_grid[conv_3d(M,i,j+1,k)]
                    + fine_grid[conv_3d(M,i+1,j+1,k)]

                    + 2 * fine_grid[conv_3d(M,i-1,j-1,k+1)]
                    + 2 * fine_grid[conv_3d(M,i,j-1,k+1)]
                    + 2 * fine_grid[conv_3d(M,i+1,j-1,k+1)] 
                    + 2 * fine_grid[conv_3d(M,i-1,j,k+1)] 
                    + 4 * fine_grid[conv_3d(M,i,j,k+1)]
                    + fine_grid[conv_3d(M,i+1,j,k+1)]
                    + fine_grid[conv_3d(M,i-1,j+1,k+1)]
                    + fine_grid[conv_3d(M,i,j+1,k+1)]
                    + fine_grid[conv_3d(M,i+1,j+1,k+1)]
                    )/64.0;
                }

            }
        }
    }
    else if (dim==2) {
        #pragma omp parallel for
        for(int i = 1; i<N_pad-1; i++){
            for(int j = 1; j<N_pad-1; j++){
                int k = 2*i;
                int l = 2*j;

                coarse_grid[i * N_pad + j] = 
                    ( 2 * fine_grid[(k+1) * M_pad + l]
                    + 2 * fine_grid[(k-1) * M_pad + l]
                    + 2 * fine_grid[k * M_pad + (l-1)] 
                    + 2 * fine_grid[k * M_pad + (l+1)] 
                    + 4 * fine_grid[k * M_pad + l]
                    + fine_grid[(k+1)*M_pad + (l+1)]
                    + fine_grid[(k+1)*M_pad + (l-1)]
                    + fine_grid[(k-1)*M_pad + (l+1)]
                    + fine_grid[(k-1)*M_pad + (l-1)])/16.0;
            }
        }
    }
    else if (dim==1){ // same as half in 1d
        #pragma omp parallel for
        for (int i=1;i<N_pad-1;i++){
            coarse_grid[i] = 0.25  * (fine_grid[2*i-1]+2*fine_grid[2*i]+fine_grid[2*i+1]);
        }
    }
}

void restriction(double *fine_grid, int M, double *coarse_grid, int N,int dim){
    restriction_full(fine_grid,M, coarse_grid,N, dim);
}

#endif

