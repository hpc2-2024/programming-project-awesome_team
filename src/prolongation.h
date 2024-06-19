#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief Prolongates a coarser grid to a finer grid using simple interpolation.
 *
 * Interpolates values from the coarse grid to the fine grid. Supports both 1D and 2D based on `dim`.
 *
 * @param coarse_grid Input coarse grid array.
 * @param N Number of internal points in the coarse grid.
 * @param fine_grid Output fine grid array.
 * @param M Number of internal points in the fine grid.
 * @param dim Dimension of the grid (1 or 2).
 *
 * @note The coarse grid should be of size (N+2)^dim and the fine grid should be of size (M+2)^dim.
 */
void prolongation_simple(double *coarse_grid, int N, double* fine_grid, int M,int dim){
    int N_pad = N + 2;
    int M_pad = M + 2;

    if (dim==2) {
        // only iterate over the inner points and interpolate them from the coarse grid
        for(int i = 1; i<M_pad-1; i++){
            for(int j = 1; j<M_pad-1; j++){
                int k = (int) (i+1)/2;
                int l = (int) (j+1)/2;

                if(i%2 == 0 && j%2==0){
                    fine_grid[i * M_pad + j] = coarse_grid[k * N_pad + l];
                }
                else if(i%2 == 0 && j%2==1){
                    fine_grid[i * M_pad + j] = 0.5 * coarse_grid[k * N_pad + l] 
                                            + 0.5 * coarse_grid[(k-1) * N_pad + l];

                }
                else if(i%2 == 1 && j%2==0){
                    fine_grid[i * M_pad + j] = 0.5 * coarse_grid[k * N_pad + l] 
                                            + 0.5 * coarse_grid[k * N_pad + (l-1)];

                }
                else if(i%2 == 1 && j%2==1){
                    fine_grid[i * M_pad + j]=  0.25 * coarse_grid[k * N_pad + l] 
                                            + 0.25 * coarse_grid[(k-1) * N_pad + l]
                                            + 0.25 * coarse_grid[k * N_pad + (l-1)]
                                            + 0.25 * coarse_grid[(k-1) * N_pad + (l-1)];
                }
            }
        }
    }
    else if (dim==1) {
        // boundery condition
        fine_grid[0]=0;
        fine_grid[N_pad-1]=0;

        for (int i = 1; i<N_pad-1;i++){
            fine_grid[2*i]=coarse_grid[i];
        }

        for (int i=1;i<N_pad;i++){
            fine_grid[2*i-1] = 0.5 * (coarse_grid[i-1]+coarse_grid[i]);
        }
    }
}

