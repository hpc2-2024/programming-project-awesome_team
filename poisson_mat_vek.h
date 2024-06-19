#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>


/*! Implementation of a matrix free multiplication with 5-star stencil
If A is the Laplace Matrix then we compute y=Ar
params:
    N: N*N is the size of the grid 
    r[]: right hand side vector with ghost layer (=> (N+2)x(N+2) vector )
    y[]: the solution vector of our matrix multiplication with ghost layer (=> (N+2)x(N+2) vector )

 */
void mfMult(int N, double r[], double y[]){
    for (int i=1;i<N+1;i++){
        for (int j=1;j<N+1;j++){
            y[(N+2)*i+j]=4*r[(N+2)*i+j]-r[(N+2)*(i+1)+j] -r[(N+2)*i+j-1]-r[(N+2)*(i-1)+j]-r[(N+2)*i+j+1];
        }
    }
}

/*!
 * @brief Applies the Poisson matrix to a vector in 1D or 2D.
 *
 * Computes the result of applying the Poisson matrix to the input vector `r`
 * and stores the result in vector `y`. Supports both 1D and 2D based on `dim`.
 *
 * @param dim Dimension of the Poisson matrix (1 or 2).
 * @param N Number of internal grid points.
 * @param r Input vector.
 * @param y Output vector.
 *
 * @note Vectors `r` and `y` should be pre-allocated:
 *       - For 2D, size is (N+2)*(N+2).
 *       - For 1D, size is N+2.
 */
void poisson_mat_vek(int dim, int N, double r[], double y[]){
    double h = 1.0/(N+1);
    if (dim==2){
        for (int i=1;i<N+1;i++){
            for (int j=1;j<N+1;j++){
                y[(N+2)*i+j]=pow(1.0/h,2)*(4*r[(N+2)*i+j]-r[(N+2)*(i+1)+j] -r[(N+2)*i+j-1]-r[(N+2)*(i-1)+j]-r[(N+2)*i+j+1]);
            }
        }
    } else if (dim==1) {
        for (int i=1; i<N+1; i++){
            y[i]=pow(1.0/h,2)*(-r[i-1]+2*r[i]-r[i+1]);
        }
    }
}