#ifndef POISSON_MAT_VEK
#define POISSON_MAT_VEK

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>


/*! Implementation of a matrix free multiplication with 5-star stencil
If A is the Laplace Matrix then we compute y=Ar
params:
    N: N*N is the size of the grid 
    r[]: right hand side vector with ghost layer (=> (N+2)x(N+2) vector )
    y[]: the solution vector of our matrix multiplication with ghost layer (=> (N+2)x(N+2) vector )

 */
void mfMult(int N, double r[], double y[]){
    #pragma omp parallel for
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
 * @param stencil9 Option for using the 9 point stencil in 2dim (value 0 or 1).
 *
 * @note Vectors `r` and `y` should be pre-allocated:
 *       - For 2D, size is (N+2)*(N+2).
 *       - For 1D, size is N+2.
 */
void poisson_mat_vek(int dim, int N, double r[], double y[], int stencil9){
    double h = 1.0/(N+1);
    double h2 = h*h;
    int N_pad = N + 2;
    int N_pad2 = N_pad*N_pad;
    if (dim==3) {
        #pragma omp parallel for
        for (int i=1;i<N+1;i++){
            for (int j=1;j<N+1;j++){
                for (int k=1; k<N+1; k++) {
                    y[N_pad2*i+j*N_pad+k] = pow(1.0/h,2)*(
                        - r[N_pad2*(i-1)+j*N_pad+k] -r[N_pad2*(i+1)+j*N_pad+k]- r[N_pad2*i+(j-1)*N_pad+k]
                        + 6*r[N_pad2*i+j*N_pad+k] 
                        - r[N_pad2*i+(j+1)*N_pad+k] -r[N_pad2*i+j*N_pad+(k+1)] - r[N_pad2*i+j*N_pad+(k-1)]);
                }
               
            }
        }
    }
    else if (dim==2){
        #pragma omp parallel for
        for (int i=1;i<N+1;i++){
            for (int j=1;j<N+1;j++){
                if (stencil9 == 1) {
                    y[(N+2)*i+j] = (1.0/6.0)* pow(1.0/(h),2) * (-0.5*r[(N+2)*i+j-1] + 3 *r[(N+2)*i+j] - 0.5 * r[(N+2)*i+j+1] - 0.5 *  r[(N+2)*(i+1)+j] -0.25 * r[N_pad*(i+1)+j+1] -0.25 * r[N_pad*(i+1)+j-1] - 0.5 * r[(N+2)*(i-1)+j] -0.25* r[N_pad*(i-1)+j+1] - 0.25 * r[N_pad*(i-1)+j-1]);
                }
                else {
                    y[(N+2)*i+j] = pow(1.0/h,2)*(4*r[(N+2)*i+j]-r[(N+2)*(i+1)+j] -r[(N+2)*i+j-1]-r[(N+2)*(i-1)+j]-r[(N+2)*i+j+1]);
                }
            }
        }
    } else if (dim==1) {
        #pragma omp parallel for
        for (int i=1; i<N+1; i++){
            y[i]=pow(1.0/h,2)*(-r[i-1]+2*r[i]-r[i+1]);
        }
    }
}

#endif