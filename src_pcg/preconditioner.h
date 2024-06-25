#ifndef PRECONDITIONER
#define PRECONDITIONER

#include "../src_mg/mg_preconditioner.h"

/*!
 * @brief Returns the Poisson matrix for a given size.
 *
 * @param a Input matrix to fill.
 * @param N Size of the matrix, i.e. a has shape (N x N)
 */
void lapl_matrix(double** a, int N){
    for (int i=0;i<N*N;i++) {
        //diagonal
        a[i][2]=4;

        if (i%N==0){
            a[i][1]=0;
        }
        else {
            a[i][1]=-1;
        }
        if (i%N==N-1){
            a[i][3]=0;
        }
        else {
            a[i][3]=-1;
        }
    }
    //outer -1 diagonals
    for (int i=N;i<N*N;i++){
        a[i][0]=-1;
    }
    for (int i=0;i<N*N-N;i++){
        a[i][4]=-1;
    }
    for (int i=0;i<N;i++){
        a[i][0]=0;
    }
    for (int i=N*N-N;i<N*N;i++){
        a[i][4]=0;
    }
}

/*!
 * @brief Compute the ILU decomposition for the Poisson matrix iteratively.
 *
 * @param a Input matrix.
 * @param N Size of the matrix, i.e. a has shape (N x N)
 * @param max_it Maximum number of iterations steps.
 */
/*!iterative ilu for laplace matrix*/
void ilu(double** a, int N,double epsilon,int max_it){
    int iteration_count = 0;
    while (iteration_count < max_it){
        iteration_count += 1;
        
        #pragma omp parallel for
        for (int i = 0;i<N*N;i++){
            if (i-N >= 0){
                a[i][0]=-1/a[i-N][2];
            }
            // we have to skip some li's on the second diagonal
            if (i%N>0) {
                a[i][1]=-1/a[i-1][2];
            }
            a[i][2]=4+a[i][0]+a[i][1];
        }
    }
}

/*!
 * @brief Given an index pair i and j of a matrix with ghost layers, compute the corresponding index in the flattened vector.
 *
 * @param i Index x-dimension.
 * @param y Index y-dimension.
 * @param N Size of the matrix. 
 */
//double_index_with_ghostlayer_to_single_index
int conv_idx(int i, int j, int N){
    return (i-1)*N+j-1;
}

/*!
 * @brief Perform the forward-solve for the 5-point stencil and a matrix with ghost layers.
 *
 * @param a Input matrix.
 * @param y Output vector.
 * @param b Right-hand side vector. 
 * @param N Number of inner points.
 */

/*! forward solve for 5 star stencil with ghostlayer in vectors*/
void forward_solve(double** a,double y[], double b[],int N){
    for (int i=1;i<N+1;i++){
        for (int j=1;j<N+1;j++){
            y[i*(N+2)+j]=b[i*(N+2)+j]-y[(i-1)*(N+2)+j]*a[conv_idx(i,j,N)][0]-y[i*(N+2)+j-1]*a[conv_idx(i,j,N)][1];
        }
    }
}

/*!
 * @brief Perform the backward-solve for the 5-point stencil and a matrix with ghost layers.
 *
 * @param a Input matrix.
 * @param x Output vector.
 * @param y Right-hand side vector. 
 * @param N Number of inner points.
 */

/*! backward solve for 5 star stencil with ghostlayer in vectors*/
void backward_solve(double** a,double x[],double y[],int N){
    for (int i=N;i>0;i--){
        for (int j=N;j>0;j--){
                      x[i*(N+2)+j]=(y[i*(N+2)+j]-a[conv_idx(i,j,N)][3]*x[i*(N+2)+j+1]-a[conv_idx(i,j,N)][4]*x[(i+1)*(N+2)+j])/a[conv_idx(i,j,N)][2];
// TB: previously was x[i*(N+2)+j]=(y[i*(N+2)+j]-a[conv_idx(i,j,N)][3]*x[i*(N+2)+j]-a[conv_idx(i,j,N)][4]*x[(i+1)*(N+2)+j])/a[conv_idx(i,j,N)][2]; there was just and index wrong
        }
    }
}

/*!
 * @brief Perform the Jacobi method matrix-free as a preconditioner.
 *
 * @param y Vector to fill.
 * @param b Right-hand side vector. 
 * @param N Number of inner points, i.e. the input is expected to have ghost layers.
 */
/*! Jacobi preconditioner */
void inv_diag(double y[], double b[], int N){
    for (int i=1;i<N+1;i++){
        for (int j=1;j<N+1;j++){
            y[i*(N+2)+j]=b[i*(N+2)+j]/(4);
        }
    }
}

/*!
 * @brief Perform the Gauss-Seidel method matrix-free as a preconditioner.
 */
/*! Gaus-Seidel preconditioner */
void gs_precon(double z[], double** lplusd, double r[], int N){
    for (int i=1;i<N+1;i++){
        for (int j=1;j<N+1;j++){
            z[i*(N+2)+j]=(r[i*(N+2)+j]-lplusd[(i-1)*N+j-1][0]*z[(i-1)*(N+2)+j]-lplusd[conv_idx(i,j,N)][1]*z[i*(N+2)+j-1])/4; // Probably some index error in here !!!

        }
    }
}

/*!
 * @brief Compute one of several preconditioners (ILU, Jacobi, Gauss-Seidel, Multi-Grid).
 *
 * @param a Pre-allocated matrix for ILU decomposition
 * @param r Residuum vector. 
 * @param z Vector for computing the backward solve.
 * @param preconditioner Specfiy which type of preconditioner to use: "1" = ILU, "2" = Jacobi, "3" = "Gauss-Seidel", "4" = Multi-Grid
 * @param smoother Specify which smoothing algorithm to use in the Multi-Grid method: "0" = Jacobi, "1" = Gauss-Seidel
 */
void apply_precon(double** a, double r[], double z[], int N, int preconditioner, int smoother){
    if (preconditioner==1) {
        //ILU(0) preconditioner
        int vec_size_ghost = (N+2)*(N+2);
        double* temp = (double *)malloc(vec_size_ghost*sizeof(double));
        null_vec(temp,vec_size_ghost);
        
        forward_solve(a,temp,r,N); // L* temp = r   
        backward_solve(a,z,temp,N); // Uz = temp

        free(temp);
    }
    else if (preconditioner==2){
        //Jacobi preconditioner
        inv_diag(z,r,N);
    }
    else if (preconditioner==3){
        gs_precon(z,a,r,N);
    }
    else if (preconditioner==4) {
        mg_precon(z, r, N, smoother);
    }
}

/*!
 * @brief Initialize one of several preconditioners (ILU, Jacobi, Gauss-Seidel, Multi-Grid).
 *
 * @param a Pre-allocated matrix for ILU decomposition
 * @param r Residuum vector. 
 * @param z Vector for computing the backward solve.
 * @param N Size, i.e. the vector we solve for has size (N * N).
 * @param preconditioner Specfiy which type of preconditioner to use: "1" = ILU, "2" = Jacobi, "3" = "Gauss-Seidel", "4" = Multi-Grid
 * @param smoother Specify which smoothing algorithm to use in the Multi-Grid method: "0" = Jacobi, "1" = Gauss-Seidel
 */
void init_preconditioner(double** a, double r[], double z[], int N, int preconditioner, int smoother){
    if (preconditioner == 1){
        lapl_matrix(a,N);
        ilu(a,N,0.001,100);
        apply_precon(a,r,z,N,preconditioner, smoother);
    }
    else if (preconditioner == 2){
        inv_diag(z,r,N);
    }
    else if (preconditioner == 3){
        lapl_matrix(a,N);
        gs_precon(z,a,r,N);
    }
    else if (preconditioner==4){
        mg_precon(z,r,N, smoother);
    }
}

#endif