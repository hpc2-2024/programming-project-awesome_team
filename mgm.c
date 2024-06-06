#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <string.h>

double abs_diff(int N, double *a, double *b){
    double diff = 0;
    for(int i = 0; i<N*N; i++){
            diff += fabs(a[i] - b[i]);
    }
    return diff / (N*N);
}

double abs_norm(int N, double *a){
    double sum = 0;
    for(int i = 0; i<N*N; i++){
            sum += fabs(a[i]);
    }
    return sum / (N*N);
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


double function(double x, double y){
    return sin(y*M_PI) * sin(x*M_PI) * 2.0*M_PI*M_PI;
}
double function_sol(double x, double y){
    return sin(y*M_PI) * sin(x*M_PI); 
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

void solve_gs(int N, double *a, double *b, double *x, double tolerance, int max_iter){
    int N2 = N*N;
    // double h = 1.0/(N+1);
    // double h2 = h*h;

    int i, j, k;
    double sum, *x_new;
    x_new = malloc(N  * sizeof(*x_new));
    double epsilon = 0.0000001;

    for (k = 0; k < max_iter; k++) {
        for (i = 0; i < N2; i++) {
            sum =  b[i];
            for (j = 0; j < N2; j++) {
                if (j != i) {
                    sum -= a[i * N2 + j] * x[j]; //probably problematic line
                }
            }
            x_new[i] = sum / (a[i * N2 + i]+epsilon);
        }

        // Check for convergence
        double diff = abs_diff(N, x, x_new) / abs_norm(N, x);
        for(int i = 0; i<N2; i++){
            x[i] = x_new[i];
        }

        // check if the algorithm has converged
        if(diff < tolerance){
            printf("Solution Gauss-Seidel: \n");
            // print_matrix(N, x);
            printf("Iterations needed: %d \n\n", k);
            
            free(x_new);
            break;
        }
    }

}

void solve_mf_gs(int N, double *B, double *X, double tolerance, int max_iter){
    int N2 = N*N;
    double h = 1.0/(N+1);
    double h2 = h*h;
    // double epsilon = 0.0000001;

    double *X_new;
    X_new = malloc(N2 * sizeof(*X_new));

    for(int k = 1; k<max_iter; k++){
        for(int i = 0; i<N2; i++){
            double sum = 0;
            if(i==0){
                sum = -X[i+1] - X[i+N];
                // int deb = 0 + sum;
            }
            else if(1 <= i <N){
                sum = -X[i+1] - X[i+N] - X[i-1];
                // int deb;
            }
            else if(N <= i < N2-N){
                sum = -X[i+1] - X[i+N] - X[i-1] - X[i-N];
                // int deb;
            }
            else if(N2-N <= i < N2-1){
                sum = -X[i+1] - X[i-N] - X[i-1];
                // int deb;
            }
            else if(i==N2-1){
                sum = -X[i-1] - X[i-N];
                // int deb;
            }
            else{
                sum = 69;
            }
            X_new[i] = (h2 * B[i] - sum) / 4 ;
        }

        // compute the difference of the old and new x
        // update old x
        double diff = abs_diff(N, X, X_new) / abs_norm(N, X);
        for(int i = 0; i<N2; i++){
            X[i] = X_new[i];
        }

        // check if the algorithm has converged
        if(diff < tolerance){
            free(X_new);

            printf("Solution MF-Gauss-Seidel: \n");
            print_matrix(N, X);
            printf("Iterations needed: %d \n\n", k);

            break;
        }
    }
}

void get_stencil(int N2, double *A){
    int N = sqrt(N2);
    double h = 1.0/(N+1);
    double h2 = h*h;

    for(int i = 0; i<N2; i++){
        for(int j = 0; j<N2;j++){
            if(i == j){
                A[i * N2 + j] = 4.0 / h2;
            }
            else if((i == j-1) || (j == i-1) || (i == j-N) || (j == i-N)){
            // else if(j == i - 1){
                A[i * N2 + j] = -1.0 / h2;
            }
            else{
                A[i * N2 + j] = 0.0;
            }
        }
    }
}

int main(int argc, char **argv){
    int N = 3;
    int N2 = N*N;
    double h = 1.0/(N+1);
    double h2 = h*h;

    double *X, *X1, *B, *S, *A;
    X = malloc(N2 * sizeof(*X));
    X1 = malloc(N2 * sizeof(*X1));
    B = malloc(N2 * sizeof(*B));
    S = malloc(N2 * sizeof(*S));
    A = malloc(N2 * N2 * sizeof(*A));

    get_stencil(N2, A);

    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            B[i * N + j] = function(i*h, j*h);
            S[i * N + j] = function_sol(i*h, j*h);
            X[i * N + j] = 0;
            X1[i * N + j] = 0;
        }
    }

    int size_slice = 10;

    // print_matrix(N, B);
    // print_matrix(N, S);
    // print_matrix(N, X);
    // print_matrix(N, X1);
    // print_matrix(N2, A);

    solve_gs(N, A, B, X, 0.0001, 100);

    // solve_mf_gs(N, B, X1, 0.0001, 1000);

    // printf("Exact Solution: \n");
    // print_matrix(N, S);

    double rel_diff = abs_diff(N, X, S) / abs_norm(N, S);
    printf("Relative Absolute Difference: %f \n", rel_diff);

    free(X);
    free(X1);
    free(A);
    free(B);
    free(S);

    return 0;

}