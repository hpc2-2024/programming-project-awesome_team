#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <string.h>
#include "utils.h"

double function(double x, double y){
    return sin(y*M_PI) * sin(x*M_PI) * 2.0*M_PI*M_PI;
}
double function_sol(double x, double y){
    return sin(y*M_PI) * sin(x*M_PI); 
}
void solve_gs(int N, double *A, double *B, double *X, double tolerance, int max_iter){
    /*
    N is the number of all elements in the input variable *X to solve for
    */

    int M = sqrt(N);
    double h = 1.0/(M+1);
    double h2 = h*h;
    double epsilon = 0.000001;

    double *X_new;
    X_new = malloc(N * sizeof(*X_new));

    int k;
    for (k = 0; k < max_iter; k++) {
        for (int i = 0; i < N; i++) {
            double sum = 0;
            for (int j = 0; j < N; j++) {
                if (j != i) {
                    sum += A[i * N + j] * X[j];
                }
            }
            X_new[i] = (B[i] - sum) / (A[i * N + i]+epsilon);
        }

        double diff = abs_diff(N, X, X_new);

        for(int i = 0; i<N; i++){
            X[i] = X_new[i];
        }

        if(diff < tolerance){
            break;
        }
    }
    printf("Gauss-Seidel converged in %d iterations.\n", k);
    free(X_new);
}

void solve_mf_gs(int N, double *B, double *X, double tolerance, int max_iter){
    int M = sqrt(N);
    double h = 1.0/(M+1);
    double h2 = h*h;

    double *X_new;
    X_new = malloc(N * sizeof(*X_new));
    int k;
    for (k = 0; k < max_iter; k++) {
        for (int i = 0; i < N; i++) {
            double sum = 0;

            if(i==0){
                sum = -X[i+1] - X[i+M];
            }
            else if(1 <= i && i < M){
                sum = -X[i+1] - X[i+M] - X[i-1];
                // int deb;
            }
            else if(M <= i && i < N-M){
                sum = -X[i+1] - X[i+M] - X[i-1] - X[i-M];
                // int deb;
            }
            else if(N-M <= i && i < N-1){
                sum = -X[i+1] - X[i-M] - X[i-1];
                // int deb;
            }
            else if(i==N-1){
                sum = -X[i-1] - X[i-M];
                // int deb;
            }
            else{
                sum = 69;
            }
            sum = sum / h2;
            X_new[i] = (B[i] - sum) / (4 / h2);
        }

        double diff = abs_diff(N, X, X_new);

        for(int i = 0; i<N; i++){
            X[i] = X_new[i];
        }

        // check if the algorithm has converged
        if(diff < tolerance){
            break;
        }
    }
    printf("Matrix-Free Gauss-Seidel converged in %d iterations.\n", k);
    free(X_new);
}

// X1 is for regular Gauss-Seidel, X2 for matrix-free Gauss-Seidel
void debug_gs(int N, double *A, double *B, double *X1, double *X2, double tolerance, int max_iter){
    int M = sqrt(N);
    double h = 1.0/(M+1);
    double h2 = h*h;
    double epsilon = 0.000001;

    double *X1_new, *X2_new;
    X1_new = malloc(N * sizeof(*X1_new));
    X2_new = malloc(N * sizeof(*X2_new));

    for (int k = 0; k < max_iter; k++) {
        for (int i = 0; i < N; i++) {
            double sum1 = 0, sum2 = 0;

            for (int j = 0; j < N; j++) {
                if (j != i) {
                    sum1 += A[i * N + j] * X1[j];
                    int temp = 0;
                }
            }
            X1_new[i] = (B[i] - sum1) / (A[i * N + i]+epsilon);

            if(i==0){
                sum2 = -X2[i+1] - X2[i+M];
            }
            else if(1 <= i && i <M){
                sum2 = -X2[i+1] - X2[i+M] - X2[i-1];
                // int deb;
            }
            else if(M <= i && i < N-M){
                sum2 = -X2[i+1] - X2[i+M] - X2[i-1] - X2[i-M];
                // int deb;
            }
            else if(N-M <= i && i < N-1){
                sum2 = -X2[i+1] - X2[i-M] - X2[i-1];
                // int deb;
            }
            else if(i==N-1){
                sum2 = -X2[i-1] - X2[i-M];
                // int deb;
            }
            else{
                sum2 = 69;
            }
            sum2 = sum2 / h2;
            // printf("Iteration %d, Index %d, sum1=%f\n", k, i, sum1);
            // printf("Iteration %d, Index %d, sum2=%f\n\n", k, i, sum2);
            X2_new[i] = (B[i] - sum2) / (4 / h2);

            int deb = M + 1;

        }

        double diff1 = abs_diff(N, X1, X1_new) / abs_norm(N, X1);
        double diff2 = abs_diff(N, X2, X2_new) / abs_norm(N, X2);

        for(int i = 0; i<N; i++){
            X1[i] = X1_new[i];
            X2[i] = X2_new[i];
        }

        // print_matrix(M, X1);
        // print_matrix(M, X2);

        int deb = M + 1;

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

void test_GaussSeidel(){
    int N = 40;
    int N2 = N*N;
    double h = 1.0/(N+1);
    double h2 = h*h;

    double *X2, *X1, *B, *S, *A;
    X2 = malloc(N2 * sizeof(*X2));
    X1 = malloc(N2 * sizeof(*X1));
    B = malloc(N2 * sizeof(*B));
    S = malloc(N2 * sizeof(*S));
    A = malloc(N2 * N2 * sizeof(*A));

    get_stencil(N2, A);

    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            B[i * N + j] = function(i*h, j*h);
            S[i * N + j] = function_sol(i*h, j*h);
            X2[i * N + j] = 0;
            X1[i * N + j] = 0;
        }
    }

    int size_slice = 5;

    // print_matrix(N, B);
    // print_matrix(N, S);
    // print_matrix(N, X);
    // print_matrix(N, X1);
    // print_matrix(N2, A);

    solve_mf_gs(N2, B, X2, 0, 100);
    solve_gs(N2, A, B, X1, 0, 100);

    double rel_diff_mf = abs_diff(N2, X2, S) / abs_norm(N2, S); 
    double rel_diff = abs_diff(N2, X1, S) / abs_norm(N2, S); 
    
    printf("Exact Solution: \n");
    print_matrix_slice(N, size_slice, S);

    printf("Relative Error Gauss-Seidel:%f \n", rel_diff);
    printf("Solution Gauss-Seidel: \n");
    print_matrix_slice(N, size_slice, X1);

    printf("Relative Error matrix-free Gauss-Seidel:%f \n", rel_diff_mf);
    printf("Solution matrix-free Gauss-Seidel: \n");
    print_matrix_slice(N, size_slice, X2);

    free(X2);
    free(X1);
    free(A);
    free(B);
    free(S);
}

int main(int argc, char **argv){
    test_GaussSeidel();

    return 0;
}