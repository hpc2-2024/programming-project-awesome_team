#include <stdlib.h>
#include <stdio.h>
#include <math.h> 
#include "utils.h"

/* INPUT: A - a pointer to an array having dimension N x N 
 *        Tol - small tolerance number to detect failure when the matrix is near degenerate
 * OUTPUT: Apply Dolittle algorithm to matrix A, i.e. L has only 1's on the diagonal. We overwrite the values in A  with the values of L and U such that the output *         A = (L-I) + U 
 *         The permutation matrix is not stored as a matrix, but in an integer vector P of size N 
 */
void decompose_LUP(double *A, int N, double Tol, int *P) {
    int i, j, k, imax; 
    double maxA, *ptr, absA;

    for (int i = 0; i < N; i++)
        P[i] = i; //Unit permutation matrix, P[N] initialized with N

    for (i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k * N + i])) > maxA) { 
                maxA = absA;
                imax = k;
            }

        if (maxA < Tol){
            printf("LUP decomposition failed, matrix is degenerate. \n");
            exit(EXIT_FAILURE);
        }

        if (imax != i) {
            //pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;

            //pivoting rows of A
            for(int k = 0; k<N; k++){
                double val = A[i * N + k];
                A[i * N + k] = A[imax * N + k];
                A[imax * N + k] = val; 
            }
        }

        for (j = i + 1; j < N; j++) {
            A[j * N + i] /= A[i * N + i];

            for (k = i + 1; k < N; k++)
                A[j * N + k] -= A[j * N + i] * A[i * N + k];
        }
    }
}

/* INPUT: A,P filled in LUP_decompose; b - rhs vector; N - dimension
 * OUTPUT: x - solution vector of A*x=b
 */
void solve_LUP(double *A, int *P, double *b, int N, double *x) {

    for (int i = 0; i < N; i++) {
        x[i] = b[P[i]];

        for (int k = 0; k < i; k++)
            x[i] -= A[i * N + k] * x[k];
    }

    for (int i = N - 1; i >= 0; i--) {
        for (int k = i + 1; k < N; k++)
            x[i] -= A[i * N + k] * x[k];

        x[i] /= A[i * N + i];
    }
}

void test_LUP(){
    int N = 3;
    double A[] = {0, 5, 1, 4, 2, 1, 2, 7, 9};
    double b[3] = {1,2,3};
    double x[3];
    int P[3];

    double tol = 0.00001;

    decompose_LUP(A, N, tol, P);

    print_matrix(3, A);
    print_vector_int(3, P);

    solve_LUP(A, P, b, 3, x);

    print_vector_double(3, b);
    print_vector_double(3, x);

}

int main(){
    test_lup();

    return 0;
}