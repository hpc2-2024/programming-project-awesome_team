#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <stdbool.h>
#include "utils.h"

/*! The function f of the exercise sheet*/
double fun(double x, double y){
    return sin(y*M_PI)*sin(x*M_PI)*2.0*M_PI*M_PI;
}

double fun_solution(double x, double y){
    return sin(x*M_PI)*sin(y*M_PI);
}

void init_b(int N, double *b){
    int N_inner = N-2; // to compute h, we need the number of inner points 
    double h = 1.0/(N_inner+1);

    for (int i = 0; i<N; i++) {
        for (int j = 0; j<N; j++){
            if(i==0 || i==N-1 || j==0 || j==N-1){
                b[i * N + j] = 0.0;
            }
            else{
            b[i * N +j]=fun(i*h,j*h)*h*h; // TB: it is better to shift the h on the rhs, so you do not have a matirx that scale with h
            }
        }
    }
}

void init_s(int N, double *s){
    int N_inner = N-2; // to compute h, we need the number of inner points 
    double h = 1.0/(N_inner+1);

    for (int i = 0; i<N; i++) {
        for (int j = 0; j<N; j++){
            if(i==0 || i==N-1 || j==0 || j==N-1){
                s[i * N + j] = 0.0;
            }
            else{
            s[i * N +j]=fun_solution(i*h,j*h)*h*h; // TB: it is better to shift the h on the rhs, so you do not have a matirx that scale with h
            }
        }
    }
}

void print_matrix(int N , double *matrix){
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            printf("%f      ", matrix[i*N+j]);
        }
        printf("\n");
    }
    printf("\n");
}

void fill_val(int N, double *matrix, double val){
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            matrix[i * N + j] = val;
        }
    }
}

void fill_zeros(int N, double *matrix){
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            matrix[i * N + j] = 0.0;
        }
    }
}


// expects a padded coarse grid of size (N+2) x (N+2) and an empty, padded fine grid of size (M+2) x (M+2)
void prolongation(double *coarse_grid, int N, double* fine_grid, int M){
    int N_pad = N + 2;
    int M_pad = M + 2;

    // only iterate over the inner points and interpolate them from the coarse grid
    for(int i = 0; i<M_pad; i++){
        for(int j = 0; j<M_pad; j++){
            int k = (int) (i+1)/2;
            int l = (int) (j+1)/2;

            if(i == 0 || i == M_pad-1 || j == 0 || j == M_pad-1){
                fine_grid[i * M_pad + j] = 0.0; // set the boundary points to 0
            }
            else if(i%2 == 1 && j%2==1){
                fine_grid[i * M_pad + j] = coarse_grid[k * N_pad + l];
            }
            else if(i%2 == 0 && j%2==1){
                fine_grid[i * M_pad + j] = 0.5 * coarse_grid[k * N_pad + l] + 0.5 * coarse_grid[(k+1) * N_pad + l];

            }
            else if(i%2 == 1 && j%2==0){
                fine_grid[i * M_pad + j] = 0.5 * coarse_grid[k * N_pad + l] + 0.5 * coarse_grid[k * N_pad + (l+1)];

            }
            else if(i%2 == 0 && j%2==0){
                fine_grid[i * M_pad + j]=  0.25 * coarse_grid[k * N_pad + l] 
                                        + 0.25 * coarse_grid[(k+1) * N_pad + l]
                                        + 0.25 * coarse_grid[k * N_pad + (l+1)]
                                        + 0.25 * coarse_grid[(k+1) * N_pad + (l+1)];
            }
        }
    }
}

void test_prolongation(){
    int N = 3;
    int N_pad = N + 2;

    int M = 2 * N - 1; // 2 * N - 1 + 2 to account for padding 
    int M_pad = M + 2;

    double b_3[] =  {0, 0, 0, 0, 0,
                     0, 1, 2, 1, 0,
                     0, 1, 3, 1, 0,
                     0, 1, 2, 2, 0,
                     0, 0, 0, 0, 0};

    print_matrix(5, b_3);

    double *b_5;
    b_5 = malloc(M_pad * M_pad * sizeof(*b_5));

    prolongation(b_3, N, b_5, M);
    print_matrix(M_pad, b_5);
}

void restriction_simple(double *fine_grid, int M, double *coarse_grid, int N){
    int M_pad = M+2;
    int N_pad = N+2;

    for(int i = 0; i<N_pad; i++){
        for(int j = 0; j<N_pad; j++){
            int k = 2*i-1;
            int l = 2*j-1;

            // set boundary points to 0
            if(i == 0 || i == N_pad-1 || j == 0 || j == N_pad-1){
                coarse_grid[i * N_pad + j] = 0.0; // set the boundary points to 0
            } 
            else{
                coarse_grid[i * N_pad + j] = fine_grid[k * M_pad + l];
            }

        }
    }
}

void restriction_half(double *fine_grid, int M, double *coarse_grid, int N){
    int M_pad = M+2;
    int N_pad = N+2;

    for(int i = 0; i<N_pad; i++){
        for(int j = 0; j<N_pad; j++){
            int k = 2*i-1;
            int l = 2*j-1;

            // set boundary points to 0
            if(i == 0 || i == N_pad-1 || j == 0 || j == N_pad-1){
                coarse_grid[i * N_pad + j] = 0.0; // set the boundary points to 0
            } 
            else{
                coarse_grid[i * N_pad + j] = 0.125  * (fine_grid[(k+1) * M_pad + l] + fine_grid[(k-1) * M_pad + l] + fine_grid[k * M_pad + (l-1)] + fine_grid[k * M_pad + (l+1)]) + 0.5 * fine_grid[k * M_pad + l];
            }
        }
    }
}

void test_restriction(){
    int N = 3;
    int  M = 2 * N - 1;

    int N_pad = N + 2;
    int M_pad = M + 2;

    double b_M[] = {0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      
                    0.000000,      1.000000,      1.500000,      2.000000,      1.500000,      1.000000,      0.000000,      
                    0.000000,      1.000000,      1.750000,      2.500000,      1.750000,      1.000000,      0.000000,      
                    0.000000,      1.000000,      2.000000,      3.000000,      2.000000,      1.000000,      0.000000,      
                    0.000000,      1.000000,      1.750000,      2.500000,      2.000000,      1.500000,      0.000000,      
                    0.000000,      1.000000,      1.500000,      2.000000,      2.000000,      2.000000,      0.000000,      
                    0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000};

    print_matrix(M_pad, b_M);

    double *b_N, *b_O;
    b_N = malloc(N_pad * N_pad * sizeof(*b_N));
    fill_zeros(N_pad, b_N);

    restriction_simple(b_M, M, b_N, N);
    print_matrix(N_pad, b_N);

}

void test_prolongation_restriction(){
    int N = 3;
    int M = 2 * N - 1;
    int O = 2 * M - 1;

    int N_pad = N + 2;
    int M_pad = M + 2;
    int O_pad = O + 2;

    double b_N[] = {0, 0, 0, 0, 0,
                    0, 1, 2, 1, 0,
                    0, 1, 3, 1, 0,
                    0, 1, 2, 2, 0,
                    0, 0, 0, 0, 0};

    printf("Matrix N1:\n");
    print_matrix(N_pad, b_N);

    double *b_M, *b_O;
    b_M = malloc(M_pad * M_pad * sizeof(*b_M));
    b_O = malloc(O_pad * O_pad * sizeof(*b_O));
    // fill_zeros(M_pad, b_M);
    // fill_zeros(O_pad, b_O);

    prolongation(b_N, N, b_M, M);
    printf("Matrix M1:\n");
    print_matrix(M_pad, b_M);

    prolongation(b_M, M, b_O, O);
    printf("Matrix O2:\n");
    print_matrix(O_pad, b_O);

    // fill_zeros(M_pad, b_M);
    restriction_simple(b_O, O, b_M, M);
    printf("Matrix M2:\n");
    print_matrix(M_pad, b_M);

    // fill_zeros(N_pad, b_N);
    restriction_simple(b_M, M, b_N, N);
    printf("Matrix N2:\n");
    print_matrix(N_pad, b_N);
}

int main(){
    int N = 5; 
    int N_pad = N + 2;

    double *b, *s;
    b = malloc(N_pad * N_pad * sizeof(*b));
    s = malloc(N_pad * N_pad * sizeof(*s));

    init_b(N_pad, b);
    init_s(N_pad, s);

    printf("Right-Hand Side b: \n");
    print_matrix(N_pad, b);

    printf("Solution s: \n");
    print_matrix(N_pad, s);

    // test_prolongation_restriction();
    return 0;
}