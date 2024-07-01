#ifndef UTILS
#define UTILS

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>

/**
 * @brief Calculates the dot product of two vectors.
 * 
 * @param v First vector.
 * @param w Second vector.
 * @param size Size of the vectors.
 * @return double Dot product of v and w.
 */
double dot(double v[], double w[], int size) {
    double sum = 0;
    int i;

    #pragma omp parallel for private(i) reduction(+:sum) 
    for (i=0;i<size;i++){
        sum += v[i]*w[i];
    } 

    return sum;
}

/**
 * @brief Computes the Euclidean norm (L2 norm) of a vector.
 * 
 * @param arr Input vector.
 * @param arrSize Size of the vector.
 * @return double Euclidean norm of the vector.
 */
double norm(double arr[], int arrSize){
    double sol;
    sol = dot(arr,arr, arrSize);
    return sqrt(sol);
}

///// FUNCTIONS FOR FILLING VECTORS

/**
 * @brief Fills an array with a specified value.
 * 
 * @param N Size of the array.
 * @param array Pointer to the array to be filled.
 * @param val Value to fill the array with.
 */
void fill_val(int N, double array[], double val){
    for (int i=0; i<N; i++){
        array[i] = val;
    }
}

/**
 * @brief Initializes a vector to zero.
 * 
 * @param array Pointer to the array to be initialized.
 * @param arrSize Size of the array.
 */
void null_vec(double array[], int arrSize){
    for (int i=0;i<arrSize;i++){
        array[i]=0;
    }
}

/**
 * @brief Fills a vector with random numbers in the interval (0,1), preserving ghost layers.
 * 
 * @param x Pointer to the vector to be filled with random numbers (size (N+2)^dim.
 * @param N Size of the vector (excluding ghost layers).
 * @param dim Dimension of the vector (1 or 2).
 */
void rand_vec(double x[], int N, int dim){
    int N_pad = N+2;
    int seed = 123456;
    for (int i = 1;i<N+1;i++) {

        if (dim==3){
            for (int j = 1;j<N+1;j++){
                for (int k=1; k<N+1; k++) {
                    // randomly initialize x with values in (0,1)
                    srand(seed + i+j+k);
                    double r = (double)rand() / (double)RAND_MAX;
                    x[N_pad*N_pad*i+j*N_pad +k]=r;
                }
            }
        }
        else if (dim==2) {
            for (int j = 1;j<N+1;j++){
                // randomly initialize x with values in (0,1)
                srand(seed + i+j);
                double r = (double)rand() / (double)RAND_MAX;
                x[(N+2)*i+j]=r;
            }
        }
        else if (dim==1) {
            // randomly initialize x with values in (0,1)
            srand(seed + i);
            double r = (double)rand() / (double)RAND_MAX;
            x[i]=r;
        }

    }
}

/////////// Multigrid specific utils

/**
 * @brief Computes the number of inner points in the next finer grid level.
 * 
 * @param N Number of inner points in the current grid level.
 * @return int Number of inner points in the finer grid level.
 */
int dim_finer(int N){
    return N*2 + 1;
}

/**
 * @brief Computes the number of inner points in the next coarser grid level.
 * 
 * @param M Number of inner points in the current grid level.
 * @return int Number of inner points in the coarser grid level.
 */
int dim_coarser(int M){
    return (M-1)/2;
}

/**
 * @brief Computes the vector size for a grid with optional ghost layers.
 * 
 * @param N Number of inner points in one dimension.
 * @param dim Number of dimensions (1 or 2).
 * @param ghostlayer Boolean indicating if ghost layers should be included (0 or 1).
 * @return int Size of the vector.
 */
int get_vec_size(int N, int dim, int ghostlayer){
    if (ghostlayer!=0){
        N = N+2;
    }
    int vec_size = pow(N,dim);
    return vec_size;
}

/**
 * @brief Allocates memory for a multigrid hierarchy.
 * 
 * @param N Number of inner points in the finest grid.
 * @param levels Number of levels in the multigrid hierarchy.
 * @param dim Number of dimensions (1 or 2).
 * @return double** Pointer to the allocated multigrid array the array is filled with zeros.
 */
double** allocate_multigrid(int N, int levels, int dim) {
    double** grid = (double**)malloc(levels * sizeof(double*));
    int Nlevel = N;
    int vec_size_pad = get_vec_size(N, dim, 1);

    for (int i = levels - 1; i >= 0; i--){
        grid[i] = (double*)malloc( vec_size_pad * sizeof(double));
        null_vec(grid[i], vec_size_pad );

        Nlevel = dim_coarser(Nlevel);
        vec_size_pad = get_vec_size(Nlevel, dim, 1);
    }
    return grid;
}

/**
 * @brief Frees the memory allocated for a multigrid hierarchy.
 * 
 * @param grid Pointer to the multigrid array to be freed.
 * @param levels Number of levels in the multigrid hierarchy.
 */
void free_multigrid(double **grid, int levels) {
    for (int i = 0; i < levels; i++){
        free(grid[i]);
    }
    free(grid);
}

/////////////////////




///// FUNCTIONS FOR PRINTING VECTORS AND MATRICES /////////////

void print_matrix(int N , double *matrix){
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            printf("%f      ", matrix[i*N+j]);
        }
        printf("\n");
    }
    printf("\n");
}

/**
 * @brief Prints a vector with optional ghost layer.
 * 
 * @param arr Pointer to the vector to be printed (always 1d arrays).
 * @param N Size of the vector (excluding ghost layer if present).
 * @param dimension Dimension of the vector (1 for 1D, 2 for 2D).
 * @param ghostlayer Flag indicating if the vector includes a ghost layer (1 for true, 0 for false).
 * @param name Name or identifier of the vector for display purposes.
 */
void print_vector(double arr[], int N, int dimension, int ghostlayer, char name[]){
    int n = N;
    if (ghostlayer==1){
        n+=2;
    }
    printf("\n %s \n",name);
    for (int i=0; i<n; i++) {
        if (dimension == 2) {
            for (int j=0;j<n;j++){
                printf("%f ",arr[i*n+j]);
            }
        }
        else if (dimension == 1) {
            printf("%f \n", arr[i]);
        }
        else if (dimension == 3) {
            for (int j=0; j<n; j++) {
                for (int k=0; k<n; k++){
                    printf("%f ", arr[i*n*n+j*n+k]);
                }
                printf("\n");
            }
            printf("\n");
        }

        printf("\n");
    }
    printf("\n");
}

/*! Displaying a vector with ghostlayer*/
void vec_print(int N, double vec[], char name[]){
    printf("\n %s \n",name);
    for (int i=0;i<N+2;i++){
        for (int j=0;j<N+2;j++){
            printf("%f ",vec[(N+2)*i+j]);
        }
        printf("\n");
    }
    printf("\n");
}

/*! Even better than vec_print */
void amazing_vec_print(int N, double vec[], char name[]) {
    printf("\n %s \n",name);
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            printf("%f ",vec[(N)*i+j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_1dim(int N, double vec[],char name[]){
    printf("\n %s \n",name);
    for (int i=0;i<N;i++){
        printf("%f\n",vec[i]);
    }
}

void print_2dim(int N, double vec[][5],char name[]){
    printf("\n %s \n",name);
        for (int i=0;i<N*N;i++){
            for (int j=0;j<5;j++){
                printf("%f ",vec[i][j]);
        }
        printf("\n");
    }
}

/////////////////////////////

// axpy calculations with scalar y
void axpy_scalar_y(double solution[], double a, double x[], double y, int N) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        solution[i] = a * x[i] + ((y == 0) ? 0 : y);
    }
}

// axpy calculations with vector y
void axpy_vector_y(double solution[], double a, double x[], double y[], int N) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        solution[i] = a * x[i] + y[i];
    }
}

// overloaded axpy function
void axpy(double solution[], double a, double x[], double y[], int N) {
    if (y == NULL) { // If y is NULL, treat it as a scalar
        axpy_scalar_y(solution, a, x, 0, N);
    } else { // Otherwise, treat y as a vector
        axpy_vector_y(solution, a, x, y, N);
    }
}
// axpy calculations with vector y
void axpby(double solution[], double a, double x[],double b, double y[], int N) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        solution[i] = a * x[i] + b*y[i];
    }
}

#endif