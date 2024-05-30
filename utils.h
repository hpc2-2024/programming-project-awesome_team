#include <stdio.h>
#include <stdbool.h>

/*! Calculating the dot (scalar) product of 2 vectors */
double dot(double v[], double w[], int size) {
    double sum = 0;
    int i;

    #pragma omp parallel for private(i) reduction(+:sum) 
    for (i=0;i<size;i++){
        sum += v[i]*w[i];
    } 

    return sum;
}

/*! fills array with zeros */
void null_vec(double array[], int arrSize){
    for (int i=0;i<arrSize;i++){
        array[i]=0;
    }
}

/*! 
fills array with ghostlayer with random numbers in (0,1) intervall
the borders will be kept 0
*/
void rand_vec(double x[], int N){
    int seed = 123456;
    for (int i = 1;i<N+1;i++) {
        for (int j = 1;j<N+1;j++){
            // randomly initialize x with values in (0,1)
            srand(seed + i);
            double r = (double)rand() / (double)RAND_MAX;
            x[(N+2)*i+j]=r;
        }

    }
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