#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <string.h>

/*! The function f of the exercise sheet*/
double fun(double x, double y){
    return sin(y*M_PI)*sin(x*M_PI)*2.0*M_PI*M_PI;
}

double fun_solution(double x, double y){
    return sin(x*M_PI)*sin(y*M_PI);
}


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

/*! Implementation of a matrix free multiplication with 5-star stencil*/
void mfMult(int N, double r[], double y[], double h){
    for (int i=1;i<N+1;i++){
        for (int j=1;j<N+1;j++){
            y[(N+2)*i+j]=4*r[(N+2)*i+j]-r[(N+2)*(i+1)+j] -r[(N+2)*i+j-1]-r[(N+2)*(i-1)+j]-r[(N+2)*i+j+1];
            y[(N+2)*i+j]=(1/(h*h))*y[(N+2)*i+j];
        }
    }
}

/*! Displaying a vector */
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



int main(int argc, char** argv){
    // initilize variables
    int N = 30; // N^2 is the number of inner points in our lattice
    int N2 = (N+2)*(N+2);
    double epsilon = 0.0000001;
    double h = 1.0/(N+1);
    // init x0
    double alpha = 0;
    double beta=0;
    double err0;
    double *x,*p,*r,*b,*m;

    x=(double *)malloc(vec_size_ghost*sizeof(double));
    null_vec(x,vec_size_ghost);
    rand_vec(x,N); // random start vector x_0 

    p=(double *)malloc(vec_size_ghost*sizeof(double));
    null_vec(p,vec_size_ghost);

    r=(double *)malloc(vec_size_ghost*sizeof(double));
    null_vec(r,vec_size_ghost);

    b=(double *)malloc(vec_size_ghost*sizeof(double));
    null_vec(b,vec_size_ghost);

    m=(double *)malloc(vec_size_ghost*sizeof(double));
    null_vec(m,vec_size_ghost);

    // fill ghost layer with zeros (and everything else also 0)
    for (int i = 0;i<N+2;i++) {
        for (int j = 0;j<N+2;j++){
            x[(N+2)*i+j]=0;
            p[(N+2)*i+j]=0;
            r[(N+2)*i+j]=0;
            b[(N+2)*i+j]=0;
            m[(N+2)*i+j]=0;
        }

    }

    int seed = 123456;
    // inner points of x_0,b
    for (int i = 1;i<N+1;i++) {
        for (int j = 1;j<N+1;j++){
            // randomly initialize x with values in (0,1)
            srand(seed + i);
            double r = (double)rand() / (double)RAND_MAX;
            x[(N+2)*i+j]=r;

            b[(N+2)*i+j]=fun(i*h,j*h);
        }

    }

    // Calculate r_0
    mfMult(N,x,r,h);
    #pragma omp parallel for
    for (int i=0;i<(N+2)*(N+2);i++){
        r[i]-=b[i];
        p[i]=r[i]*(-1);
    }
    err0 = sqrt(dot(r,r,N2));
    
    int number_of_iterations = 0;
    do
    {
        number_of_iterations++;

        mfMult(N,p,m,h);
        double dot_rk_rk = dot(r,r,N2);
        alpha=dot_rk_rk/dot(p,m,N2);

        // update x,r
        #pragma omp parallel for
        for (int i=0;i<(N+2)*(N+2);i++){
            x[i]+=alpha*p[i];
            r[i]+=alpha*m[i];
        }
        //update p
        beta = dot(r,r,N2)/dot_rk_rk;
        #pragma omp parallel for
        for (int i=0;i<(N+2)*(N+2);i++){
            p[i]=-r[i]+beta*p[i];
        }

    } while (sqrt(dot(r,r,N2))/err0 >= epsilon);

    // Compute the relative absolute difference 
    double abs_diff = 0.0;
    double abs_sol = 0.0;
    for (int i = 1;i<N+1;i++) {
        for (int j = 1;j<N+1;j++){
            double x_sol = fun_solution(i*h,j*h);
            abs_diff += fabs(x[(N+2)*i+j] - x_sol);
            abs_sol += x_sol;
        }
    }
    printf("Run %d iterations.\n",number_of_iterations);
    printf("Relative absolute difference between exact and approx. solution: %f%%\n", 100.0 * abs_diff / abs_sol);


    // free allocated memory
    free(x);
    free(p);
    free(r);
    free(b);
    free(m);
    return 0;
}