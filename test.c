#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "src_mg/utils.h"
#include "src_mg/mg_solver.h"
// #include "./src_mg/mg_solver.h"


/*! The function f of the exercise sheet*/
double fun(double x, double y){
    return sin(y*M_PI)*sin(x*M_PI)*2.0*M_PI*M_PI;
}

double fun_solution(double x, double y){
    return sin(x*M_PI)*sin(y*M_PI);
}

void init_b(double b[],int N){
    double h = 1.0/(N+1);
    // inner points of x_0,b
    for (int i = 1;i<N+1;i++) {
        for (int j = 1;j<N+1;j++){
            b[(N+2)*i+j]=fun(i*h,j*h)*h*h; // TB: it is better to shift the h on the rhs, so you do not have a matirx that scale with h
        }
    }
}
// void test_norm(){
//     printf("Testing norm funtion - START \n");
//     // TC1
//     double x[4] = {1,1,1,1};
//     if (norm(x,4)==2) {
//         printf("TC 1 successful\n");
//     } 

//     printf("Test norm function END\n\n");
// }

// void test_prolongation(){
//     printf("Testing  prolongation - START \n");
//     double u_2h[] = {   0,0,0,
//                         0,4,0,
//                         0,0,0   };
//     int N_2h=1;

//     double u_h[5*5];
//     null_vec(u_h,25);

//     int N_h = 3;

//     double expected_u[] = {
//         0,0,0,0,0,
//         0,1,2,1,0,
//         0,2,4,2,0,
//         0,1,2,1,0,
//         0,0,0,0,0
//     };

//     prolongation(u_2h,N_2h,u_h,N_h);
//     //vec_print(N_h,u_h,"u_h prolongation");

//     int fail = 0;
//     for (int i = 0;i<25;i++){
//         if (u_h[i]!=expected_u[i]){
//             fail = 1;
//         }
//     }
//     if (fail==0){
//         printf("TC1 successful\n");
//     }
//     else {
//         printf("TC1 failed\n");
//     }
//     printf("Test prolongation - END\n\n");
    
// }

// void test_restriction(){
//     printf("Testing  restriction - START \n");
//     double expected_u[] = {
//         0,0,0,
//         0,2.5,0,
//         0,0,0   };
//     int N_2h=1;

//     double u_2h[9];
//     null_vec(u_2h,9);

//     int N_h = 3;
//     double u_h[] = {
//         0,0,0,0,0,
//         0,1,1,1,0,
//         0,1,4,1,0,
//         0,1,1,1,0,
//         0,0,0,0,0
//     };

//     restriction(u_h,N_h, u_2h,N_2h);

//     int fail = 0;
//     for (int i = 0;i<9;i++){
//         if (u_2h[i]!=expected_u[i]){
//             fail = 1;
//         }
//     }
//     if (fail==0){
//         printf("TC1 successful\n");
//     }
//     else {
//         printf("TC1 failed\n");
//     }
//     printf("Test restriction - END\n\n");
    
// }

// void test_gaussian_elemination(){
//     printf("Testing gaussian elemination - START \n");

//     //TC1
//     double A[] = {
//         1,2,
//         3,4
//     };
//     double u[2];
//     double f[] = {5,5};

//     gaussian_elimination(A,f,u,2);
//     print_1dim(2,u,"u gaussian TC1");  

//     double expected_u[]={-5,5};
//     if (u[0]!=expected_u[0] || u[1]!=expected_u[1]){
//         printf("TC1 failed\n");
//     }
//     else {
//         printf("TC1 successful\n");
//     }
//     printf("Testing gaussian elemination - END \n\n");
// }

void test_jacobi_smoothing(){
    printf("Testing jacobi smoothing - START \n");

    printf("TC1 - Convergences with only Jacobi\n");
    int N = 3;
    int vec_size = (N+2)*(N+2);

    double u[vec_size];
    null_vec(u,vec_size);
    rand_vec(u,N,2);

    double b[vec_size];
    null_vec(b,vec_size);
    init_b(b,N);
    
    double r[vec_size];
    null_vec(r,vec_size);

    for (int i = 0; i<10000; i++){
        smooth_jacobi(u,b,N,1,2,0);
    }

    // calculate new residual
    poisson_mat_vek(2, N, u, r, 0);                //Au
    axpy(r, -1, r, b, vec_size);    //r= f-Au

    double err = norm(r,vec_size);
    printf("Error in jacobi: %f\n",err);


    printf("Testing jacobi smoothing - END \n\n");
}

void test_gauss_smoothing(){
    printf("Testing gauss-seidel smoothing - START \n");

    printf("TC1 - Convergences with only Gauss-Seidel\n");
    int N = 49;
    int vec_size = (N+2)*(N+2);

    double u[vec_size];
    null_vec(u,vec_size);
    rand_vec(u,N,2);

    double b[vec_size];
    null_vec(b,vec_size);
    init_b(b,N);
    
    double r[vec_size];
    null_vec(r,vec_size);

    for (int i = 0; i<10000; i++){
        smooth_gauss_seidel(u,b,N,1,2);
    }

    // calculate new residual
    mfMult(N, u, r);                //Au
    axpy(r, -1, r, b, vec_size);    //r= f-Au

    double err = norm(r,vec_size);
    printf("Error in Gauss-Seidel: %f\n",err);


    printf("Testing gauss-seidel smoothing - END \n\n");
}




int main(){
    // test_norm();
    // test_prolongation();
    // test_restriction();
    // test_gaussian_elemination();
    test_jacobi_smoothing();
    test_gauss_smoothing();
    return 0;
}