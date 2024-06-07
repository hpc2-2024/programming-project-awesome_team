#include <stdio.h>
#include <stdbool.h>
#include "utils.h"
#include "mg_solver.h"

void test_norm(){
    printf("Testing norm funtion - START \n");
    // TC1
    double x[4] = {1,1,1,1};
    if (norm(x,4)==2) {
        printf("TC 1 successful\n");
    } 

    printf("Test norm function END\n\n");
}

void test_prolongation(){
    printf("Testing  prolongation - START \n");
    double u_2h[] = {   0,0,0,
                        0,4,0,
                        0,0,0   };
    int N_2h=1;

    double u_h[5*5];
    null_vec(u_h,25);

    int N_h = 3;

    double expected_u[] = {
        0,0,0,0,0,
        0,1,2,1,0,
        0,2,4,2,0,
        0,1,2,1,0,
        0,0,0,0,0
    };

    prolongation(u_2h,N_2h,u_h,N_h);
    //vec_print(N_h,u_h,"u_h prolongation");

    int fail = 0;
    for (int i = 0;i<25;i++){
        if (u_h[i]!=expected_u[i]){
            fail = 1;
        }
    }
    if (fail==0){
        printf("TC1 successful\n");
    }
    else {
        printf("TC1 failed\n");
    }
    printf("Test prolongation - END\n\n");
    
}

void test_restriction(){
    printf("Testing  restriction - START \n");
    double expected_u[] = {
        0,0,0,
        0,2.5,0,
        0,0,0   };
    int N_2h=1;

    double u_2h[9];
    null_vec(u_2h,9);

    int N_h = 3;
    double u_h[] = {
        0,0,0,0,0,
        0,1,1,1,0,
        0,1,4,1,0,
        0,1,1,1,0,
        0,0,0,0,0
    };

    restriction(u_h,N_h, u_2h,N_2h);

    int fail = 0;
    for (int i = 0;i<9;i++){
        if (u_2h[i]!=expected_u[i]){
            fail = 1;
        }
    }
    if (fail==0){
        printf("TC1 successful\n");
    }
    else {
        printf("TC1 failed\n");
    }
    printf("Test restriction - END\n\n");
    
}

void test_gaussian_elemination(){
    printf("Testing gaussian elemination - START \n");

    //TC1
    double A[] = {
        1,2,
        3,4
    };
    double u[2];
    double f[] = {5,5};

    gaussian_elimination(A,f,u,2);
    print_1dim(2,u,"u gaussian TC1");  

    double expected_u[]={-5,5};
    if (u[0]!=expected_u[0] || u[1]!=expected_u[1]){
        printf("TC1 failed\n");
    }
    else {
        printf("TC1 successful\n");
    }
    printf("Testing gaussian elemination - END \n\n");
}
int main(){
    test_norm();
    test_prolongation();
    test_restriction();
    test_gaussian_elemination();

    return 0;
}