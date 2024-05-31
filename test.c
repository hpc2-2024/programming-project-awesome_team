#include <stdio.h>
#include <stdbool.h>
#include "utils.h"

void test_norm(){
    printf("Testing norm funtion - START \n");
    // TC1
    double x[4] = {1,1,1,1};
    if (norm(x,4)==2) {
        printf("TC 1 successful\n");
    } 

    printf("Test norm function END\n");
}

int main(){
    test_norm();

    return 0;
}