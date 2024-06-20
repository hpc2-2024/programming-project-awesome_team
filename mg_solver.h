#ifndef MG_SOLVER
#define MG_SOLVER

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "src/vcycle.h"
#include "src/poisson_mat_vek.h"


void mg_solve(double** u, double **f, int N, int levels,int v, int dim){
    int iter_max = 200;
    int iter = 0;
    double err;
    double epsilon=0.0001;

    int vec_size = get_vec_size(N,dim,1);

    double *r;
    r = malloc(vec_size * sizeof(*r));
    null_vec(r, vec_size);

    int debug = 0;

    if (debug==1){
        printf("\n Debug Mg method \n");
        printf("Grid size N: %d\n",N);
        iter_max = 1;
        printf("Max v cycle iterations: %d\n",iter_max);
        
        vec_print(N,u[levels-1],"u_0");
        vec_print(N,f[levels-1],"f_0");
        vec_print(N,r,"r_0");
    }
    
    do {
        iter += 1;

        // Perform a V-cycle to update the solution
        v_cycle(u, f, N, levels, v, dim, debug);


        // clean up f
        int Nlevel = dim_coarser(N);
        for (int i=levels-2;i>=0;i--){
            int vec_size_level = get_vec_size(Nlevel,dim,1);
            null_vec(f[i], vec_size_level );
            Nlevel=dim_coarser(Nlevel);
        }

        // Update residual
        poisson_mat_vek(dim,N, u[levels - 1], r);                //r=Au
        axpy(r, -1, r, f[levels - 1], vec_size);    //r= f-r

        err = norm(r, vec_size);
        printf("Error: %f\n",err);
        if (iter>iter_max){
            printf("Multigrid solve stopped after %d iter without convergence.", iter);
            break;
        }
    } while (err > epsilon);
    
    free(r);
    // solution:
}

#endif