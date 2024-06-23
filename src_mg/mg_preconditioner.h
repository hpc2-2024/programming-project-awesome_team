#include "utils.h"
#include "vcycle.h"

int calc_number_of_levels(int N){
    // Calculate level
    int level = 1;
    while ((N-1) % 2 == 0 && N>1) {
        N = (N-1)/2;
        level++;
    }
    return level;
}

void mg_precon(double* z, double* r, int N_start){
    int mg_iter = 2;
    int dimension = 2;
    int N = N_start;
    int vec_size = get_vec_size(N, dimension,1);
    int v = 3; // smoothing steps
    int use_stencil9 = 0;
    int debug = 0;
    
    int levels = calc_number_of_levels(N);

    double** u = allocate_multigrid(N, levels, dimension);
    double** f = allocate_multigrid(N, levels, dimension);

    axpy(f[levels-1],1,r,0, vec_size); // f[levels -1] = r
    
    for (int i=0; i<mg_iter; i++) {
        //V_cycle
        v_cycle(u, f, N, levels, v, dimension, use_stencil9, debug);
    }

    axpy(z, 1, u[levels-1], 0, vec_size); // z = u[levels -1] 

}