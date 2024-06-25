#include "utils.h"
#include "vcycle.h"

/**
 * @brief Calculates the number of levels for a multigrid hierarchy.
 *
 * This function calculates the number of levels in a multigrid hierarchy given
 * the initial grid size `N`. It returns the maximal amount of levels possible.
 *
 * @param N The initial grid size.
 * @return The number of levels in the multigrid hierarchy.
 */
int calc_number_of_levels(int N){
    // Calculate level
    int level = 1;
    while ((N-1) % 2 == 0 && N>1) {
        N = (N-1)/2;
        level++;
    }
    return level;
}

/**
 * @brief Performs multigrid preconditioning.
 *
 * This function applies multigrid preconditioning to the input vector `r` and stores the result in `z`.
 *
 * @param[out] z Output vector to store the result of the preconditioning.
 * @param[in] r Input vector to which preconditioning is applied.
 * @param[in] N_start Initial grid size.
 */
void mg_precon(double* z, double* r, int N_start, int smoother){
    int mg_iter = 1;
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
        v_cycle(u, f, N, levels, v, dimension, use_stencil9, debug, smoother);
    }

    axpy(z, 1, u[levels-1], 0, vec_size); // z = u[levels -1] 

}