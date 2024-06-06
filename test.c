    
    // // test prolongation and restriction
    // for(int i = 1; i<N; i++){
    //     for(int j = 1; j<N; j++){
    //         X[i * N + j] = 0;
    //         B[i * N + j] = function(i*h, j*h);
    //         S[i * N + j] = function_sol(i*h, j*h);
    //     }
    // }

    // solve_gauss_seidel(N, B, X, 0.0005, 10000);

    // print_matrix(N, X);
    // print_matrix(N, S);

    // int dim_coarse = 5;
    // double coarse_grid[] = {0,0,0,0,0,
    //                         0,3,2,3,0,
    //                         0,4,3,5,0,
    //                         0,3,5,4,0,
    //                         0,0,0,0,0,};

    // int dim_fine = (dim_coarse-2)*2 + 2;

    // double *grid;
    // grid = malloc(N * sizeof(*grid));

    // for(int i; i<N; i++){
    //     for(int j; j<N; j++){
    //         grid[i * N + j] = 1;
    //     }
    // }
    // double coarse_grid[] = {1.0, 2.0,
                            // 2.0, 3.0};

    // double fine_grid[8*8];

    // print_matrix(dim_coarse, coarse_grid);

    // prolongation(dim_coarse, fine_grid, coarse_grid);

    // restriction(dim_fine, fine_grid, coarse_grid);

    // print_matrix(dim_coarse, coarse_grid);
