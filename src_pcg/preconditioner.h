#ifndef PRECONDITIONER
#define PRECONDITIONER

void lapl_matrix(double a[][5], int N){

    for (int i=0;i<N*N;i++) {
        //diagonal
        a[i][2]=4;

        if (i%N==0){
            a[i][1]=0;
        }
        else {
            a[i][1]=-1;
        }
        if (i%N==N-1){
            a[i][3]=0;
        }
        else {
            a[i][3]=-1;
        }
    }
    //outer -1 diagonals
    for (int i=N;i<N*N;i++){
        a[i][0]=-1;
    }
    for (int i=0;i<N*N-N;i++){
        a[i][4]=-1;
    }
    for (int i=0;i<N;i++){
        a[i][0]=0;
    }
    for (int i=N*N-N;i<N*N;i++){
        a[i][4]=0;
    }
}

/*!iterative ilu for laplace matrix*/
void ilu(double a[][5], int N,double epsilon,int max_it){
    int iteration_count = 0;
    while (iteration_count < max_it){
        iteration_count += 1;
        
        #pragma omp parallel for
        for (int i = 0;i<N*N;i++){
            if (i-N >= 0){
                a[i][0]=-1/a[i-N][2];
            }
            // we have to skip some li's on the second diagonal
            if (i%N>0) {
                a[i][1]=-1/a[i-1][2];
            }
            a[i][2]=4+a[i][0]+a[i][1];
        }

    }
}

//double_index_with_ghostlayer_to_single_index
int conv_idx(int i, int j, int N){
    return (i-1)*N+j-1;
}

/*! forward solve for 5 star stencil with ghostlayer in vectors*/
void forward_solve(double a[][5],double y[], double b[],int N){
    for (int i=1;i<N+1;i++){
        for (int j=1;j<N+1;j++){
            y[i*(N+2)+j]=b[i*(N+2)+j]-y[(i-1)*(N+2)+j]*a[conv_idx(i,j,N)][0]-y[i*(N+2)+j-1]*a[conv_idx(i,j,N)][1];
        }
    }
}

/*! backward solve for 5 star stencil with ghostlayer in vectors*/
void backward_solve(double a[][5],double x[],double y[],int N){
    for (int i=N;i>0;i--){
        for (int j=N;j>0;j--){
                      x[i*(N+2)+j]=(y[i*(N+2)+j]-a[conv_idx(i,j,N)][3]*x[i*(N+2)+j+1]-a[conv_idx(i,j,N)][4]*x[(i+1)*(N+2)+j])/a[conv_idx(i,j,N)][2];
// TB: previously was x[i*(N+2)+j]=(y[i*(N+2)+j]-a[conv_idx(i,j,N)][3]*x[i*(N+2)+j]-a[conv_idx(i,j,N)][4]*x[(i+1)*(N+2)+j])/a[conv_idx(i,j,N)][2]; there was just and index wrong
        }
    }
}

/*!
Jacobi preconditioner 
*/
void inv_diag(double y[], double b[], int N){
    for (int i=1;i<N+1;i++){
        for (int j=1;j<N+1;j++){
            y[i*(N+2)+j]=b[i*(N+2)+j]/(4);
        }
    }
}

/*!

!!! Probably some index error in here !!!

Gaus-Seidel preconditioner
*/
void gs_precon(double z[], double lplusd[][5], double r[], int N){
    for (int i=1;i<N+1;i++){
        for (int j=1;j<N+1;j++){
            z[i*(N+2)+j]=(r[i*(N+2)+j]-lplusd[(i-1)*N+j-1][0]*z[(i-1)*(N+2)+j]-lplusd[conv_idx(i,j,N)][1]*z[i*(N+2)+j-1])/4;
        }
    }
}

/*! 
ILU(0) preconditioner
params: 
        a[][5] lu decomposition of matrix of poisson equation
        r[] vector we solve for
        temp[] temporary vector we use as Ux (Ax = LUx = Ly =  b)
        z[] right hand side 
        N is the grid size
solution is overwriten in z
*/
void apply_precon(double a[][5],double r[], double z[],int N, int preconditioner){
    if (preconditioner==1) {

        int vec_size_ghost = (N+2)*(N+2);
        double* temp = (double *)malloc(vec_size_ghost*sizeof(double));
        null_vec(temp,vec_size_ghost);
        
        forward_solve(a,temp,r,N); // L* temp = r   
        backward_solve(a,z,temp,N); // Uz = temp

        free(temp);
    }
    else if (preconditioner==2){
        inv_diag(z,r,N);
    }
    else if (preconditioner==3){
        gs_precon(z,a,r,N);
    }
}

void init_preconditioner(double a[][5], double r[], double z[], int N, int preconditioner){
    if (preconditioner == 1){
        lapl_matrix(a,N);
        ilu(a,N,0.001,100);
        apply_precon(a,r,z,N,preconditioner);
    }
    else if (preconditioner == 2){
        inv_diag(z,r,N);
    }
    else if (preconditioner == 3){
        lapl_matrix(a,N);
        gs_precon(z,a,r,N);
    }
}

#endif