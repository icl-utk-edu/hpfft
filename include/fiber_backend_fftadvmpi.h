/*
* -----------------
* FFTADVMPI backend
* -----------------
*/
#ifndef HPFFT_BACKEND_FFTADVMPI_H
#define HPFFT_BACKEND_FFTADVMPI_H

#include <stdio.h>

#if defined(HPFFT_ENABLE_FFTADVMPI)

#include <fftw3-mpi.h>
#include <mpi.h>
#include <math.h>
#include <complex.h>

//=================== Initialization (if required) ============================
int init_fftadvmpi(int option){
    return(0);
}

#define min(x, y) (((x) < (y)) ? (x) : (y))

void decompose(int N, int M, int p, int *n, int *s)
{
    int q = N / M;
    int r = N % M;
    *n = q + (r > p);
    *s = q * p + min(r, p);
}


void subcomm(MPI_Comm comm , int ndims , MPI_Comm subcomms[ndims])
{
    MPI_Comm comm_cart;
    int nprocs , dims[ndims], periods[ndims], remdims[ndims];
    for (int i = 0; i < ndims; i++)
        { dims[i] = periods[i] = remdims[i] = 0; }
    
    MPI_Comm_size(comm , &nprocs);
    MPI_Dims_create(nprocs , ndims , dims);
    MPI_Cart_create(comm , ndims , dims , periods , 1, &comm_cart);

    for (int i = 0; i < ndims; i++) {
        remdims[i] = 1;
        MPI_Cart_sub(comm_cart , remdims , &subcomms[i]);
        remdims[i] = 0;
    }
    MPI_Comm_free(&comm_cart);
}



void subarray(MPI_Datatype datatype ,
            int ndims ,
            int sizes[ndims],
            int axis ,
            int nparts ,
            MPI_Datatype subarrays[nparts])
{
    int subsizes[ndims], substarts[ndims], n, s;
    for (int i = 0; i < ndims; i++){
        subsizes[i] = sizes[i];
        substarts[i] = 0;
    }        
    for (int p = 0; p < nparts; p++) {
        decompose(sizes[axis], nparts , p, &n, &s);
        subsizes[axis] = n; substarts[axis] = s;
        MPI_Type_create_subarray(
        ndims , sizes , subsizes , substarts ,
        MPI_ORDER_C , datatype , &subarrays[p]);
        MPI_Type_commit(&subarrays[p]);
    }
}


void exchange(MPI_Comm comm ,
                MPI_Datatype datatype ,
                int ndims ,
                int sizesA[ndims],
                void *arrayA ,
                int axisA ,
                int sizesB[ndims],
                void *arrayB ,
                int axisB)
{
    int nparts;
    MPI_Comm_size(comm , &nparts);
    MPI_Datatype subarraysA[nparts], subarraysB[nparts];
    subarray(datatype , ndims , sizesA , axisA , nparts , subarraysA);
    subarray(datatype , ndims , sizesB , axisB , nparts , subarraysB);
    int counts[nparts], displs[nparts];
    for (int p = 0; p < nparts; p++)
        { counts[p] = 1; displs[p] = 0; }
        MPI_Alltoallw(arrayA , counts , displs , subarraysA ,
        arrayB , counts , displs , subarraysB , comm);
        for (int p = 0; p < nparts; p++) {
            MPI_Type_free(&subarraysA[p]);
            MPI_Type_free(&subarraysB[p]);
        }
}


// Helper function to compute local sizes
static int lsz(int N, MPI_Comm comm)
{
    int size , rank , n, s;
    MPI_Comm_size(comm , &size);
    MPI_Comm_rank(comm , &rank);
    decompose(N, size , rank , &n, &s);
    return n;
}

// Helper macros
#define product(n) (n[0]*n[1]*n[2])
#define allocate(t,n) malloc(product(n)*sizeof(t))
#define deallocate(a) free(a)

//=====================  Complex-to-Complex transform =========================

void compute_z2z_fftadvmpi( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, int *fftadvmpi_options, double *timer)
{

    printf("Calling FFTADVMPI \n");
    MPI_Datatype T = MPI_C_DOUBLE_COMPLEX;

    MPI_Comm P[2];
    subcomm(MPI_COMM_WORLD , 2, P);

    int N[3];
    N[0] = fftadvmpi_options[1];
    N[1] = fftadvmpi_options[2];
    N[2] = fftadvmpi_options[3];

    int local_nx, local_ny, local_nz;
    local_nx = inbox_high[0]-inbox_low[0]+1;
    local_ny = inbox_high[1]-inbox_low[1]+1;
    local_nz = inbox_high[2]-inbox_low[2]+1;

    int embed_x[]= {0,local_nx}; // size of x horizontal pencil
    int embed_y[]= {0,local_ny}; // size of y horizontal pencil
    int embed_z[]= {0,local_nz}; // size of z horizontal pencil

    void *plan_1, *plan_2, *plan_3;
    plan_1 = fftw_plan_many_dft(1, &local_nx, local_ny*local_nz, NULL, embed_x, 1, local_nx, NULL, embed_x, 1, local_nx, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_2 = fftw_plan_many_dft(1, &local_ny, local_nx, NULL, embed_y, local_nx, 1, NULL, embed_y, local_nx, 1, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_3 = fftw_plan_many_dft(1, &local_nz, local_nx*local_ny, NULL, embed_z, local_nx*local_ny, 1, NULL, embed_z, local_nx*local_ny, 1, FFTW_FORWARD, FFTW_ESTIMATE);


    void *temp  = calloc(32, 2*sizeof(double));

    // fftw_execute_dft(plan_1, (fftw_complex *) in, (fftw_complex *) out);

    // for(int i=0;i<local_nz;++i)
    //     fftw_execute_dft(plan_2, (fftw_complex *) in + i*local_nx*local_ny , (fftw_complex *) out + i*local_nx*local_ny);


    fftw_execute_dft(plan_1, (fftw_complex *) in, (fftw_complex *) temp);

    for(int i=0;i<local_nz;++i)
        fftw_execute_dft(plan_2, (fftw_complex *) temp + i*local_nx*local_ny , (fftw_complex *) temp + i*local_nx*local_ny);

    
    // int sizes_pre_remap[3]  = {local_nz, local_ny, local_nx};
    // int sizes_post_remap[3] = {local_nx, local_ny, local_nz};
    
    int sizes_pre_remap[3]  = {2,4,4};
    int sizes_post_remap[3] = {4,2,4};
    exchange(comm, T, 3, sizes_pre_remap , temp , 1, sizes_post_remap , out , 0);

    // int sizes_pre_remap[3]  = {4,4,2};
    // int sizes_post_remap[3] = {4,2,4};
    // // exchange(P[1], T, 3, sizes_pre_remap , temp , 2, sizes_post_remap , out , 1);
    // exchange(comm, T, 3, sizes_pre_remap , temp , 2, sizes_post_remap , out, 1);

    /*

    exchange(P[1], T, 3, sizesA , out , 2, sizesB , out , 1);

    for(int i=0;i<N[1];++i)
        fftw_execute_dft(plan_2, out + i*N[0]*N[2] , out + i*N[0]*N[2]);
    exchange(P[0], T, 3, sizesB , out , 1, sizesC , out , 0);

    fftw_execute_dft(plan_3, out, out);
    exchange(P[0], T, 3, sizesC , out , 0, sizesB , out , 1);

    fftw_destroy_plan(plan_1);
    fftw_destroy_plan(plan_2);
    fftw_destroy_plan(plan_3);

    MPI_Comm_free(&P[0]);
    MPI_Comm_free(&P[1]);
    */


}

//=====================  Real-to-Complex transform =========================

void compute_d2z_fftadvmpi( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, int *fftadvmpi_options, double *timer)
{

}

//=====================  Complex-to-Real transform =========================


void compute_z2d_fftadvmpi( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, int *fftadvmpi_options, double *timer)
{

  
}


#else
int init_fftadvmpi(int option)
{}

void compute_z2z_fftadvmpi( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, int *fftadvmpi_options, double *timer)
{
    printf("Calling FFTADVMPI ++++++++++++++++\n");


}

void compute_d2z_fftadvmpi( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, int *fftadvmpi_options, double *timer)
{}

void compute_z2d_fftadvmpi( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, int *fftadvmpi_options, double *timer)
{}

#endif


#endif  //! HPFFT_BACKEND_FFTADVMPI_H
