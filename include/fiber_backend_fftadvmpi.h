/*
* -----------------
* FFTADVMPI backend
* -----------------
*/
#ifndef FIBER_BACKEND_FFTADVMPI_H
#define FIBER_BACKEND_FFTADVMPI_H

#include <stdio.h>

#if defined(FIBER_ENABLE_FFTADVMPI)

#include <fftw3-mpi.h>
#include <mpi.h>
#include <math.h>
#include <complex.h>

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

    printf("\n");
    printf("\n");
    for(int i=0; i<32; i++)
        // printf(" %g+%gi \t ", in[i].r, in[i].i);
        printf(" %g \t ", (double* out)[i]);

    printf("\n");
    printf("\n");

    int N[3] = {4, 4, 4};

    MPI_Comm P[2];
    subcomm(comm, 2, P);

    int me, me_g;
    int num_ranks, num_ranks_g;

    MPI_Comm_rank(comm, &me_g);
    MPI_Comm_size(comm, &num_ranks_g);

    for (int i=0; i<2; i++)
    {
        MPI_Comm_rank(P[i], &me);
        MPI_Comm_size(P[i], &num_ranks);
    
        printf("---->>>> global[%d]/%d,  P[%d] ++++ %d/%d \n", me_g, num_ranks_g, i, me, num_ranks);
    }

    // Define elementary MPI datatype
    MPI_Datatype T = MPI_C_DOUBLE_COMPLEX;

    int sizesA[3] = {lsz(N[0],P[0]), lsz(N[1],P[1]), N[2]};
    int sizesB[3] = {lsz(N[0],P[0]), N[1], lsz(N[2],P[1])};
    int sizesC[3] = {N[0], lsz(N[1],P[0]), lsz(N[2],P[1])};
    
    printf("Input 0 ++++ %d ++++ %d \n", P[0], P[1]);
    printf("Input kia A ++++ %d ++++ %d ++++ %d \n", sizesA[0], sizesA[1], sizesA[2]);
    printf("Input kia B ++++ %d ++++ %d ++++ %d \n", sizesB[0], sizesB[1], sizesB[2]);
    printf("Input kia C ++++ %d ++++ %d ++++ %d \n", sizesC[0], sizesC[1], sizesC[2]);

    int inembed[]= {0,N[0]}; // size of 1 horizontal pencil
    int onembed[]= {0,N[0]};

    void *plan_1, *plan_2, *plan_3;
    plan_1 = fftw_plan_many_dft(1, &N[0], N[1]*2, NULL, inembed, 1, N[0],        NULL, onembed, 1, N[0], FFTW_FORWARD, FFTW_ESTIMATE);
    plan_2 = fftw_plan_many_dft(1, &N[1], N[1],   NULL, inembed, N[1], 1,        NULL, onembed, N[1], 1, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_3 = fftw_plan_many_dft(1, &N[2], N[1]*2, NULL, inembed, N[1]*N[2]/2, 1, NULL, onembed, N[1]*N[2]/2, 1, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute_dft(plan_1, in, out);

    // printf("\n");
    // printf("A-Alan kia [%d] \n", me);
    // printf("\n");
    // for(int i=0; i<32; i++)
    //     printf(" %g \t ", out[i]);
    //     // printf(" %g+%gi \t ", out[i].r, out[i].i);
    // printf("\n");
    // printf("\n");


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


#endif  //! FIBER_BACKEND_FFTADVMPI_H