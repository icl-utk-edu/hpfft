/*
* ---------------
* DECOM2D backend
* ---------------
*/
#ifndef FIBER_BACKEND_DECOM2D_H
#define FIBER_BACKEND_DECOM2D_H

#include <stdio.h>

#if defined(FIBER_ENABLE_2DECOMP)
#include <decomp_2d_iface.h>

//=================== Initialization (if required) ============================
int init_decomp2d(int physical, int nx, int ny, int nz, int p_row, int p_col){
    // plan create
    decomp_2d_init( nx, ny, nz*2, p_row, p_col );
    decomp_2d_fft_init( physical );
    return(0);
}

int finalize_decomp2d(){
    // Delete plan ...
    decomp_2d_finalize();
    return(0);
}

//=====================  Complex-to-Complex transform =========================
void compute_z2z_decomp2d( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, int *decomp2d_options, double *timer)
{
    // Plan creation ...
    int physical  = PHYSICAL_IN_X; // XXX has to be the same as the init phase
    int lxsize    = 0;
    int lysize    = 0;
    int lzsize    = 0;
    int direction = 0;

    direction = ( decomp2d_options[0] == 0 ? FORWARD : BACKWARD );

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();

    decomp_2d_get_local_sizes( physical, &lxsize, &lysize, &lzsize );

    MPI_Barrier(comm);
    timer[0] += MPI_Wtime();

    // FFT execution ...
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

    // compute FFT 
    decomp_2d_fft_3d_c2c( lxsize, lysize, lzsize, (double complex*)in, 
        lxsize, lysize, lzsize, (double complex*)out, direction );

    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

}

//=====================  Real-to-Complex transform =========================


void compute_d2z_decomp2d( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, int *decomp2d_options, double *timer)
{
    printf( "SKIP %s\n", __FUNCTION__ );
    return ;

    // Plan creation ...
    // void *plan;

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();

    // plan create

    MPI_Barrier(comm);
    timer[0] += MPI_Wtime();



    // FFT execution ...
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

    // compute FFT 

    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

    // Delete plan ...

    printf("DECOM2D: Get Real-to-Complex using a wrapper to functions in fft_test_r2c folder \n");
    MPI_Abort(comm, 1);
}

//=====================  Complex-to-Real transform =========================

void compute_z2d_decomp2d( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, int *decomp2d_options, double *timer)
{

// Missing!


}


#else
int init_decomp2d(int physical, int nx, int ny, int nz, int p_row, int p_col)
{}

int finalize_decomp2d()
{}

void compute_z2z_decomp2d( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, int *decomp2d_options, double *timer)
{}

void compute_d2z_decomp2d( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, int *decomp2d_options, double *timer)
{}

void compute_z2d_decomp2d( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, int *decomp2d_options, double *timer)
{}

#endif


#endif  //! FIBER_BACKEND_DECOM2D_H
