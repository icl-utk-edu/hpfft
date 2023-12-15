/*
* ---------------
* FFTMPI backend
* ---------------
*/
#ifndef HPFFT_BACKEND_FFTMPI_H
#define HPFFT_BACKEND_FFTMPI_H

#include <stdio.h>

#if defined(HPFFT_ENABLE_FFTMPI)
#include "fft3d_wrap.h"



//=================== Initialization (if required) ============================
int init_fftmpi(int option){
    return(0);
}


//=====================  Complex-to-Complex transform =========================

void compute_z2z_fftmpi( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, int *fftmpi_options, double *timer)
{
    printf("Complex-to-Complex with FFTMPI. \n");
    // setup
    int nx,ny,nz;
    nx = ny = nz = 4;
    int precision = 2;

    int permute = 0;
    int fftsize; // FFT buffer size returned by FFT setup
    int sendsize,recvsize;


    // Plan definition
    void *plan;

    int cflag = 0;   // point/all/combo
    int eflag = 0;   // pencil/brick
    int pflag = 2;   // array/ptr/memcpy
    int rflag = 0;   // remap only, no 1d FFTs

    int tuneflag = 0;
    int tflag = 0;   // provide more timing details

    fft3d_create(comm, precision, &plan);
    fft3d_set(plan, "collective", cflag);
    fft3d_set(plan, "exchange", eflag);
    fft3d_set(plan, "pack", pflag);
    fft3d_set(plan, "remaponly", rflag);
    
    int scale = 0;
    fft3d_set(plan, "scale", scale);

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();

    if (!tuneflag) {
        fft3d_setup(plan, nx, ny, nz, 
                    inbox_low[0], inbox_high[0], inbox_low[1], inbox_high[1], inbox_low[2], inbox_high[2], 
                    outbox_low[0], outbox_high[0], outbox_low[1], outbox_high[1], outbox_low[2], outbox_high[2], 
                    permute, &fftsize, &sendsize, &recvsize);
    } else {
        printf("Tuning functions added.");
    }
    MPI_Barrier(comm);

    printf("given size by plan = %d \n", fftsize);

    timer[0] += MPI_Wtime();


    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

    if(fftmpi_options[0]==0)
        fft3d_compute(plan, (double *)in , (double *)in , -1); // -1 is for Forward
    else
        fft3d_compute(plan, (double *)in , (double *)in , 1); // 1 is for Backward

    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

    // Delete plan
    fft3d_destroy(plan);
}

//=====================  Real-to-Complex transform =========================


void compute_d2z_fftmpi( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, int *fftmpi_options, double *timer)
{
    printf("Real-to-Complex transform is not available for FFTMPI. \n");
    timer[0] = -1;
    timer[1] = -1;
    // MPI_Abort(comm, 1);
}

//=====================  Complex-to-Real transform =========================

void compute_z2d_fftmpi( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, int *fftmpi_options, double *timer)
{
    printf("Complex-to-Real transform is not available for FFTMPI. \n");
    timer[0] = -1;
    timer[1] = -1;
    // MPI_Abort(comm, 1);
}

#else
int init_fftmpi(int option)
{}

void compute_z2z_fftmpi( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, int *fftmpi_options, double *timer)
{}

void compute_d2z_fftmpi( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, int *fftmpi_options, double *timer)
{}

void compute_z2d_fftmpi( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, int *fftmpi_options, double *timer)
{}

#endif

#endif  //! HPFFT_BACKEND_FFTMPI_H
