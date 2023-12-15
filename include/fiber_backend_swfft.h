/*
* ---------------
* SWFFT backend
* ---------------
*/
#ifndef HPFFT_BACKEND_SWFFT_H
#define HPFFT_BACKEND_SWFFT_H

#include <stdio.h>

#if defined(HPFFT_ENABLE_SWFFT)
#include "hpfft_swfft.h"

#include <fftw3.h>

//=================== Initialization (if required) ============================
int init_swfft(int option){
    return(0);
}

//=====================  Complex-to-Complex transform =========================

void compute_z2z_swfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, int *swfft_options, double *timer) {

    // Plan creation
    // void *plan;
    void *dfft, *dist;
    int i, dim3d[3];

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();

    // plan create
    for (i = 0; i < 3; ++i)
        dim3d[i] = inbox_high[i] - inbox_low[i];

    SWFFT_Distribution_new(comm, dim3d, 0);
    dfft = SWFFT_Dfft_new(dist);

    SWFFT_makePlans(dfft, in, out, in, out, FFTW_MEASURE);

    MPI_Barrier(comm);
    timer[0] += MPI_Wtime();

    // FFT execution
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

    // compute FFT 

    SWFFT_forward(dfft, in);
    SWFFT_backward(dfft, in);

    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

    // Delete plan
    SWFFT_Dfft_delete(dfft);
    SWFFT_Distribution_delete(dist);

    printf("SWFFT: Get Complex-to-Complex using TestDfft.cpp or TestFDfft.f90  \n");
    MPI_Abort(comm, 1);
}

//=====================  Real-to-Complex transform =========================


void compute_d2z_swfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, int *swfft_options, double *timer)
{
    printf("Real-to-Complex transform is not available for SWFFT. \n");
    MPI_Abort(comm, 1);
    timer[0] = -1;
    timer[1] = -1;
}

//=====================  Complex-to-Real transform =========================

void compute_z2d_swfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, int *swfft_options, double *timer)
{
    printf("Complex-to-Real transform is not available for SWFFT. \n");
    MPI_Abort(comm, 1);
    timer[0] = -1;
    timer[1] = -1;
}

#else

int init_swfft(int option)
{}

void compute_z2z_swfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, int *swfft_options, double *timer)
{}

void compute_d2z_swfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, int *swfft_options, double *timer)
{}

void compute_z2d_swfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, int *swfft_options, double *timer)
{}

#endif

#endif  //! HPFFT_BACKEND_SWFFT_H
