/*
* ---------------
* NB3DFFT backend
* ---------------
*/
#ifndef FIBER_BACKEND_NB3DFFT_H
#define FIBER_BACKEND_NB3DFFT_H

#include <stdio.h>
// #include "nb3dfft.h"

//=====================  Complex-to-Complex transform =========================

void compute_z2z_nb3dfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, double *timer)
{
    printf("Complex-to-Complex transform is not available for nb3dFFT. \n");
    MPI_Abort(comm, 1);
    timer[0] = -1;
    timer[1] = -1;
}

//=====================  Real-to-Complex transform =========================


void compute_d2z_nb3dfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, double *timer)
{

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

    printf("NB3DFFT: Get Real-to-Complex using a wrapper to functions in test_nb3dfft_simple.F90 \n");
    MPI_Abort(comm, 1);
}




#endif  //! FIBER_BACKEND_NB3DFFT_H