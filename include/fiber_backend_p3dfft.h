/*
* ---------------
* P3DFFT backend
* ---------------
*/
#ifndef FIBER_BACKEND_P3DFFT_H
#define FIBER_BACKEND_P3DFFT_H

#include <stdio.h>
#include "p3dfft.h"

//=====================  Complex-to-Complex transform =========================

void compute_z2z_p3dfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, double *timer)
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

    printf("P3DFFT: Get Complex-to-Complex using sample/C/test3d_c2c.c  \n");
    MPI_Abort(comm, 1);
}

//=====================  Real-to-Complex transform =========================


void compute_d2z_p3dfft( int const inbox_low[3], int const inbox_high[3],
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

    printf("P3DFFT: Get Real-to-Complex using sample/C/test3d_r2c.c  \n");
    MPI_Abort(comm, 1);
}




#endif  //! FIBER_BACKEND_P3DFFT_H