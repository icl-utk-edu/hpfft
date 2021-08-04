/*
* ---------------
* FFTW backend
* ---------------
*/
#ifndef FIBER_BACKEND_FFTW_H
#define FIBER_BACKEND_FFTW_H

#include <stdio.h>
#include "fftw3.h"

//=====================  Complex-to-Complex transform =========================

void compute_z2z_fftw( int const inbox_low[3], int const inbox_high[3],
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

    printf("FFTW: Get Complex-to-Complex using the tester available in C. \n");
    MPI_Abort(comm, 1);
}

//=====================  Real-to-Complex transform =========================


void compute_d2z_fftw( int const inbox_low[3], int const inbox_high[3],
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

    printf("FFTW: Get Real-to-Complex using the tester available in C. \n");
    MPI_Abort(comm, 1);
}




#endif  //! FIBER_BACKEND_FFTW_H