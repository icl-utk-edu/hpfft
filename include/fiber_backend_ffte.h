/*
* ---------------
* FFTE backend
* ---------------
*/
#ifndef FIBER_BACKEND_FFTE_H
#define FIBER_BACKEND_FFTE_H

#include <stdio.h>
#include "ffte.h"

//=====================  Complex-to-Complex transform =========================

void compute_z2z_ffte( int const inbox_low[3], int const inbox_high[3],
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

    printf("Benchmarking FFTE: Get Complex-to-Complex backend using testmfft3d2v.c \n");
    MPI_Abort(comm, 1);
}

//=====================  Real-to-Complex transform =========================


void compute_d2z_ffte( int const inbox_low[3], int const inbox_high[3],
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


    printf("Benchmarking FFTE: Get Real-to-Complex backend using testmrfft3d2v.c  \n");
    MPI_Abort(comm, 1);
}



#endif  //! FIBER_BACKEND_FFTE_H