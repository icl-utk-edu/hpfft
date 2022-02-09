/*
* ---------------
* FFTW backend
* ---------------
*/
#ifndef FIBER_BACKEND_FFTWPP_H
#define FIBER_BACKEND_FFTWPP_H

#include <stdio.h>

#if defined(FIBER_ENABLE_FFTWPP)
#include <fftw++.h>


/*!
 * \ingroup CPU_libraries
 * \addtogroup fiber_cpu Backend fftw++
 *
 * Wrappers and template specializations related to the FFTW backend.
 * Requires CMake option:
 * \code
 *  -D Heffte_ENABLE_FFTWPP=ON
 * \endcode
 * Flag list:
 * 
 */

//=================== Initialization (if required) ============================
int init_fftwpp(int option){
    fftw_mpi_init();
    return(0);
}

//=====================  Complex-to-Complex transform =========================

void compute_z2z_fftwpp( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, int *fftwpp_options, double *timer)
{
   
    // FFT plan
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

    // FFT Z2Z execution
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

}

//=====================  Real-to-Complex transform =========================

void compute_d2z_fftwpp( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, int *fftwpp_options, double *timer)
{

    // FFT plan
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

    // FFT R2Z execution
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

}

//=====================  Complex-to-Real transform =========================


void compute_z2d_fftwpp( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, int *fftwpp_options, double *timer)
{

    // FFT plan
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

    // FFT Z2R execution
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

}


#else

int init_fftwpp(int option)
{}


void compute_z2z_fftwpp( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, int *fftwpp_options, double *timer)
{}

void compute_d2z_fftwpp( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, int *fftwpp_options, double *timer)
{}

void compute_z2d_fftwpp( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, int *fftwpp_options, double *timer)
{}

#endif


#endif  //! FIBER_BACKEND_FFTWPP_H