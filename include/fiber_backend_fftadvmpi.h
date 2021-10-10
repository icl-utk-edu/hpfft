/*
* ---------------
* FFTADVMPI backend
* ---------------
*/
#ifndef FIBER_BACKEND_FFTADVMPI_H
#define FIBER_BACKEND_FFTADVMPI_H

#include <stdio.h>

#if defined(FIBER_ENABLE_FFTADVMPI)


//=====================  Complex-to-Complex transform =========================

void compute_z2z_fftadvmpi( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, int fftw_switch, double *timer)
{
    

}

//=====================  Real-to-Complex transform =========================

void compute_d2z_fftadvmpi( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, int fftw_switch, double *timer)
{

}

//=====================  Complex-to-Real transform =========================


void compute_z2d_fftadvmpi( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, int fftw_switch, double *timer)
{

  
}


#else

void compute_z2z_fftadvmpi( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, double *timer)
{}

void compute_d2z_fftadvmpi( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, double *timer)
{}

void compute_z2d_fftadvmpi( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, double *timer)
{}

#endif


#endif  //! FIBER_BACKEND_FFTADVMPI_H