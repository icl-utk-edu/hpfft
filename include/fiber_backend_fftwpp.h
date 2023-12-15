/*
* ---------------
* FFTW backend
* ---------------
*/
#ifndef HPFFT_BACKEND_FFTWPP_H
#define HPFFT_BACKEND_FFTWPP_H

#include <stdio.h>

#if defined(HPFFT_ENABLE_FFTWPP)
#include <cfftw++.h>


/*!
 * \ingroup CPU_libraries
 * \addtogroup hpfft_cpu Backend fftw++
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

    unsigned int nthreads = 1;
    set_fftwpp_maxthreads(nthreads);
    return(0);
}

//=====================  Complex-to-Complex transform =========================

void compute_z2z_fftwpp( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, int *fftwpp_options, double *timer)
{
    int nx, ny, nz;
    nx = fftwpp_options[1];
    ny = fftwpp_options[2];
    nz = fftwpp_options[3];

    int niter = fftwpp_options[5];

    // FFT plan
    void *plan_z2z;

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();
        printf("FFTW++ plan creation \n");
        plan_z2z = fftwpp_plan_3d(nx, ny, nz, comm);
    MPI_Barrier(comm);
    timer[0] = +MPI_Wtime();

    // FFT Z2Z execution
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();
	for(int i = 0; i < niter; i++){
        if (fftwpp_options[0]==0)
	    	fftwpp_execute(plan_z2z, in, out, -1);  // Execute forward
        if (fftwpp_options[0]==1)
	    	fftwpp_execute(plan_z2z, in, out, 1);   // Execute backward
	}
    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

	fftwpp_destroy_plan(plan_z2z);
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
{
    return(0);
}


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


#endif  //! HPFFT_BACKEND_FFTWPP_H
