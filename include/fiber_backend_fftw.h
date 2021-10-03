/*
* ---------------
* FFTW backend
* ---------------
*/
#ifndef FIBER_BACKEND_FFTW_H
#define FIBER_BACKEND_FFTW_H

#include <stdio.h>

#if defined(FIBER_ENABLE_FFTW)
#include <fftw3-mpi.h>

//=====================  Complex-to-Complex transform =========================

void compute_z2z_fftw( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, double *timer)
{
    int niter = 1;
    int nx = inbox_high[0]-inbox_low[0];
    int ny = inbox_high[1]-inbox_low[1];
    int nz = inbox_high[2]-inbox_low[2];

    // Plan creation ...
    void *plan_z2z;
    int tuning_option = 0

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();
    if (tuning_option==0)
        plan_z2z = fftw_mpi_plan_dft_c2c_3d(nx, ny, nz, in, out, comm, FFTW_ESTIMATE);
    if (tuning_option==1)
        plan_z2z = fftw_mpi_plan_dft_c2c_3d(nx, ny, nz, in, out, comm, FFTW_MEASURE);
    timer[0] += MPI_Wtime();

    // FFT execution
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

	// for(int i = 0; i < niter; i++){
		fftw_execute(plan_z2z);
	// }

    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

    // Delete plan
	fftw_destroy_plan(plan_z2z);
}

//=====================  Real-to-Complex transform =========================


void compute_d2z_fftw( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, double *timer)
{

    int niter = 1;
    int nx = inbox_high[0]-inbox_low[0];
    int ny = inbox_high[1]-inbox_low[1];
    int nz = inbox_high[2]-inbox_low[2];

    // Plan creation ...
    void *plan_d2z;
    int tuning_option = 0

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();
    if (tuning_option==0)
        plan_d2z = fftw_mpi_plan_dft_r2c_3d(nx, ny, nz, in, out, comm, FFTW_ESTIMATE);
    if (tuning_option==1)
        plan_d2z = fftw_mpi_plan_dft_r2c_3d(nx, ny, nz, in, out, comm, FFTW_MEASURE);
    timer[0] += MPI_Wtime();

    // FFT execution
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

	// for(int i = 0; i < niter; i++){
		fftw_execute(plan_d2z);
	// }

    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

    // Delete plan
	fftw_destroy_plan(plan_d2z);
}

//=====================  Complex-to-Real transform =========================


void compute_z2d_fftw( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, double *timer)
{
    int niter = 1;
    int nx = inbox_high[0]-inbox_low[0];
    int ny = inbox_high[1]-inbox_low[1];
    int nz = inbox_high[2]-inbox_low[2];

    // Plan creation ...
    void *plan_z2d;
    int tuning_option = 0

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();
    if (tuning_option==0)
        plan_z2d = fftw_mpi_plan_dft_c2r_3d(nx, ny, nz, in, out, comm, FFTW_ESTIMATE);
    if (tuning_option==1)
        plan_z2d = fftw_mpi_plan_dft_c2r_3d(nx, ny, nz, in, out, comm, FFTW_MEASURE);
    timer[0] += MPI_Wtime();

    // FFT execution
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

	// for(int i = 0; i < niter; i++){
		fftw_execute(plan_z2d);
	// }

    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

    // Delete plan
	fftw_destroy_plan(plan_z2d);

}


#else

void compute_z2z_fftw( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, double *timer)
{}

void compute_d2z_fftw( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, double *timer)
{}

void compute_z2d_fftw( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, double *timer)
{}

#endif


#endif  //! FIBER_BACKEND_FFTW_H