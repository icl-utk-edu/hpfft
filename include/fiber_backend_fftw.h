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
                  void const *in, void *out, int fftw_switch, double *timer)
{
    // FFTW tuning flags: 
    // FFTW_MEASURE is more time consuming than FFTW_ESTIMATE, it runs and measures the execution time of several FFTs in order to find the
    // best way to compute the transform, given the FFT size and resources.
    
    int niter = 1;

    ptrdiff_t local_nx, local_x_start, nx, ny, nz;

    // nx = inbox_high[0] - inbox_low[0] + 1;
    // ny = inbox_high[1] - inbox_low[1] + 1;
    // nz = inbox_high[2] - inbox_low[2] + 1;

    // Global size, to come as input later on
    nx=ny=nz=4;
    printf(" global size: %d \t %d \t %d  \n" , nx, ny, nz);

    local_nx = inbox_high[0] - inbox_low[0] + 1;
    local_x_start = 0;

    printf("local_nx = %d \n", local_nx);
    printf("local_nx_start = %d \n", local_x_start);


    // Plan creation ...
    void *plan_z2z;
    fftw_mpi_init();

    printf(" flags: %d \t %d \t %d  \n" , FFTW_FORWARD, FFTW_ESTIMATE, FFTW_MEASURE);
    printf("Tuning flag: %d \n", fftw_switch);

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();
    if (fftw_switch==0)
        plan_z2z = fftw_mpi_plan_dft_3d(nx, ny, nz, in, out, comm, FFTW_FORWARD, FFTW_ESTIMATE);
    if (fftw_switch==1)
        // plan_z2z = fftw_mpi_plan_dft_3d(nx, ny, nz, in, out, comm, FFTW_FORWARD, FFTW_MEASURE);
        plan_z2z = fftw_mpi_plan_dft_3d(nx, ny, nz, in, out, comm, FFTW_BACKWARD, FFTW_ESTIMATE);
    timer[0] += MPI_Wtime();

    // FFT execution
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

	// for(int i = 0; i < niter; i++){
		fftw_execute(plan_z2z);
		// fftw_mpi_execute_dft(plan_z2z, (fftw_complex *) in, (fftw_complex *) out);
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
                  double const *in, void *out, int fftw_switch, double *timer)
{

    int niter = 1;
    ptrdiff_t local_nx, local_x_start, nx, ny, nz;
    nx=ny=nz=4;
    local_nx = inbox_high[0] - inbox_low[0] + 1;
    local_x_start = 0;

    // nx = inbox_high[0] - inbox_low[0] + 1;
    // ny = inbox_high[1] - inbox_low[1] + 1;
    // nz = inbox_high[2] - inbox_low[2] + 1;

    // Plan creation ...
    void *plan_d2z;
    
    fftw_mpi_init();

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();

    printf(" flags: %d \t %d \t %d  \n" , FFTW_FORWARD, FFTW_ESTIMATE, FFTW_MEASURE);
    printf("Tuning flag: %d \n", fftw_switch);

    if (fftw_switch==0)
        plan_d2z = fftw_mpi_plan_dft_r2c_3d(nx, ny, nz, in, out, comm, FFTW_ESTIMATE);
    if (fftw_switch==1)
        plan_d2z = fftw_mpi_plan_dft_r2c_3d(nx, ny, nz, in, out, comm, FFTW_MEASURE);
    timer[0] += MPI_Wtime();


    // FFT execution
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

	// for(int i = 0; i < niter; i++){
        printf("Executing : \n");
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
                  void const *in, double *out, int fftw_switch, double *timer)
{

    int niter = 1;
    ptrdiff_t local_nx, local_x_start, nx, ny, nz;
    nx=ny=nz=4;
    local_nx = inbox_high[0] - inbox_low[0] + 1;
    local_x_start = 0;

    // nx = inbox_high[0] - inbox_low[0] + 1;
    // ny = inbox_high[1] - inbox_low[1] + 1;
    // nz = inbox_high[2] - inbox_low[2] + 1;

    // Plan creation ...
    void *plan_z2d;
    fftw_mpi_init();

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();
    if (fftw_switch==0)
        plan_z2d = fftw_mpi_plan_dft_c2r_3d(nx, ny, nz, in, out, comm, FFTW_ESTIMATE);
    if (fftw_switch==1)
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