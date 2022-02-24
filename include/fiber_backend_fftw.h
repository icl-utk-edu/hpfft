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


/*!
 * \ingroup CPU_libraries
 * \addtogroup fiber_cpu Backend fftw3
 *
 * Wrappers and template specializations related to the FFTW backend.
 * Requires CMake option:
 * \code
 *  -D Heffte_ENABLE_FFTW=ON
 * \endcode
 * Flag:
 * The fftw_options[0] flag works as follows:
 * For C2C transforms: it chooses between FFTW_FORWARD and FFTW_BACKWARD flags
 * For R2C transforms: it chooses between FFTW_ESTIMATE and FFTW_MEASURE flags
 * nx = fftw_options[1]
 * ny = fftw_options[2]
 * nz = fftw_options[3]
 */

//=================== Initialization (if required) ============================
int init_fftw(int physical, int nx, int ny, int nz, int p_row, int p_col){
    fftw_mpi_init();
    return(0);
}

int finalize_fftw(){
    return(0);
}

//=====================  Complex-to-Complex transform =========================

void compute_z2z_fftw( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, int *fftw_options, double *timer)
{
    // FFTW tuning flags: 
    // FFTW_MEASURE is more time consuming than FFTW_ESTIMATE, it runs and measures the execution time of several FFTs in order to find the
    // best way to compute the transform, given the FFT size and resources.
    
    int niter = 1;
    ptrdiff_t nx, ny, nz;
    nx = fftw_options[1];
    ny = fftw_options[2];
    nz = fftw_options[3];
    nx=ny=nz=4;

    // Global size, to come as input later on
    // Plan creation ...
    void *plan_z2z;

    // printf(" flags: %d \t %d \t %d  \n" , FFTW_FORWARD, FFTW_ESTIMATE, FFTW_MEASURE);
    // printf("FFTW first tuning flag: %d \n", fftw_options[0]);

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();
    if (fftw_options[0]==0){
        printf("Forward 3-D C2C transform using FFTW\n");
        plan_z2z = fftw_mpi_plan_dft_3d(nx, ny, nz, in, out, comm, FFTW_FORWARD, FFTW_ESTIMATE);
    }        
    if (fftw_options[0]==1){
        printf("Backward 3-D C2C transform using FFTW\n");
        // plan_z2z = fftw_mpi_plan_dft_3d(nx, ny, nz, in, out, comm, FFTW_FORWARD, FFTW_MEASURE);
        plan_z2z = fftw_mpi_plan_dft_3d(nx, ny, nz, in, out, comm, FFTW_BACKWARD, FFTW_ESTIMATE);
    }        
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
                  double const *in, void *out, int *fftw_options, double *timer)
{

    int niter = 1;
    ptrdiff_t local_nx, local_x_start, nx, ny, nz;
    nx=ny=nz=4;
    local_nx = inbox_high[0] - inbox_low[0] + 1;
    local_x_start = 0;

    // Plan creation ...
    void *plan_d2z;
    
    fftw_mpi_init();

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();

    printf(" flags: %d \t %d \t %d  \n" , FFTW_FORWARD, FFTW_ESTIMATE, FFTW_MEASURE);
    printf("First tuning flag: %d \n", fftw_options[0]);

    if (fftw_options[0]==0)
        plan_d2z = fftw_mpi_plan_dft_r2c_3d(nx, ny, nz, in, out, comm, FFTW_ESTIMATE);
    if (fftw_options[0]==1)
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
                  void const *in, double *out, int *fftw_options, double *timer)
{

    int niter = 1;
    ptrdiff_t local_nx, local_x_start, nx, ny, nz;
    nx=ny=nz=4;
    local_nx = inbox_high[0] - inbox_low[0] + 1;
    local_x_start = 0;

    // Plan creation ...
    void *plan_z2d;
    fftw_mpi_init();

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();
    if (fftw_options[0]==0)
        plan_z2d = fftw_mpi_plan_dft_c2r_3d(nx, ny, nz, in, out, comm, FFTW_ESTIMATE);
    if (fftw_options[0]==1)
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

int init_fftw(int physical, int nx, int ny, int nz, int p_row, int p_col)
{}

int finalize_fftw()
{}


void compute_z2z_fftw( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, int *fftw_options, double *timer)
{}

void compute_d2z_fftw( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, int *fftw_options, double *timer)
{}

void compute_z2d_fftw( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, int *fftw_options, double *timer)
{}

#endif


#endif  //! FIBER_BACKEND_FFTW_H
