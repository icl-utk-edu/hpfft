/*
* ---------------
* AccFFT backend
* ---------------
*/
#ifndef FIBER_BACKEND_ACCFFT_H
#define FIBER_BACKEND_ACCFFT_H

// #include <accfft.h>


//=====================  Complex-to-Complex transform =========================

void compute_z2z_accfft( int const inbox_low[3], int const inbox_high[3],
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

    printf("Benchmarking AccFFT: Get Complex-to-Complex creating a wrapper within AccFFT \n");

}    


// void compute_d2z_accfft( int const * , int const *, int const * , int const * , MPI_Comm const , double const *, void *);
void compute_d2z_accfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, double *timer)
{

    printf("Pending D2Z AccFFT \n");

	// int c_dims[2] = { 0 };
	// MPI_Dims_create(nprocs, 2, c_dims);

    // print(" C =  %d -- %d", c_dims[0], c_dims[1]);

	// accfft_create_comm(comm, c_dims, &c_comm);

    // // Plan definition
	// accfft_plan * plan = accfft_plan_dft_3d_r2c(n, in, (double*) out, c_comm, ACCFFT_MEASURE);

    // int status = accfft_plan_create(accfft_BACKEND_FFTW, inbox_low, inbox_high, NULL, outbox_low, outbox_high, NULL, comm, NULL, &plan);


	// accfft_execute_r2c(plan, data, data_hat);

    // status = accfft_plan_destroy(plan);
    // accfft_destroy_plan(plan);
	// accfft_cleanup();

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

}


void compute_z2d_accfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, double *timer)
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

    printf("Pending AccFFT Z2D \n");

}



#endif  //! FIBER_BACKEND_ACCFFT_H