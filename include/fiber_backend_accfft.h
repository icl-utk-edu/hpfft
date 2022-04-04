/*
* ---------------
* AccFFT backend
* ---------------
*/
#ifndef FIBER_BACKEND_ACCFFT_H
#define FIBER_BACKEND_ACCFFT_H


#if defined(FIBER_ENABLE_ACCFFT)
#include <accfft_c.h>


//=================== Initialization (if required) ============================
int init_accfft(int option){
    return(0);
}


//=====================  Complex-to-Complex transform =========================

void compute_z2z_accfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, int *accfft_options, double *timer)
{

    int n[3] = {accfft_options[1], accfft_options[2], accfft_options[3]};

    // Plan creation ...
    void *plan;

    // -------- Create plan for FORWARD transform ----------
    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();

    accfft_create_plan(n, in, out, comm, 2, &plan); // Using flag: ACCFFT_MEASURE=2

    MPI_Barrier(comm);
    timer[0] += MPI_Wtime();
    // --------------------------------------------------------


    // ------------- Execute FORWARD transform ----------------
    heffte_forward_z2z(plan, in, out, Heffte_SCALE_NONE);

    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

    for (int i=0; i<accfft_options[8]; ++i)
        accfft_compute(plan, in, out, -1);

    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();
    // --------------------------------------------------------


    // Delete plan ...
    accfft_destroy_plan_c(plan);
}    


// void compute_d2z_accfft( int const * , int const *, int const * , int const * , MPI_Comm const , double const *, void *);
void compute_d2z_accfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, int *accfft_options, double *timer)
{

    /*

    printf("R2C AccFFT CPU support \n");

	int c_dims[2] = { 0 };
	MPI_Dims_create(nprocs, 2, c_dims);
    print(" C =  %d -- %d", c_dims[0], c_dims[1]);
	accfft_create_comm(comm, c_dims, &c_comm);

    // Plan creation ...
    // void *plan;

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();

    // plan create
	accfft_plan * plan = accfft_plan_dft_3d_r2c(n, in, (double*) out, c_comm, ACCFFT_MEASURE);
    int status = accfft_plan_create(accfft_BACKEND_FFTW, inbox_low, inbox_high, NULL, outbox_low, outbox_high, NULL, comm, NULL, &plan);

    MPI_Barrier(comm);
    timer[0] += MPI_Wtime();



    // FFT execution ...
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();
    // compute FFT 
	accfft_execute_r2c(plan, data, data_hat);

    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

    // Delete plan ...
    status = accfft_plan_destroy(plan);
	accfft_cleanup();
    */
}


void compute_z2d_accfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, int *accfft_options, double *timer)
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

#else
int init_accfft(int option)
{}

void compute_z2z_accfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, int *accfft_options, double *timer)
{}

void compute_d2z_accfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, int *accfft_options, double *timer)
{}

void compute_z2d_accfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, int *accfft_options, double *timer)
{}

#endif