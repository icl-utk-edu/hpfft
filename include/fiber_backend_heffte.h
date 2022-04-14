/*
* ---------------
* heFFTe backend
* ---------------
*/
#ifndef FIBER_BACKEND_HEFFTE_H
#define FIBER_BACKEND_HEFFTE_H

#include <stdio.h>

//=====================  Complex-to-Complex transform =========================

#if defined(FIBER_ENABLE_HEFFTE)
#include "heffte.h"

//=================== Initialization (if required) ============================
int init_heffte(int option){
    return(0);
}

void compute_z2z_heffte( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, int *heffte_options, double *timer)
{
    // Plan definition
    heffte_plan plan;

    MPI_Barrier(comm);
    double t;
    int status;

    t = -MPI_Wtime();
    
    switch (heffte_options[4])
    {
    //   case 0:
        // Forward 3-D C2C transform using heFFTe with STOCK backend
    //     status = heffte_plan_create(Heffte_BACKEND_STOCK, inbox_low, inbox_high, NULL, outbox_low, outbox_high, NULL, comm, NULL, &plan);
    //     break;

      case 1:
        // Forward 3-D C2C transform using heFFTe with FFTW backend
        status = heffte_plan_create(Heffte_BACKEND_FFTW, inbox_low, inbox_high, NULL, outbox_low, outbox_high, NULL, comm, NULL, &plan);
        break;
        
    //   case 2:
        // Forward 3-D C2C transform using heFFTe with MKL backend
    //     status = heffte_plan_create(Heffte_BACKEND_MKL, inbox_low, inbox_high, NULL, outbox_low, outbox_high, NULL, comm, NULL, &plan);
    //     break;

      case 3:
        // Forward 3-D C2C transform using heFFTe with CUFFT backend
        status = heffte_plan_create(Heffte_BACKEND_CUFFT, inbox_low, inbox_high, NULL, outbox_low, outbox_high, NULL, comm, NULL, &plan);
        break;        

    //   case 4:
    //     printf("Forward 3-D C2C transform using heFFTe with ROCFFT backend \n");
    //     status = heffte_plan_create(Heffte_BACKEND_ROCFFT, inbox_low, inbox_high, NULL, outbox_low, outbox_high, NULL, comm, NULL, &plan);
    //     break;    

      default:
        printf("ERROR: Invalid heFFTe backend!\n");
        break;
    }
    t += MPI_Wtime();

    MPI_Barrier(comm);
    MPI_Reduce(&t, &timer[0], 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    int workspace_size = heffte_size_workspace(plan);
    // printf("Workspace size = %d \n", heffte_size_workspace(plan));
    double *workspace = malloc(2 * workspace_size * sizeof(double));

    if (status != Heffte_SUCCESS){
        printf("Failed at heffte_plan_create() with error code: %d\n", status);
        MPI_Abort(comm, 1);
    }

    MPI_Barrier(comm);
    // heffte_forward_z2z(plan, in, out, Heffte_SCALE_NONE);
    heffte_forward_z2z_buffered(plan, in, out, workspace, Heffte_SCALE_NONE);

    // FFT execution
    t = -MPI_Wtime();
    
    if(heffte_options[0] == 0){
        // printf("heFFTe Forward Transform \n");
        for (int i=0; i<heffte_options[8]; ++i){
            // heffte_forward_z2z(plan, in, in, Heffte_SCALE_NONE);
            heffte_forward_z2z_buffered(plan, in, out, workspace, Heffte_SCALE_NONE);
            // printf("Iteration %d, timing %g \n", i, MPI_Wtime()+t);
        }
        t += MPI_Wtime();

        MPI_Barrier(comm);
        MPI_Reduce(&t, &timer[1], 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        // out = in;
    }

    if(heffte_options[0] == 1){
        // printf("heFFTe Backward Transform \n");
        // heffte_backward_z2z(plan, in, out, Heffte_SCALE_NONE);
        heffte_backward_z2z_buffered(plan, in, out, workspace, Heffte_SCALE_NONE);
    }    

    status = heffte_plan_destroy(plan);
    if (status != Heffte_SUCCESS){
        printf("Failed at heffte_plan_destroy() with error code: %d\n", status);
        MPI_Abort(comm, 1);
    }
}

//=====================  Real-to-Complex transform =========================


// void compute_d2z_heffte( int const * , int const *, int const * , int const * , MPI_Comm const , double const *, void *);
void compute_d2z_heffte( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, int *heffte_options, double *timer)
{

    // Plan definition
    heffte_plan plan;

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();
    int status;
    if(heffte_options[0] == 1){
        printf("Forward 3-D R2C transform using heFFTe with CUFFT backend \n");
        // status = heffte_plan_create_r2c(Heffte_BACKEND_CUFFT, inbox_low, inbox_high, NULL, outbox_low, outbox_high, NULL, 0, comm, NULL, &plan);
        status = heffte_plan_create(Heffte_BACKEND_CUFFT, inbox_low, inbox_high, NULL, outbox_low, outbox_high, NULL, comm, NULL, &plan);
    }        
    else{        
        printf("Forward 3-D R2C transform using heFFTe with FFTW backend \n");
        // status = heffte_plan_create_r2c(Heffte_BACKEND_FFTW, inbox_low, inbox_high, NULL, outbox_low, outbox_high, NULL, 0, comm, NULL, &plan);
        status = heffte_plan_create(Heffte_BACKEND_FFTW, inbox_low, inbox_high, NULL, outbox_low, outbox_high, NULL, comm, NULL, &plan);
    }        
    MPI_Barrier(comm);
    timer[0] = +MPI_Wtime();

    if (status != Heffte_SUCCESS){
        printf("Failed at heffte_plan_create() with error code: %d\n", status);
        MPI_Abort(comm, 1);
    }

    printf("size in = %d \n" , heffte_size_inbox(plan) );
    printf("size ou = %d \n" , heffte_size_outbox(plan) );
    printf("size wr = %d \n" , heffte_size_workspace(plan) );
    printf("backend = %d \n" , heffte_get_backend(plan)  );
    printf("r2c = %d \n" , heffte_is_r2c(plan)  );

    // FFT execution
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();
    // compute
    heffte_forward_d2z(plan, in, out, Heffte_SCALE_NONE);
    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

    status = heffte_plan_destroy(plan);
    if (status != Heffte_SUCCESS){
        printf("Failed at heffte_plan_destroy() with error code: %d\n", status);
        MPI_Abort(comm, 1);
    }
}

void compute_z2d_heffte( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, int *heffte_options, double *timer)
{
    // Plan definition
    heffte_plan plan;
    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();

    int status = heffte_plan_create(Heffte_BACKEND_FFTW, inbox_low, inbox_high, NULL, outbox_low, outbox_high, NULL, comm, NULL, &plan);
    MPI_Barrier(comm);
    timer[0] = +MPI_Wtime();

    if (status != Heffte_SUCCESS){
        printf("Failed at heffte_plan_create() with error code: %d\n", status);
        MPI_Abort(comm, 1);
    }

    // FFT execution
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();
    // compute
    heffte_backward_z2d(plan, in, out, Heffte_SCALE_FULL);
    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

    status = heffte_plan_destroy(plan);
    if (status != Heffte_SUCCESS){
        printf("Failed at heffte_plan_destroy() with error code: %d\n", status);
        MPI_Abort(comm, 1);
    }
}


// ========================================================================================

#else
int init_heffte(int option)
{}

void compute_z2z_heffte( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, void *out, int *heffte_options, double *timer)
{}

void compute_d2z_heffte( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  double const *in, void *out, int *heffte_options, double *timer)
{}

void compute_z2d_heffte( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
                  void const *in, double *out, int *heffte_options, double *timer)
{}
#endif

#endif  //! FIBER_BACKEND_HEFFTE_H