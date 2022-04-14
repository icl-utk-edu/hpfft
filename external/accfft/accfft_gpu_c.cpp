/* ----------------------------------------------------------------------
    C interface to AccFFT library, 3-D FFT functions for GPU hardware
    Author: Alan Ayala, aayala@icl.utk.edu

    Add this file into accfft/include folder.
    Allows calling AccFFT from a C program.
------------------------------------------------------------------------- */

#include <mpi.h>
#include <stdio.h>

#include "accfft_gpu.h"
#include "accfft_c.h"

#include <cuda_runtime_api.h>


/* ----------------------------------------------------------------------
    Create AccFFT plan GPU-case
------------------------------------------------------------------------- */

void accfft_create_plan_gpu(int * n, double * data, double * data_out,
		MPI_Comm c_comm, unsigned flags, void ** plan)
{
    accfft_plan_gpu * ptr = accfft_plan_dft_3d_c2c_gpu(n, (Complex *) data, (Complex *) data_out, c_comm, ACCFFT_MEASURE);
    *plan = (void *) ptr;
}

/* ----------------------------------------------------------------------
    Execute AccFFT plan GPU-case
------------------------------------------------------------------------- */

void accfft_compute_gpu(void * ptr, double * data, double * data_out, int flag)
{
    accfft_plan_gpu * plan = (accfft_plan_gpu *) ptr;
    accfft_execute_c2c_gpu(plan, flag, (Complex *) data, (Complex *) data_out);
}


/* ----------------------------------------------------------------------
    Delete plan GPU-case
------------------------------------------------------------------------- */

void accfft_destroy_plan_gpu_c(void * ptr)
{
    accfft_plan_gpu *plan = (accfft_plan_gpu *) ptr;
	accfft_destroy_plan_gpu(plan);
}
