/* ----------------------------------------------------------------------
    C interface to AccFFT library, 3-D FFT functions for CPU hardware
    Author: Alan Ayala, aayala@icl.utk.edu

    Add this file into accfft/include folder.
    Allows calling AccFFT from a C program.
------------------------------------------------------------------------- */

#include <mpi.h>
#include <stdio.h>

#include "accfft.h"
#include "accfft_common.h"
#include "accfft_c.h"

/* ----------------------------------------------------------------------
    Initialize AccFFT C-Interface 
------------------------------------------------------------------------- */

int accfft_init_c(){
    printf(" Calling AccFFT C-Interface \n");
    return 0;
}


/* ----------------------------------------------------------------------
    Create AccFFT plan
------------------------------------------------------------------------- */

void accfft_create_plan(int * n, Complex * data, Complex * data_out,
		MPI_Comm c_comm, unsigned flags, void ** plan)
{
    accfft_plan * ptr = accfft_plan_dft_3d_c2c(n, data, data, c_comm, ACCFFT_MEASURE);
    *plan = (void *) ptr;
}

/* ----------------------------------------------------------------------
    Execute AccFFT plan
------------------------------------------------------------------------- */

void accfft_compute(void * ptr, Complex * data, Complex * data_out, int flag)
{
    accfft_plan * plan = (accfft_plan *) ptr;
    accfft_execute_c2c(plan, flag, data, data_out);
}


/* ----------------------------------------------------------------------
    Delete plan
------------------------------------------------------------------------- */

void accfft_destroy_plan_c(void * ptr)
{
    accfft_plan *plan = (accfft_plan *) ptr;
	accfft_destroy_plan(plan);
}