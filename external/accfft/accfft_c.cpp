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

int accfft_init_c()
{
    return 0;
}


/* ----------------------------------------------------------------------
    Create AccFFT Communicator
------------------------------------------------------------------------- */

int accfft_create_comm_c(MPI_Comm in_comm, int * c_dims, MPI_Comm *c_comm) 
{
	int nprocs, procid;
	MPI_Comm_rank(in_comm, &procid);
	MPI_Comm_size(in_comm, &nprocs);

	if (c_dims[0] * c_dims[1] != nprocs) {
		c_dims[0] = 0;
		c_dims[1] = 0;
		MPI_Dims_create(nprocs, 2, c_dims);
	}

	// Create Cartesian Communicator
	int period[2], reorder;
	int coord[2];
	period[0] = 0;
	period[1] = 0;
	reorder = 1;

	MPI_Cart_create(in_comm, 2, c_dims, period, reorder, c_comm);

    return 0;
}


/* ----------------------------------------------------------------------
    Create AccFFT plan
------------------------------------------------------------------------- */

void accfft_create_plan(int * n, double * data, double * data_out,
		MPI_Comm c_comm, unsigned flags, void ** plan)
{
    accfft_plan * ptr = accfft_plan_dft_3d_c2c(n, (Complex *) data, (Complex *) data, c_comm, ACCFFT_MEASURE);
    *plan = (void *) ptr;
}

/* ----------------------------------------------------------------------
    Execute AccFFT plan
------------------------------------------------------------------------- */

void accfft_compute(void * ptr, double * data, double * data_out, int flag)
{
    
    accfft_plan * plan = (accfft_plan *) ptr;
    accfft_execute_c2c(plan, flag, (Complex *) data, (Complex *) data_out);
}


/* ----------------------------------------------------------------------
    Delete plan
------------------------------------------------------------------------- */

void accfft_destroy_plan_c(void * ptr)
{
    accfft_plan *plan = (accfft_plan *) ptr;
	accfft_destroy_plan(plan);
}
