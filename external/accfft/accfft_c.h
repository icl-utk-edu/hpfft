/* ----------------------------------------------------------------------
    C interface to AccFFT library, 3-D FFT functions
    Author: Alan Ayala, aayala@icl.utk.edu

    Add this file into accfft/include folder.
    Allows calling AccFFT from a C program.
------------------------------------------------------------------------- */

#ifndef ACCFFT_C_H
#define ACCFFT_C_H

#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <fftw3.h>

typedef double Complex[2];

/* ifdefs allow this file to be included in a C program */

#ifdef __cplusplus
extern "C" {
#endif


/* ----------------------------------------------------------------------
    Initialize AccFFT C-Interface, return 0 if successful.
------------------------------------------------------------------------- */

int accfft_init_c();

/* ----------------------------------------------------------------------
    Create AccFFT plan CPU-case
------------------------------------------------------------------------- */

void accfft_create_plan(int * n, Complex * data, Complex * data_out,
		MPI_Comm c_comm, unsigned flags, void ** plan);


/* ----------------------------------------------------------------------
    Create AccFFT plan GPU-case
------------------------------------------------------------------------- */

void accfft_plan_dft_3d_c2c_gpu(int * n, Complex * data, Complex * data_out,
		MPI_Comm c_comm, unsigned flags, void ** plan);



/* ----------------------------------------------------------------------
    Execute AccFFT plan CPU-case
------------------------------------------------------------------------- */

void accfft_compute(void * ptr, Complex * data, Complex * data_out, int flag);


/* ----------------------------------------------------------------------
    Execute AccFFT plan GPU-case
------------------------------------------------------------------------- */

void accfft_compute_gpu(void * ptr, Complex * data, Complex * data_out, int flag);



/* ----------------------------------------------------------------------
    Delete plan CPU-case
------------------------------------------------------------------------- */

void accfft_destroy_plan_c(void * ptr);


/* ----------------------------------------------------------------------
    Delete plan GPU-case
------------------------------------------------------------------------- */

void accfft_destroy_plan_gpu_c(void * ptr);



#ifdef __cplusplus
}
#endif


#endif
