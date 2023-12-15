
/*
* Harness software for benchmarking parallel FFT libraries.
* Benchmark: R2C 3-D transform, double precision, 2 MPI ranks.
* Autor: Alan Ayala - ICL, UTK.
*/

// Use this program to benchmark R2C FFT performance for all libraries

#include "hpfft_backends.h"
#include "hpfft_utils.h"

int main(int argc, char** argv){

    // Choose backend type
    int my_backend  = heffte;

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int me;
    MPI_Comm_rank(comm, &me);

    int num_ranks;
    MPI_Comm_size(comm, &num_ranks);

    if (num_ranks != 2){
        if (me == 0) printf(" example is set to 2 mpi ranks, exiting \n");
        return 0;
    }

    int i;
    // Initial configuration setup
    int box_low[3]  = {0, 0, 0};
    int box_high[3] = {3, 3, 3};

    if (me == 0) // split across the last dimension
        box_high[2] = 1;
    else
        box_low[2] = 2;

    int size_inbox  = 32;
    int size_outbox = 32;

    double *input  = malloc(size_inbox * sizeof(double));
    hpfft_complex *output = calloc(size_outbox, sizeof(hpfft_complex));

    // Data Initialization
    for(i=0; i<size_inbox; i++)
        input[i] = (double) i;

    double timer[20];

    // ********************************
    // Compute forward (D2Z) transform
    // ********************************

    if (me == 0){
        printf("\t\t______________________________________________________ \n");
        printf("\t\t FFT Infraestructure Benchmark for Exascale Research   \n");
        printf("\t\t           Fixed-size 3-D FFT Benchmark                \n");
        printf("\t\t______________________________________________________ \n");

        printf("\t\tFFT size   : 4x4x4         \n");
        printf("\t\tPrecision  : DOUBLE        \n");
        printf("\t\tComputing  : R2C 1 FORWARD and 1 BACKWARD \n");
        // printf("\t\tNum. Iter. : %d \n", n_iterations);

        printf("\t\t_________________________________________________ \n");
        printf("\t\t%10s \t %12s \t %12s \n", "LIBRARY", "Time Plan (s)", "Time FFT (s)");
        printf("\t\t_________________________________________________ \n");
    }


//  Benchmark multiple libraries at the time, some libs still in progress
    // for (int i = 0; i < 9; i++)
    // {
    //     hpfft_execute_z2z[i].function(box_low, box_high, box_low, box_high, comm, input, output, timer);
    //     if (me == 0){
    //         printf("\t\t%10s \t %6.3e \t %6.3e \n", backends[i], timer[0], timer[1]);
    //         printf("\t\t------------------------------------------------- \n");
    //     }        
    // }

    // Data deallocation 
    free(input);
    free(output);

    MPI_Finalize();

    return 0;
}



