
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

    hpfft_complex *input  = calloc(size_inbox, sizeof(hpfft_complex));

    // Data Initialization
    for(i=0; i<size_inbox; i++){
        input[i].r = (double) i;
        input[i].i = (double) 0*i;
    }        

    input[0].r = 1;
    input[0].i = 2;

    // for(i=0; i<size_inbox; i++) 
    //     printf("  %g+%gi  \t ", input[i].r, input[i].i);
    // printf("\n");        
    // printf("\n");      

    double timer[20];

    // ********************************
    // Compute forward (D2Z) transform
    // ********************************

    double c2c_time[18];

    for(i=0; i<9; i++) {   // Benchmark all 9 backends 
        if(i>=4 && i<=7 ){ // currently missing these backends
            c2c_time[2*i] = c2c_time[2*i+1] = -1;
        }
        else{
            hpfft_execute_z2z[i].function(box_low, box_high, box_low, box_high, comm, input, input, 0, timer);
            c2c_time[2*i] = timer[0]; 
            c2c_time[2*i+1] = timer[1];
        }
    }


    if (me == 0){
        printf("\t\t______________________________________________________ \n");
        printf("\t\t FFT Infraestructure Benchmark for Exascale Research   \n");
        printf("\t\t           Fixed-size 3-D FFT Benchmark                \n");
        printf("\t\t______________________________________________________ \n");

        printf("\t\tFFT size   : 4x4x4         \n");
        printf("\t\tPrecision  : DOUBLE        \n");
        printf("\t\tComputing  : C2C 1 FORWARD \n");

        // ********************************
        // Compute forward (Z2Z) transform
        // ********************************

        printf("\t\t_________________________________________________ \n");
        printf("\t\t%10s \t %12s \t %12s \n", "LIBRARY", "Time Plan (s)", "Time FFT (s)");
        printf("\t\t_________________________________________________ \n");

    }

    for(i=0; i<9; i++)
        if(me==0)
            printf(" \t\t %10s \t %g \t  %g \n ", backends[i], c2c_time[2*i], c2c_time[2*i+1]);

    // // Output after forward
    // for(i=0; i<size_outbox; i++) 
    //     printf(" %g+%gi \t ", input[i].r, input[i].i);
    // printf("\n");


    free(input);

    MPI_Finalize();

    return 0;
}



