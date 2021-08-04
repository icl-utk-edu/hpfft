/*
* Harness software for benchmarking parallel FFT libraries.
* Benchmark: C2C 3-D transform, double precision, 2 MPI ranks.
* Autor: Alan Ayala - ICL, UTK.
*/

// Use this program to verify the correct integration of a third-party library
// make -j; mpirun -n 2 ./test3D_CPU_C2C

#include "fiber_backends.h"
#include "fiber_utils.h"

int main(int argc, char** argv){

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

    fiber_complex *input  = calloc(size_outbox, sizeof(fiber_complex));
    fiber_complex *output = calloc(size_outbox, sizeof(fiber_complex));


    for(i=0; i<size_inbox; i++){
        input[i].r = (double) i;
        input[i].i = (double) 0*i;
    }        

    input[0].r = 1;
    input[0].i = 2;

    double timer[20];

    if (me == 0){
        printf("\t\t_________________________________________________________ \n");
        printf("\t\t\t  FFT Infraestructure Benchmark for Exascale Research   \n");
        printf("\t\t\t           Fixed-size 3-D FFT Benchmark                 \n");
        printf("\t\t_________________________________________________________ \n");

        printf("\t\tFFT size   : 4x4x4         \n");
        printf("\t\tPrecision  : DOUBLE        \n");
        printf("\t\tComputing  : C2C 1 FORWARD and 1 BACKWARD \n");
        // printf("\t\tNum. Iter. : %d \n", n_iterations);

        printf("\t\t_________________________________________________ \n");
        printf("\t\t%10s \t %12s \t %12s \n", "LIBRARY", "Time Plan (s)", "Time FFT (s)");
        printf("\t\t_________________________________________________ \n");
    }

    /*
    * Benchmark multiple libraries at the time, some libs still in progress
    for (int i = 0; i < 9; i++)
    {
        fiber_execute_z2z[i].function(box_low, box_high, box_low, box_high, comm, input, output, timer);
        if (me == 0){
            printf("\t\t%10s \t %6.3e \t %6.3e \n", backends[i], timer[0], timer[1]);
            printf("\t\t------------------------------------------------- \n");
        }        
    }
    */

    /*
    * Or benchmark a single library:
    */
    fiber_execute_z2z[heffte].function(box_low, box_high, box_low, box_high, comm, input, output, 0, timer);
    // fiber_execute_z2z[p3dfft].function(box_low, box_high, box_low, box_high, comm, input, output);

    /*
    * We computed the solution by hand, so the library has to give the same result, which is:
    for(i=0; i<size_outbox; i++) 
        printf(" %g + %g i", output[i].r, output[i].i);
    printf("\n");
    */

    // Check correctness after forward
    // for(i=0; i<size_inbox; i++) 
    //     printf(" %g + %gi  | ", output[i].r, output[i].i);
    // printf("\n");

    if (me == 1){
        int pass = 0;
        if (fabs(output[0].r + 510.0) > 1.E-11 || fabs(output[0].i - 4.0) > 1.E-11)
            pass = 1;
        for(i=1; i<16; i++)
            if (fabs(output[i].r-2) > 1.E-11 || fabs(output[i].i-4) > 1.E-11)
                pass = 1;
        for(i=16; i<32; i++)
            if (fabs(output[i].r) > 1.E-11 || fabs(output[i].i) > 1.E-11)
                pass = 1;
        if (pass){
            printf("The computed transform deviates by more than the tolerance.\n");
            MPI_Abort(comm, 1);
        }
    }

    if(me == 0)
        printf("\t\t%s library computed a correct forward C2C 3-D transform \n ", backends[heffte] );


/* 
    // Compute backwad transform
    for(i=0; i<size_inbox; i++){
        input[i].r = 0.0;
        input[i].i = 0.0;
    }        
    fiber_execute_z2z[heffte].function(box_low, box_high, box_low, box_high, comm, output, input, 1, timer);

    // Check correctness after backward

    printf("\n");
    printf("\n");

    for(i=0; i<size_inbox; i++) 
        printf(" %g + %gi  | ", input[i].r, input[i].i);
    printf("\n");

    double err = 0.0;
    for(i=1; i<size_inbox; i++)
        if (fabs(input[i].r - (double) i - 1) > err)
            err = fabs(input[i].r - (double) i - 1);

    err = fmax( fabs(input[0].r - (double) 1) , err);
    err = fmax( fabs(input[0].r - (double) 2) , err);


    // Print errors
    if (me == 0) printf("rank 0 computed error: %1.6le\n", err);
    MPI_Barrier(comm);
    if (me == 1) printf("rank 1 computed error: %1.6le\n", err);
*/

    // Data deallocation 
    free(input);
    free(output);

    MPI_Finalize();

    return 0;
}
