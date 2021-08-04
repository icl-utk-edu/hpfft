/*
* Harness software for benchmarking parallel FFT libraries.
* Benchmark: R2C 3-D transform, double precision, 2 MPI ranks.
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

    double *input  = malloc(size_inbox * sizeof(double));
    fiber_complex *output = calloc(size_outbox, sizeof(fiber_complex));

    for(i=0; i<size_inbox; i++)
        input[i] = (double) i;

    double timer[20];
    fiber_execute_d2z[heffte].function(box_low, box_high, box_low, box_high, comm, input, output, 0, timer);

    // Error computation after forward
    if (me == 1){
        int pass = 0;
        if (fabs(output[0].r + 512.0) > 1.E-11 || fabs(output[0].i) > 1.E-11)
            pass = 1;
        for(i=1; i<size_outbox; i++)
            if (fabs(output[i].r) > 1.E-11 || fabs(output[i].i) > 1.E-11)
                pass = 1;
        if (pass){
            printf("The computed transform deviates by more than the tolerance.\n");
            MPI_Abort(comm, 1);
        }
    }
    
    if(me == 0)
        printf("\t\t%s library computed a correct forward R2C 3-D transform \n ", backends[heffte] );

    /* 
    // Compute backwad transform
    for(i=0; i<size_inbox; i++) input[i] = 0.0;

    // Backward execution
    fiber_execute_z2d[heffte].function(box_low, box_high, box_low, box_high, comm, output, input, timer);

    // Error computation after backward
    double err = 0.0;
    for(i=0; i<size_inbox; i++)
        if (fabs(input[i] - (double) i) > err)
            err = fabs(input[i] - (double) i);

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
