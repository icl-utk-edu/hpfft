/*
* Harness software for benchmarking parallel FFT libraries.
* Benchmark: R2C 3-D transform, double precision, 2 MPI ranks.
* Autor: Alan Ayala - ICL, UTK.
------------------------------------------------------------------------------
Use this program to verify the correct integration of a third-party library
make clean; make -j; mpirun -n 2 ./test3D_CPU_C2C <library>
*/

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

    if (me == 0){
        printf("\t\t______________________________________________________ \n");
        printf("\t\t FFT Infraestructure Benchmark for Exascale Research   \n");
        printf("\t\t           Fixed-size 3-D FFT Benchmark                \n");
        printf("\t\t______________________________________________________ \n");
        printf("\t\tLibrary    : %s         \n", argv[1]);
        printf("\t\tFFT size   : 4x4x4         \n");
        printf("\t\tPrecision  : DOUBLE        \n");
        printf("\t\tComputing  : R2C 1 FORWARD and 1 BACKWARD \n");
    }

    // Get backend type from user
    int my_backend  = fiber_get_backend(argv[1]);

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

    // Data Initialization
    for(i=0; i<size_inbox; i++)
        input[i] = (double) i;

    for(i=0; i<size_inbox; i++) 
        printf("  %g \t ", input[i]);
    printf("\n");        
    printf("\n");        

    double timer[20];

    // ********************************
    // Compute forward (D2Z) transform
    // ********************************
    fiber_execute_d2z[my_backend].function(box_low, box_high, box_low, box_high, comm, input, output, 0, timer);

    // Output after forward
    for(i=0; i<size_outbox; i++) 
        printf(" %g+%gi  \t ", output[i].r, output[i].i);
    printf("\n");
    printf("\n");        

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
    
    if(me == 1)
        printf("\t\t%s library computed a correct forward R2C 3-D transform \n \n ", backends[my_backend] );


    // ********************************
    // Compute backwad (Z2D) transform
    // ********************************
    for(i=0; i<size_inbox; i++) input[i] = 0.0;

    // Backward execution
    fiber_execute_z2d[my_backend].function(box_low, box_high, box_low, box_high, comm, output, input, timer);

    // Output after backward
    for(i=0; i<size_inbox; i++) 
        printf("  %g \t ", input[i]);
    printf("\n");
    printf("\n");        

    // Error computation after backward
    double err = 0.0;
    for(i=0; i<size_inbox; i++)
        if (fabs(input[i] - (double) i) > err)
            err = fabs(input[i] - (double) i);

    printf("%s: rank %d computed error |X - ifft(fft(X)) |: %1.6le\n", backends[my_backend], me, err);

    // Data deallocation 
    free(input);
    free(output);

    MPI_Finalize();

    return 0;
}
