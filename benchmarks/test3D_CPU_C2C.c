/*
* Harness software for benchmarking parallel FFT libraries.
* Benchmark: C2C 3-D transform, double precision, 2 MPI ranks.
* Autor: Alan Ayala - ICL, UTK.
------------------------------------------------------------------------------
Use this program to verify the correct integration of a third-party library.
All libraries must pass this error validation tester.
make clean; make -j; mpirun -n 2 ./test3D_CPU_C2C <library>
*/

#include "fiber_backends.h"
#include "fiber_utils.h"

int main(int argc, char** argv){
   
    // Get backend type from user
    int my_backend  = fiber_get_backend(argv[1]);

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
        printf("\t\tComputing  : C2C 1 FORWARD and 1 BACKWARD \n");
    }

    int i;
    // Initial configuration setup
    int box_low[3]  = {0, 0, 0};
    int box_high[3] = {3, 3, 3};
    int fftsize = 64;

    if (me == 0) // split across the last dimension
        box_high[2] = 1;
    else
        box_low[2] = 2;

    int size_inbox  = 32;
    int size_outbox = 32;

    fiber_complex *input  = calloc(size_outbox, sizeof(fiber_complex));
    fiber_complex *output = calloc(size_outbox, sizeof(fiber_complex));

    // Data Initialization
    for(i=0; i<size_inbox; i++){
        input[i].r = (double) i;
        input[i].i = (double) 0*i;
    }        

    input[0].r = 1;
    input[0].i = 2;

    for(i=0; i<size_inbox; i++) 
        printf("  %g+%gi  \t ", input[i].r, input[i].i);
    printf("\n");        
    printf("\n");      

    double timer[20];

    // ********************************
    // Compute forward (Z2Z) transform
    // ********************************
    fiber_execute_z2z[my_backend].function(box_low, box_high, box_low, box_high, comm, input, input, 0, timer);

    // Output after forward
    for(i=0; i<size_outbox; i++) 
        // printf(" %g+%gi \t ", output[i].r, output[i].i);
        printf(" %g+%gi \t ", input[i].r, input[i].i);
    printf("\n");
    printf("\n");

    if (me == 1){
        int pass = 0;
        if (fabs(input[0].r + 510.0) > 1.E-11 || fabs(input[0].i - 4.0) > 1.E-11)
            pass = 1;
        for(i=1; i<16; i++)
            if (fabs(input[i].r-2) > 1.E-11 || fabs(input[i].i-4) > 1.E-11)
                pass = 1;
        for(i=16; i<32; i++)
            if (fabs(input[i].r) > 1.E-11 || fabs(input[i].i) > 1.E-11)
                pass = 1;
        if (pass){
            printf("The computed transform deviates by more than the tolerance.\n");
            MPI_Abort(comm, 1);
        }
    }

    if(me == 1)
        printf("\t\t%s library computed a correct forward C2C 3-D transform \n \n", backends[my_backend] );


    // ********************************
    // Compute backward (Z2Z) transform
    // ********************************
    // for(i=0; i<size_inbox; i++){
    //     input[i].r = 0.0;
    //     input[i].i = 0.0;
    // }        

    // fiber_execute_z2z[my_backend].function(box_low, box_high, box_low, box_high, comm, input, input, 1, timer);
    fiber_execute_z2z[8].function(box_low, box_high, box_low, box_high, comm, input, input, 1, timer);

    // Output after backward
    for(i=0; i<size_inbox; i++) {
        // Scaling
        input[i].r /= fftsize;
        input[i].i /= fftsize;
        printf(" %g+%gi  \t ", input[i].r, input[i].i);
    }        
    printf("\n");
    printf("\n");

    // Error computation after backward
    double err = 0.0;
    for(i=1; i<size_inbox; i++){
        // printf("%g+%g ----- %d \n" , input[i].r, input[i].i, i);

        if (fabs(input[i].r - (double) i ) > err){
            // printf("error at %g+%g ----- %d \n" , input[i].r, input[i].i, i);
            err = fabs(input[i].r - (double) i - 1);
        }

        if( fabs(input[i].i) > 1.E-11 )
             err = fmax( fabs(input[0].i) , err);
        // printf("%g \n" , err);
    }

    err = fmax( fabs(input[0].r - (double) 1) , err);
    err = fmax( fabs(input[0].i - (double) 2) , err);

    // Print errors
    printf("%s: rank [%d] computed error |X - ifft(fft(X)) |_{max}: %1.6le\n", backends[my_backend], me, err);

    // Data deallocation 
    free(input);
    free(output);

    MPI_Finalize();

    return 0;
}
