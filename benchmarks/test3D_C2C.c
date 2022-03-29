/*
* Harness software for benchmarking parallel FFT libraries.
* Benchmark: C2C 3-D transform, double precision
* Autor: Alan Ayala - ICL, UTK.
! USAGE:
As input, give the library you want to benchmark, the FFT size (nx,ny,nz), and the processor partition pxq.
    make clean; make -j;
    mpirun -n $NUM_RANKS ./test3D_C2C -lib <library> -backend <1D_backend> -size <nx> <ny> <nz> -pgrid <p> <q>
*/

#include "fiber_backends.h"
#include "fiber_utils.h"


int main(int argc, char** argv){
   
   // MPI initialization
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int me;
    MPI_Comm_rank(comm, &me);

    int num_ranks;
    MPI_Comm_size(comm, &num_ranks);

    // Parameters
    double timer[20];    
    int backend_options[n_options];
    char lib_name[30];
    char lib_1d_backend[30] = {"fftw"}; // FFTW is the default backend for most parallel FFT libraries
    backend_options[option_niter] = 10; // default number of iterations

    fiber_parse_options(argc, argv, backend_options, lib_name, lib_1d_backend);
    backend_options[option_fft_op] = 0;  // Fixed to benchmark a forward FFT computation
    backend_options[option_backend] = fiber_get_1d_backend(lib_1d_backend);

    if(me==0){
        printf("----------------------------------------------- \n");
        printf("Benchmarking %s library using %s backend \n", lib_name, lib_1d_backend);
        printf("----------------------------------------------- \n");
    }        
    int my_backend  = fiber_get_backend(lib_name);

    // Global grid
    int nx, ny, nz ;
    nx = backend_options[option_nx];
    ny = backend_options[option_ny];
    nz = backend_options[option_nz];

    int fftsize = nx*ny*nz;

    // select grid of processors with MPI built-in function
    // TODO: read proc_grid from CL
    int proc_grid[2] = {0,0};
    proc_grid[0] = backend_options[option_grid_p];
    proc_grid[1] = backend_options[option_grid_q];

    if (proc_grid[0]==0 || proc_grid[1]==0){
        MPI_Dims_create(num_ranks, 2, proc_grid);
        backend_options[option_grid_p] = proc_grid[0];
        backend_options[option_grid_q] = proc_grid[1];
    }
    // For slab decomposition use 1x1xnum_ranks
    // proc_grid[0] = 1;
    // proc_grid[1] = num_ranks;
    
    // create cartesian grid
    int periods[2]  = {0,0};
    int mpi_reorder = 1;

    MPI_Comm cart_comm;
    MPI_Cart_create(comm, 2, proc_grid, periods, mpi_reorder, &cart_comm);

    int proc_rank;
    MPI_Comm_rank(cart_comm, &proc_rank);
    int proc_coords[2];
    MPI_Cart_coords(cart_comm, proc_rank, 2, proc_coords);

    // definition of local FFT grid (box)
    int box_low[3]  = {0, 0, 0};
    int box_high[3] = {0, 0, 0};

    // fast dimension
    box_low[0]  =  0;
    box_high[0] =  nx-1;
    // mid dimension
    box_low[1]  =  proc_coords[0] * (ny / proc_grid[0]);
    box_high[1] =  (proc_coords[0]+1) * (ny / proc_grid[0]) - 1; 
    // slow dimension
    box_low[2]  =  proc_coords[1] * (nz / proc_grid[1]);
    box_high[2] =  (proc_coords[1]+1) * (nz / proc_grid[1]) - 1; 

    int local_fft_size = (box_high[0]-box_low[0]+1)*(box_high[1]-box_low[1]+1)*(box_high[2]-box_low[2]+1);

   printf("[%d] proc_grid 1x%dx%d - [New MPI process %d] I am located at (%d, %d).   %d-%d-%d    %d-%d-%d   local_size={%d} \n", me, proc_grid[0], proc_grid[1],  proc_rank, proc_coords[0],proc_coords[1], \
                                        box_low[0], box_low[1], box_low[2], box_high[0], box_high[1], box_high[2], local_fft_size);

    fiber_complex *input  = calloc(local_fft_size, sizeof(fiber_complex));
    fiber_complex *output = calloc(local_fft_size, sizeof(fiber_complex));

    // Data Initialization
    int i;
    for(i=0; i<local_fft_size; i++){
        input[i].r = (double) i;
        input[i].i = (double) 0*i;
    }        
    input[0].r = 1;
    input[0].i = 2;

    // for(i=0; i<local_fft_size; i++) 
    //     printf("  %g+%gi, ", input[i].r, input[i].i);
    // printf("\n");        
    // printf("\n");      

    
    // Initialize Library
    int init_option = 1; // 
    if(fiber_initialize[my_backend].function(init_option) == 0){
        if(me==0)
            printf("[%s] successfully initialized.\n", lib_name);
    }
    else {
        if(me==0)
            printf("[%s] failed to initialize.\n", lib_name);
    }
    // ********************************
    // Compute forward (Z2Z) transform
    // ********************************
    fiber_execute_z2z[my_backend].function(box_low, box_high, box_low, box_high, comm, input, output, backend_options, timer);

    if (me==0){
        printf("Time for FFT plan  = %10.6g \n", timer[0]);
        printf("Time for execution = %10.6g \n", timer[1]/backend_options[option_niter]);
    }

    // Output after forward
    // for(i=0; i<local_fft_size; i++) 
    //     printf(" %g+%gi, ", output[i].r, output[i].i);
    // printf("\n");
    // printf("\n");

    MPI_Barrier(MPI_COMM_WORLD);

    backend_options[option_fft_op] = 1; // Setting to backward for error validation
    // Using heffte for backward FFT:
    // fiber_execute_z2z[0].function(box_low, box_high, box_low, box_high, comm, output, input, backend_options, timer);
    fiber_execute_z2z[my_backend].function(box_low, box_high, box_low, box_high, comm, output, input, backend_options, timer);

    // Output after backward
    for(i=0; i<local_fft_size; i++) {
        input[i].r /= fftsize;
        input[i].i /= fftsize;
        // printf(" %g+%gi, ", input[i].r, input[i].i);
    }        
    // printf("\n");
    // printf("\n");

    // Error computation after backward
    double err = 0.0;
    for(i=1; i<local_fft_size; i++){
        if (fabs(input[i].r - (double) i ) > err)
            err = fabs(input[i].r - (double) i);

        if( fabs(input[i].i) > 1.E-11 )
            err = fmax( fabs(input[i].i) , err);
    }
    
    err = fmax( fabs(input[0].r - (double) 1) , err);
    err = fmax( fabs(input[0].i - (double) 2) , err);

    double max_err = 0.0;
    MPI_Allreduce(&err, &max_err, 1, MPI_DOUBLE, MPI_MAX, comm);

    // Print error
    // printf("%s: rank [%d] computed error |X - ifft(fft(X)) |_{max}: %1.6le\n", backends[my_backend], me, err);

    if (me==0)
        printf("%s: |X - ifft(fft(X)) |_{max}: %1.6le\n", backends[my_backend], err);

    // Data deallocation 
    free(input);
    free(output);

    MPI_Finalize();

    return 0;
}
