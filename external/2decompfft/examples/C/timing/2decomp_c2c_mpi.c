#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>
#include <decomp_2d_iface.h>
#include <string.h>
#ifdef WITH_LOGGER
#include <Log.h>
#else
#define log_print(_a, _b, _c, ...)
#endif


int main(int argc, char **argv){
  int rank  = 0;
  int size  = 1;
  int niter = 20;
  double *titer = NULL;
  double *tfwd  = NULL;
  double *tbwd  = NULL;
  double *max_tfwd = NULL;
  double *max_tbwd = NULL;
  double *min_tfwd = NULL;
  double *min_tbwd = NULL;
  int nx = 8;
  int ny = 8;
  int nz = 8;
  int p_row, p_col;
  int i, j, k;
  int physical = PHYSICAL_IN_X;
  double scaling = 0.0;
  MPI_Comm comm;

  if ( argc < 4 ) {
    fprintf ( stderr, "Error, the problem size must be given in the commandline\n" );
    return -1;
  }

  nx = atoi(argv[1]);
  ny = atoi(argv[2]);
  nz = atoi(argv[3]);

#ifdef WITH_LOGGER
  if ( argc > 4 ) {
    if ( ! strcmp( argv[4], "-v" ) ) {
      int lvl = (argc > 5) ? atoi( argv[4] ) : info;

      log_on( algo, lvl );
      log_on( alloc, details );
    //log_on( comp, details );
    }
  }
#endif
  
//nx = 32;

  p_row = 0; // The library will handle it
  p_col = 0; // The library will handle it

  // FORCE for testing
  physical = PHYSICAL_IN_X;

  scaling = 1 / (  (double)nx * ny * nz );

  MPI_Init(&argc, &argv);
  MPI_Comm_dup ( MPI_COMM_WORLD, &comm );

  log_print( algo, info, "Initialize decomp_2d_init\n" );
  decomp_2d_init( nx, ny, nz, p_row, p_col );

  MPI_Barrier( comm );
  log_print( algo, info, "Initialize decomp_2d_fft_init\n" );
  decomp_2d_fft_init(physical);

  MPI_Comm_rank ( comm, &rank );
  MPI_Comm_size ( comm, &size );

  titer     = (double*) malloc ( niter * sizeof(double));
  tfwd      = (double*) malloc ( niter * sizeof(double));
  tbwd      = (double*) malloc ( niter * sizeof(double));
  max_tfwd  = (double*) malloc ( niter * sizeof(double));
  max_tbwd  = (double*) malloc ( niter * sizeof(double));
  min_tfwd  = (double*) malloc ( niter * sizeof(double));
  min_tbwd  = (double*) malloc ( niter * sizeof(double));

  /*get local data size and allocate*/
  int lxsize = 0;
  int lysize = 0;
  int lzsize = 0;
  int dim = 0;
  double complex *data_in  = NULL;
  double complex *data_out = NULL;
  double complex *ref = NULL;

  decomp_2d_get_local_sizes(physical, &lxsize, &lysize, &lzsize);
  log_print( alloc, details, "RETURNED lxsize=%d\tlysize=%d\tlzsize=%d\n",
      lxsize, lysize, lzsize );

//lxsize = nx;
//lysize = ny;
//lzsize = nz;
//log_print( alloc, details, "OVERWRITTEN lxsize=%d\tlysize=%d\tlzsize=%d\n",
//    lxsize, lysize, lzsize );

  size_t localsize = lxsize * lysize * lzsize * sizeof(double complex);
  log_print( alloc, info, "Local size = %zu\n", localsize );

//MPI_Barrier( comm );
//goto finalize;

  data_in = (double complex*) malloc ( localsize );
  if ( ! data_in ) {
    fprintf( stderr, "Error, cannot allocate data of size %zu\n",
      localsize);
    goto finalize;
  }
  data_out = (double complex*) malloc ( localsize );
  if ( ! data_out ) {
    fprintf( stderr, "Error, cannot allocate data of size %zu\n",
      localsize);
    goto finalize;
  }
  ref  = (double complex*) malloc ( localsize );
  if ( ! ref) {
    fprintf( stderr, "Error, cannot allocate ref of size %zu\n",
      localsize);
    goto finalize;
  }

#if 1
  /*initialize rin to some functionmy_func(x,y,z) */
  for (i = 0; i < lxsize; ++i) {
    double tmp_x = ((double)i) / nx;
    for (j = 0; j < lysize; ++j) {
      double tmp_y = ((double)j) / ny;
      for (k = 0; k < lzsize; ++k) {
        double tmp_z = ((double)k) / nz; 
      //data[(i*M + j) * N + k] = rand() / ((double) RAND_MAX ) + 0 * I;
      //ref [(i*M + j) * N + k] = data[(i*M + j) * N + k];
      //double complex entry = 3 + 2*I;
        double complex entry = tmp_x * tmp_y * tmp_z;
        data_in[(i*lysize + j) * lzsize + k] = entry;
        ref [(i*lysize + j) * lzsize + k]    = data_in[(i*lysize + j) * lzsize + k];

        log_print( comp, details, "ref[%d]=(%.2e,%.2e)\n",
          (i*lysize + j) * lzsize + k,
          creal( ref [(i*lysize + j) * lzsize + k] ),
          cimag( ref [(i*lysize + j) * lzsize + k] ) );

      }
    }
  }

  MPI_Barrier( comm );

  /*compute transforms as many times as desired*/
  double tmp    = 0.0;
  double ti     = 0.0;
  double tn     = 0.0;
  double max_tn = 0.0;
  
  // Warmup
  log_print( algo, info, "Warmup\n" );
  log_print( algo, info, "FORWARD\n" );
  decomp_2d_fft_3d_c2c( lxsize, lysize, lzsize, data_in, 
      lxsize, lysize, lzsize, data_out, FORWARD );

  MPI_Barrier( comm );
  log_print( comp, details, "data_out[%d]=(%.2e,%.2e)\n",
      0, creal( data_out [0] ), cimag( data_out [0] ) );
//goto finalize;

  log_print( algo, info, "BACKWARD\n" );
  decomp_2d_fft_3d_c2c( lxsize, lysize, lzsize, data_out, 
      lxsize, lysize, lzsize, data_in, BACKWARD );
  MPI_Barrier( comm );
  log_print( comp, details, "data_in[%d]=(%.2e,%.2e)\n",
      0, creal( data_out [0] ), cimag( data_out [0] ) );

  // Normalization
  log_print( algo, info, "NORMALIZATION\n" );
  for (i = 0; i < lxsize; ++i)
    for (j = 0; j < lysize; ++j)
      for (k = 0; k < lzsize; ++k) {
        data_in[(i*lysize + j) * lzsize + k] *= scaling;
      }
  log_print( comp, details, "AFTER Normalization: data_in[%d]=(%.2e,%.2e)\n",
      0, creal( data_out [0] ), cimag( data_out [0] ) );

  log_print( algo, info, "Main loop\n" );
  MPI_Barrier( comm );
  for ( int iter = 0; iter < niter; ++iter ) {
    ti = MPI_Wtime();
    decomp_2d_fft_3d_c2c( lxsize, lysize, lzsize, data_in, 
        lxsize, lysize, lzsize, data_out, FORWARD );
    tfwd[iter] = MPI_Wtime() - ti;

    tmp = MPI_Wtime();
    decomp_2d_fft_3d_c2c( lxsize, lysize, lzsize, data_out, 
        lxsize, lysize, lzsize, data_in, BACKWARD );
    tbwd[iter] = MPI_Wtime() - tmp;
    titer[iter] = MPI_Wtime() - ti;
    tn += titer[iter];

    // Normalization
    for (i = 0; i < lxsize; ++i)
      for (j = 0; j < lysize; ++j)
        for (k = 0; k < lzsize; ++k) {
          data_in[(i*lysize + j) * lzsize + k] *= scaling;
        }
    
    log_print( algo, info, "iter %d: fwd=%.2e\tbwd=%.2e\n",
      iter, tfwd[iter], tbwd[iter] );
    MPI_Barrier ( comm );
  }

  // Reduce operation to compute min, max and avg
  MPI_Reduce ( &tn, &max_tn, 1, MPI_DOUBLE, MPI_MAX, 0, comm );

  for ( int i = 0; i < niter; ++i ) {
    MPI_Reduce ( tfwd + i, max_tfwd + i, 1, MPI_DOUBLE, MPI_MAX, 0, comm );
    MPI_Reduce ( tbwd + i, max_tbwd + i, 1, MPI_DOUBLE, MPI_MAX, 0, comm );
    MPI_Reduce ( tfwd + i, min_tfwd + i, 1, MPI_DOUBLE, MPI_MIN, 0, comm );
    MPI_Reduce ( tbwd + i, min_tbwd + i, 1, MPI_DOUBLE, MPI_MIN, 0, comm );
  }

  // Compute local error
  double err = 0.0;
  double ierr = 0.0;
  double max_ierr = 0.0;
  for (i = 0; i < lxsize; ++i)
    for (j = 0; j < lysize; ++j)
      for (k = 0; k < lzsize; ++k) {
        double complex cref_tmp   = ref[(i*lysize + j) * lzsize + k];
        double complex cdata_tmp  = data_in[(i*lysize + j) * lzsize + k];
        log_print( comp, details, "data_in[%d]=(%.2e,%.2e)\n",
          (i*lysize + j) * lzsize + k,
          creal( cdata_tmp ),
          cimag( cdata_tmp ) );
        double tmp_r = creal(cref_tmp) - creal(cdata_tmp);
        double tmp_i = cimag(cref_tmp) - cimag(cdata_tmp);
        err += sqrt(tmp_r*tmp_r + tmp_i * tmp_i);
        double tmp_real = sqrt(tmp_r * tmp_r);
        if ( tmp_real > ierr )
          ierr = tmp_real;
      }

  MPI_Reduce ( &ierr, &max_ierr, 1, MPI_DOUBLE, MPI_MAX, 0, comm );
/*
     n1 = real(nx,mytype) * real(ny,mytype) * real(nz,mytype)
     n1 = n1 ** (1._mytype/3._mytype)
     ! 5n*log(n) flops per 1D FFT of size n using Cooley-Tukey algorithm
     flops = 5._mytype * n1 * log(n1) / log(2.0_mytype)
     ! 3 sets of 1D FFTs for 3 directions, each having n^2 1D FFTs
     flops = flops * 3._mytype * n1**2  
     flops = 2._mytype * flops / ((t1+t3)/real(NTEST,mytype))
     write(*,*) 'GFLOPS : ', flops / 1000._mytype**3
*/
  // Computation of the flops
  double fwdAdd, fwdMul, fwdFma;
  double bwdAdd, bwdMul, bwdFma;
  double n1 = pow( (double)nx * ny * nz, 1./3. );
  double flops = 5. * n1 * log(n1) / log(2);
  flops *= 3. * n1 * n1;

  if ( ! rank ) {
    printf ( "\n========================\nInput:" );
    for ( int i = 1; i < argc; ++i)
      printf ( " %s", argv[i] );
    printf ( "\n" );
    printf ( "2decomp&FFT test results\n---------\n"
        "Problem size:          %dx%dx%d\n"
        "MPI ranks:             %d\n"
        "Execution(niter=%3d):  %.8fs\n"
        "GFlops:                %.3e\n"
        "Error:                 %.2e\n",
        nx, ny, nz,
        size,
        niter, max_tn,
        (( 2 * flops ) / (max_tn / niter)) / ( 1024l * 1024l * 1024l ),
        max_ierr );

    // Print some details
    printf ( "Fwd details: max_i=[" );
    for ( int i = 0; i < niter; ++i )
      printf ( "%.6e ", max_tfwd[i] );
    printf ( "]\n" );
    printf ( "Fwd details: min_i=[" );
    for ( int i = 0; i < niter; ++i )
      printf ( "%.6e ", min_tfwd[i] );
    printf ( "]\n" );
    printf ( "Bwd details: max_i=[" );
    for ( int i = 0; i < niter; ++i )
      printf ( "%.6e ", max_tbwd[i] );
    printf ( "]\n" );
    printf ( "Bwd details: min_i=[" );
    for ( int i = 0; i < niter; ++i )
      printf ( "%.6e ", min_tbwd[i] );
    printf ( "]\n" );
  }


finalize:
//decomp_2d_fft_finalize();
//decomp_2d_finalize();
#endif

//if ( data_in ) free ( data_in );
//if ( data_out ) free ( data_out );
//if ( ref ) free ( ref );

  if ( titer ) free ( titer );
  if ( tfwd ) free ( tfwd );
  if ( tbwd ) free ( tbwd );
  if ( max_tfwd ) free ( max_tfwd );
  if ( max_tbwd ) free ( max_tbwd );
  if ( min_tfwd ) free ( min_tfwd );
  if ( min_tbwd ) free ( min_tbwd );

  MPI_Comm_free ( &comm );
  MPI_Finalize();

  return 0;
}
