#ifndef DECOMP_2D_IFACE_H
#define DECOMP_2D_IFACE_H

#include <complex.h>

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif

enum {
  PHYSICAL_IN_X=0,
  PHYSICAL_IN_Z=1,
};

enum {
  BACKWARD=0,
  FORWARD=1,
};

extern void decomp_2d_init(int nx, int ny, int nz, int p_row, int p_col);

extern void decomp_2d_fft_init(int physical);

//extern void decomp_2d_get_local_start(int dim, int *xstart, int *ystart, int *zstart);

extern void decomp_2d_get_local_sizes(int physical, int *size_0, int *size_1, int *size_2);

// Direction is either DECOMP_2D_FFT_FORWARD or DECOMP_2D_FFT_BACKWARD
//extern void decomp_2d_fft_3d_c2c( double complex  *data_in,
//                                  double complex  *data_out,
//                                  int             direction);
extern void decomp_2d_fft_3d_c2c( int nx_in, int ny_in, int nz_in,
                                  double complex  *data_in,
                                  int nx_out, int ny_out, int nz_out,
                                  double complex  *data_out,
                                  int             direction);

// For other datatypes like real, add interfaces here

extern void decomp_2d_fft_finalize();

extern void decomp_2d_finalize();

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
}
#endif

#endif
