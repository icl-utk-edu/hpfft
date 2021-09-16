#ifndef FIBER_BACKENDS_H_
#define FIBER_BACKENDS_H_


#include <mpi.h>

#include "fiber_backend_heffte.h"
#include "fiber_backend_fftmpi.h"
#include "fiber_backend_accfft.h"
// #include "fiber_backend_p3dfft.h"
// #include "fiber_backend_ffte.h"
// #include "fiber_backend_swfft.h"
// #include "fiber_backend_decomp2d.h"
// #include "fiber_backend_nb3dfft.h"
// #include "fiber_backend_fftw.h"

// Available backends
char backends[][20] = {"heFFTe", "FFTMPI", "AccFFT", "P3DFFT", "FFTE", "SWFFT", "2DECOMP&FFT", "nb3dFFT", "FFTW"};
int n_backends = 9;

// enum backend{heffte, fftmpi, accfft, p3dfft, ffte, swfft, decomp2d, nb3dfft, fftw};
enum backend{heffte, fftmpi, accfft};

// Real-to-complex transform
struct fiber_map_backend_d2z
{
    char *name;
    void (*function)( int const * , int const *, int const * , int const * , MPI_Comm const , double const *, void *, int, double *);
};

const struct fiber_map_backend_d2z fiber_execute_d2z[] = {
    { "heffte",   &compute_d2z_heffte   }
   ,{ "fftmpi",   &compute_d2z_fftmpi   }
   ,{ "accfft",   &compute_d2z_accfft   }
//    ,{ "p3dfft",   &compute_d2z_p3dfft   }
//    ,{ "ffte",     &compute_d2z_ffte     }
//    ,{ "swfft",    &compute_d2z_swfft    }
//    ,{ "decomp2d", &compute_d2z_decomp2d }
//    ,{ "nb3dfft",  &compute_d2z_nb3dfft  }
//    ,{ "fftw",     &compute_d2z_fftw     }
};


// Complex-to-real transform
struct fiber_map_backend_z2d
{
    char *name;
    void (*function)( int const * , int const *, int const * , int const * , MPI_Comm const , double const *, void *, double *);
};

const struct fiber_map_backend_z2d fiber_execute_z2d[] = {
    { "heffte",   &compute_z2d_heffte   }
   ,{ "fftmpi",   &compute_z2d_fftmpi   }
   ,{ "accfft",   &compute_z2d_accfft   }
//    ,{ "p3dfft",   &compute_z2d_p3dfft   }
//    ,{ "ffte",     &compute_z2d_ffte     }
//    ,{ "swfft",    &compute_z2d_swfft    }
//    ,{ "decomp2d", &compute_z2d_decomp2d }
//    ,{ "nb3dfft",  &compute_z2d_nb3dfft  }
//    ,{ "fftw",     &compute_z2d_fftw     }
};


struct fiber_map_backend_z2z
{
    char *name;
    void (*function)( int const * , int const *, int const * , int const * , MPI_Comm const , void const *, void *, int, double *);
};

const struct fiber_map_backend_z2z fiber_execute_z2z[] = {
    { "heffte",   &compute_z2z_heffte   }
   ,{ "fftmpi",   &compute_z2z_fftmpi   }
   ,{ "accfft",   &compute_z2z_accfft   }
//    ,{ "p3dfft",   &compute_z2z_p3dfft   }
//    ,{ "ffte",     &compute_z2z_ffte     }
//    ,{ "swfft",    &compute_z2z_swfft    }
//    ,{ "decomp2d", &compute_z2z_decomp2d }
//    ,{ "nb3dfft",  &compute_z2z_nb3dfft  }
//    ,{ "fftw",     &compute_z2z_fftw     }
};

#endif  //! FIBER_BACKENDS_H_ 
