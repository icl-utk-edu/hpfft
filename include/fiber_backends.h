#ifndef HPFFT_BACKENDS_H
#define HPFFT_BACKENDS_H

#include <mpi.h>

#include "hpfft_backend_heffte.h"
#include "hpfft_backend_fftmpi.h"
#include "hpfft_backend_accfft.h"
#include "hpfft_backend_p3dfft.h"
#include "hpfft_backend_ffte.h"
#include "hpfft_backend_swfft.h"
#include "hpfft_backend_decomp2d.h"
#include "hpfft_backend_nb3dfft.h"
#include "hpfft_backend_fftw.h"
#include "hpfft_backend_fftadvmpi.h"
#include "hpfft_backend_fftwpp.h"

// Available backends
char backends[][20] = {"heFFTe", "FFTMPI", "AccFFT", "P3DFFT", "FFTE", "SWFFT", "2DECOMP&FFT", "nb3dFFT", "FFTW", "fftadvmpi", "FFTW++"};
int n_backends = 11;

enum backend{heffte, fftmpi, accfft, p3dfft, ffte, swfft, decomp2d, nb3dfft, fftw, fftadvmpi, fftwpp};

// Backends initialization
struct hpfft_init
{
    char *name;
    int (*function)( int );
};

const struct hpfft_init hpfft_initialize[] = {
    { "heffte",     &init_heffte    },
    { "fftmpi",     &init_fftmpi    },
    { "accfft",     &init_accfft    },
    { "p3dfft",     &init_p3dfft    },
    { "ffte",       &init_ffte      },
    { "swfft",      &init_swfft     },
    { "decomp2d",   &init_decomp2d  },
    { "nb3dfft",    &init_nb3dfft   },
    { "fftw",       &init_fftw      },
    { "fftadvmpi",  &init_fftadvmpi },
    { "fftwpp",     &init_fftwpp }
};

// Real-to-complex transform
struct hpfft_map_backend_d2z
{
    char *name;
    void (*function)( int const * , int const *, int const * , int const * , MPI_Comm const , double const *, void *, int *, double *);
};

const struct hpfft_map_backend_d2z hpfft_execute_d2z[] = {
    { "heffte",     &compute_d2z_heffte    },
    { "fftmpi",     &compute_d2z_fftmpi    },
    { "accfft",     &compute_d2z_accfft    },
    { "p3dfft",     &compute_d2z_p3dfft    },
    { "ffte",       &compute_d2z_ffte      },
    { "swfft",      &compute_d2z_swfft     },
    { "decomp2d",   &compute_d2z_decomp2d  },
    { "nb3dfft",    &compute_d2z_nb3dfft   },
    { "fftw",       &compute_d2z_fftw      },
    { "fftadvmpi",  &compute_d2z_fftadvmpi },
    { "fftwpp",     &compute_d2z_fftwpp }
};


// Complex-to-real transform
struct hpfft_map_backend_z2d
{
    char *name;
    void (*function)( int const * , int const *, int const * , int const * , MPI_Comm const , void const *, double *, int *, double *);
};

const struct hpfft_map_backend_z2d hpfft_execute_z2d[] = {
    { "heffte",     &compute_z2d_heffte    },
    { "fftmpi",     &compute_z2d_fftmpi    },
    { "accfft",     &compute_z2d_accfft    },
    { "p3dfft",     &compute_z2d_p3dfft    },
    { "ffte",       &compute_z2d_ffte      },
    { "swfft",      &compute_z2d_swfft     },
    { "decomp2d",   &compute_z2d_decomp2d  },
    { "nb3dfft",    &compute_z2d_nb3dfft   },
    { "fftw",       &compute_z2d_fftw      },
    { "fftadvmpi",  &compute_z2d_fftadvmpi },
    { "fftwpp",     &compute_z2d_fftwpp }
};


// Complex-to-Complex transform
struct hpfft_map_backend_z2z
{
    char *name;
    void (*function)( int const * , int const *, int const * , int const * , MPI_Comm const , void const *, void *, int *, double *);
};

const struct hpfft_map_backend_z2z hpfft_execute_z2z[] = {
    { "heffte",     &compute_z2z_heffte    },
    { "fftmpi",     &compute_z2z_fftmpi    },
    { "accfft",     &compute_z2z_accfft    },
    { "p3dfft",     &compute_z2z_p3dfft    },
    { "ffte",       &compute_z2z_ffte      },
    { "swfft",      &compute_z2z_swfft     },
    { "decomp2d",   &compute_z2z_decomp2d  },
    { "nb3dfft",    &compute_z2z_nb3dfft   },
    { "fftw",       &compute_z2z_fftw      },
    { "fftadvmpi",  &compute_z2z_fftadvmpi },
    { "fftwpp",     &compute_z2z_fftwpp }
};

#endif  //! HPFFT_BACKENDS_H_
