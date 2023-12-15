/* Declaration of SWFFT funcitons for use in C code */

#ifndef HPFFT_SWFFT_H
#define HPFFT_SWFFT_H

#include <mpi.h>

#if defined(__cplusplus)
extern "C" {
#endif

void *SWFFT_Distribution_new(MPI_Comm comm, int *dim3d, int debug);
void *SWFFT_Dfft_new(void * dist);
void SWFFT_Dfft_delete(void *dfft);
void SWFFT_Distribution_delete(void *dist);
void SWFFT_makePlans(void *dfft, void *forward_output, void *forward_scratch, void *backward_input, void *backward_scratch, unsigned int flags);
void SWFFT_forward(void *dfft, void *input);
void SWFFT_backward(void *dfft, void *output);

#if defined(__cplusplus)
}
#endif

#endif  /* HPFFT_SWFFT_H */
