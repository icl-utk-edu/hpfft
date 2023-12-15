
#include "hpfft_swfft.h"

/* SWFFT header files */
#include "distribution.h"
#include "Dfft.hpp"
#include "Distribution.hpp"

#include <fftw3.h>

extern "C" {

void*
SWFFT_Distribution_new(MPI_Comm comm, int *dim3d, int debug) {
    return new hacc::Distribution(comm, dim3d, debug);
}

void
SWFFT_Distribution_delete(void *dist) {
    delete static_cast<hacc::Distribution*>(dist);
}

void*
SWFFT_Dfft_new(void * dist) {
    return new hacc::Dfft(*static_cast<hacc::Distribution*>(dist));
}

void
SWFFT_Dfft_delete(void *dfft) {
    delete static_cast<hacc::Dfft*>(dfft);
}

void
SWFFT_makePlans(void *dfft, void *forward_output, void *forward_scratch, void *backward_input, void *backward_scratch, unsigned int flags) {
    hacc::Dfft *this_ = static_cast<hacc::Dfft*>(dfft);
    this_->makePlans((std::complex<double>*)forward_output, (std::complex<double>*)forward_scratch, (std::complex<double>*)backward_input, (std::complex<double>*)backward_scratch, flags);
}

void
SWFFT_forward(void *dfft, void *input) {
    hacc::Dfft *this_ = static_cast<hacc::Dfft*>(dfft);
    this_->forward((std::complex<double>*)input);
}

void
SWFFT_backward(void *dfft, void *output) {
    hacc::Dfft *this_ = static_cast<hacc::Dfft*>(dfft);
    this_->backward((std::complex<double>*)output);
}

}
