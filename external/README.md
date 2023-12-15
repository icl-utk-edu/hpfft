
* * *

External Interfaces
=====

HPFFT is built in C language to support most libraries in the literature. Use this folder to add C-interface for libraries that do not support wrappers by default.


* * *

AccFFT
============

[AccFFT](http://www.accfft.org/) is one of the few libraries with GPU (NVIDIA) support. However, it does not provide a C-interface; hence, we created C-wrappers for FFT functions for both CPU and GPU kernels. To use the C-interface, do as follows:

* [Clone AccFFT](https://github.com/amirgholami/accfft)
* Use the files in `external/accfft/`: add the `.cpp` files into the cloned AccFFT `src/` folder, and the `.h` files into the `include/` folder.
* Overwrite the AccFFT `CMakeLists` using `external/accfft/CMakeLists.txt`. Modify lines `129-131` according to the available GPU hardware.

2DECOMP&FFT
============

The [2DECOMP&FFT](http://www.2decomp.org/) library is written in Fortran and no longer being maintained. Creating a C-interface requires some extra work, which we elaborate in [2Decomp C-Interface](https://github.com/cayrols/2decompFFT_c_iface/tree/3957f639890c2e9a8113d27b49c5cf75eb060d95).
