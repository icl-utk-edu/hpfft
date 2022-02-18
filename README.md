![FBI_banner](https://bitbucket.org/aayala32/logos/raw/de08df3c3bfd9d595587bf840f31afcb45d6019c/fiber.png)

**FFT Benchmarking Initiative**

**Innovative Computing Laboratory**

**University of Tennessee**



* * *

About
=====

The FFT Infrastructure Benchmark for Exascale Research (FIBER) provides a framework for Fast Fourier Transform (FFT) benchmarks targeting exascale computing systems. It evaluates performance and scalability of distributed FFTs on different architectures. Furthermore, it analyzes the effect on applications that directly depend on FFTs. It can also stress and test the overall network of a supercomputer, give an indication on bisection bandwidth, noise, and other network and MPI collectives limitations that are of interest to many other ECP applications.


The current harness software puts together FFT libraries supporting distributed 3-D complex-to-complex and real-to-complex FFTs.


* * *

Setting up
==========

Create a folder; e.g., `Benchmarks_FFT`, and install the FFT libraries to benchmark; or load them as modules.

~~~
-- Benchmarks_FFT
        |-- heFFTe
        |-- fftMPI
        |-- AccFFT
        |-- P3DFFT
        |-- FFTE
        |-- SWFFT
        |-- 2DECOMP&FFT
        |-- nb3dFFT
        |-- FFTW
        |-- FFTW++
~~~

Current libraries targeted by FIBER:
- CPU support: [fftMPI](https://lammps.github.io/fftmpi/), [SWFFT](https://xgitlab.cels.anl.gov/hacc/SWFFT), 
[P3DFFT](https://github.com/sdsc/p3dfft.3),
[nb3dFFT](https://gitlab.jsc.fz-juelich.de/goebbert/nb3dfft),
[2DECOMP&FFT](http://www.2decomp.org/download.html), [FFTW](http://www.fftw.org/), [FFTW++](fftwpp.sourceforge.net/)

- CPU-GPU support: [heFFTe](https://bitbucket.org/icl/heffte), [AccFFT](https://github.com/amirgholami/accfft),   [FFTE](http://www.ffte.jp/)


Compilation
===========

Next clone this repository and create  build folder, and execute the `cmake` commands.
In the following example, we install FIBER with heFFTe and fftMPI backends:

~~~
mkdir build; cd $_
build/
cmake -DFIBER_FFT_LIB_DIRS="/home/Benchmarks_FFT/fftmpi/src;/home/heffte/build/lib"
-DFIBER_FFT_INCLUDE_DIRS="/home/Benchmarks_FFT/fftmpi/src;/home/heffte/build/include"
-DFIBER_ENABLE_HEFFTE=ON -DFIBER_ENABLE_FFTMPI=ON
-DMPI_DIR=/sw/openmpi/4.0.0/ .. 
make -j
~~~

List the `lib` and `include` folders of libraries to test, respectively, in `FIBER_FFT_LIB_DIRS` and `FIBER_FFT_INCLUDE_DIRS`.

Testing integration
===================

Run tests as follows:
~~~
cd build/benchmarks
mpirun -n 2 ./test3D_CPU_C2C <library>
mpirun -n 2 ./test3D_CPU_R2C <library>
~~~

If FIBER was build linked to GPU enabled libraries:
~~~
cd build/benchmarks
mpirun -n 2 ./test3D_GPU_C2C <gpu_library>
mpirun -n 2 ./test3D_GPU_R2C <gpu_library>
~~~

Running benchmarks
===================
~~~
cd build/benchmarks
mpirun -n $NUM_RANKS ./test3D_C2C -lib <library> -backend <1D_backend> -size <nx> <ny> <nz> -pgrid <p> <q>
~~~


where `library` has to be replaced by one of the nine available libraries, provided user has it installed.
Once a parallel FFT library has been correctly integrated to heFFTe, running these benchmarks should report a correct validation output.


Documentation
=============

* Installation and a Doxygen documentation will be available shortly.

* * *

Getting Help
============

For assistance with the FIBER project, email *fiber@icl.utk.edu* or start a GitHub issue. 

Contributions are very welcome, please create a pull request.

Resources
=========


* Visit the [FIBER website](http://icl.utk.edu/fiber/) for more information about the HeFFTe project.
* Visit the [ECP website](https://exascaleproject.org) to find out more about the DOE Exascale Computing Initiative.

* * *

Acknowledgments
===============

This research was supported by the United States Exascale Computing Project.

* * *

License
=======

    Copyright (c) 2022, University of Tennessee
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
          notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
          notice, this list of conditions and the following disclaimer in the
          documentation and/or other materials provided with the distribution.
        * Neither the name of the University of Tennessee nor the
          names of its contributors may be used to endorse or promote products
          derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL UNIVERSITY OF TENNESSEE BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
