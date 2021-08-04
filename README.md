<!-- ![FBI_banner](https://bitbucket.org/aayala32/logos/raw/4fd193bfb8e06939c1d8ca4d3ffc389fee50f7f6/FBI_logo.png) -->

**FFT Benchmarking Initiative**

**Innovative Computing Laboratory**

**University of Tennessee**


* * *

[TOC]

* * *

About
=====

The FFT Benchmarking Initiative (FBI) provides a framework for Fast Fourier Transform (FFT) benchmarks targeting exascale computing systems. It evaluates performance and scalability of distributed FFTs on different architectures. Furthermore, it analyzes the effect on applications that directly depend on FFTs. It can also stress and test the overall network of a supercomputer, give an indication on bisection bandwidth, noise, and other network and MPI collectives limitations that are of interest to many other ECP applications.


The current harness software allows to compute 3-D complex-to-complex and real-to-complex FFTs.


* * *

Compilation & first experiment
==============================

~~~
mkdir build; cd $_
build/
    cmake ..
    make -j
~~~

*Integrated libraries*: heFFTe, FFTMPI and AccFFT.

*In progress*: P3DFFT, FFTE, SWFFT, 2DECOMP&FFT, nb3dFFT, FFTW

Running examples:
~~~
mpirun -n 2 ./test3D_CPU_C2C
mpirun -n 2 ./test3D_CPU_R2C
~~~


Documentation
=============

* Installation and a Doxygen documentation will be available shortly.

* * *

Getting Help
============

For assistance with FBI, email *aayala@icl.utk.edu* or start a GitHub issue. 

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

    Copyright (c) 2021, University of Tennessee
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
