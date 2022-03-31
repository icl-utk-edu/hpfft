# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os
from spack import *

FFT_LIBS = (
    'heffte',
    'fftw',
    'ffte',
    'accfft',
)


class Fiber(CMakePackage, CudaPackage):
    """High Performance Fast Fourier Transform benchmark."""

    homepage = "https://fiber.icl.utk.edu/"
    git      = "https://github.com/icl-utk-edu/fiber"

    maintainers = ['G-Ragghianti', 'luszczek']

    version('master', branch='master')

    variant('fft', default='none', values=FFT_LIBS,
            multi=True, description='Supported FFT libraries')
    #variant('onlycomplex', default=False)

    depends_on('mpi')

    # Iterate over the fftlib variant for simple dependencies
    for fft in FFT_LIBS:
        depends_on(fft, when='fft={0}'.format(fft))

    depends_on('heffte+fftw+cuda', when='fft=heffte')

    def cmake_args(self):
        args = []
        for fft in self.spec.variants['fft'].value:
            args.extend(["-DFIBER_ENABLE_{0}=ON".format(fft.upper())])
            # These don't seem to be required within spack 
            #args.extend("-DFIBER_FFT_INCLUDE_DIRS={0}".format(self.spec[fft].prefix.include))
            #args.extend("-DFIBER_FFT_LIB_DIRS={0}".format(self.spec[fft].prefix.lib))
        return args

    def install(self, spec, prefix):
        mkdirp(prefix.lib)
        mkdirp(prefix.bin)
        # Copy the generated binaries to the installation prefix
        with working_dir(os.path.join(self.build_directory, 'benchmarks')):
            install('../libfiber.so', prefix.lib)
            for f in os.listdir('.'):
                if os.path.isfile(f) and os.access(f, os.X_OK):
                    install(f, prefix.bin)
