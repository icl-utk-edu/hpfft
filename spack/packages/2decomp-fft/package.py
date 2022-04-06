# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os
from spack import *


class _2decompFft(MakefilePackage):
    """The 2DECOMP&FFT library is a software framework in Fortran to build
large-scale parallel applications. It is designed for applications using
three-dimensional structured mesh and spatially implicit numerical
algorithms."""

    homepage = "http://www.2decomp.org/"
    url      = "http://www.2decomp.org/download/2decomp_fft-1.5.847.tar.gz"

    maintainers = ['G-Ragghianti']

    version('1.5.847', sha256='b137d7cf9b771de0a174d1e6c4ff0e48b4a84b51ff6c62140b4f522092e6784f')

    variant(
        'backend', default='generic', description='FFT backend',
        values=('generic', 'fftw3', 'mkl'), multi=False
    )

    depends_on('mpi')
    depends_on('mkl', when='backend=mkl')
    depends_on('fftw', when='backend=fftw3')

    #@property
    #def build_targets(self):
    #    return ['-j1']

    def edit(self, spec, prefix):
        makefile = FileFilter('src/Makefile')
        makefile.filter('include Makefile.inc', '')
        #with working_dir('src'):
        #    touch('Makefile.inc')

    def build(self, spec, prefix):
        make_args = ['-j1', # This avoids a race condition in the Makefile
                     'FFT={0}'.format(spec.variants['backend'].value),
                     'LDFLAGS=-O3',
                     'CC=mpicc',
                     'F90=mpif90',
                     'IFORT=gfortran',
                     'OPTIONS=-DDOUBLE_PREC'
                    ]
        f90flags = 'F90FLAGS=-O3 -fcray-pointer -cpp'
        if 'backend=fftw3' in spec:
            f90flags += ' -I{0}'.format(spec['fftw'].prefix.include)
        if 'backend=mkl' in spec:
            f90flags += ' -I{0}'.format(spec['mkl'].prefix.include)
            make_args.append('MKL_PATH={0}/mkl'.format(spec['mkl'].prefix))
        make_args.append(f90flags)
        make('lib', *make_args)

    def install(self, spec, prefix):
        # Avoid using the default "make install" because it
        # requires building the example codes.
        mkdirp(prefix.include, prefix.lib)
        install('include/*.mod', prefix.include)
        install('lib/lib*.a', prefix.lib)
