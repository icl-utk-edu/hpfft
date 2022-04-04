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

    depends_on('mpi')

    @property
    def build_targets(self):
        return ['-j1']

    def edit(self, spec, prefix):
        with working_dir('src'):
            copy('Makefile.inc.x86', 'Makefile.inc')

    def install(self, spec, prefix):
        make('install', 'prefix={0}'.format(prefix))
