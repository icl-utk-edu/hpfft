#!/bin/bash --login

MPI=$1
FFT=$2

source /etc/profile
set +x
set -e
trap 'echo "# $BASH_COMMAND"' DEBUG
shopt -s expand_aliases

export HOME=`pwd`
git clone https://github.com/spack/spack ../spack || true
cp spack/modules.yaml spack/packages.yaml spack/upstreams.yaml ../spack/etc/spack/
source ../spack/share/spack/setup-env.sh
module load gcc@7
spack compiler find
spack repo add `pwd`/spack/ || true
spack uninstall -a -y --dependents $FFT fiber || true
spack install --fresh cmake fftw
spack dev-build --fresh fiber@master fft=$FFT ^$MPI

# Run the tests
spack load fiber
mpirun -n 2 test3D_C2C -lib $FFT -backend fftw -size 4 4 4 -pgrid 1 2
