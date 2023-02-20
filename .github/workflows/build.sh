#!/bin/bash --login

MPI=$1
FFT=$2

source /etc/profile
set +x
set -e
trap 'echo "# $BASH_COMMAND"' DEBUG
shopt -s expand_aliases

# Set up dependencies with spack
OLDHOME=$HOME
export HOME=`pwd`
git clone https://github.com/spack/spack $OLDHOME/spack || true
source $OLDHOME/spack/share/spack/setup-env.sh
module load gcc@7
spack compiler find
spack repo add `pwd`/spack/ || true
spack uninstall -a -y --dependents $FFT || true
spack env activate --temp
spack add cmake cuda fftw $MPI $FFT
spack install --fail-fast
spack load --first cmake cuda fftw $MPI $FFT

# Build the project
mkdir -p build && cd build
FFT_DIR=`spack location -i $FFT`
FFTW_DIR=`spack location -i fftw`
MPI_DIR=`spack location -i $MPI`
export CPATH=$FFT_DIR/include:$FFTW_DIR/include:$MPI_DIR/include
export LIBRARY_PATH=$FFT_DIR/lib:$FFTW_DIR/lib:$MPI_DIR/lib
export LD_LIBRARY_PATH=$LIBRARY_PATH
LIBNAME=$FFT
cmake -DFIBER_ENABLE_${LIBNAME^^}=ON -DMPI_DIR=$MPI_DIR ..
make VERBOSE=1

# Run the tests
cd benchmarks
ldd test3D_C2C
mpirun -n 2 test3D_C2C -lib $FFT -backend fftw -size 4 4 4 -pgrid 1 2

