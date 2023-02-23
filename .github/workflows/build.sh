#!/bin/bash --login

MPI=$1
FFT=$2
COMPILER=$3

source /etc/profile
set +x
set -e
trap 'echo "# $BASH_COMMAND"' DEBUG
shopt -s expand_aliases

export HOME=`pwd`
source ../spack/share/spack/setup-env.sh
spack env activate --temp
spack add cmake fftw $MPI $FFT %$COMPILER
spack add $COMPILER
spack install
spack load

# Build the project
mkdir -p build && cd build
FFT_DIR=`spack location -i $FFT %$COMPILER`
FFTW_DIR=`spack location -i fftw %$COMPILER`
MPI_DIR=`spack location -i $MPI %$COMPILER`
export CPATH=$FFT_DIR/include:$FFTW_DIR/include:$MPI_DIR/include
echo CPATH=$FFT_DIR/include:$FFTW_DIR/include:$MPI_DIR/include
export LIBRARY_PATH=$FFT_DIR/lib:$FFTW_DIR/lib:$MPI_DIR/lib
echo LIBRARY_PATH=$FFT_DIR/lib:$FFTW_DIR/lib:$MPI_DIR/lib
export LD_LIBRARY_PATH=$LIBRARY_PATH
LIBNAME=${FFT^^}
cmake -DFIBER_ENABLE_$LIBNAME=ON ..
make VERBOSE=1

# Run the tests
cd benchmarks
ldd test3D_C2C
mpirun -n 2 ./test3D_C2C -lib $FFT -backend fftw -size 4 4 4 -pgrid 1 2

