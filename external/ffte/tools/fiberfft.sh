#! /usr/bin/env bash

ver=7.0
bffte=buildffte
bfiber=buildfiber

date

for mod in xl/16.1.1-10 spectrum-mpi/10.4.0.3-20210112 cmake/3.18.4
do
  module load $mod
done

mkdir /tmp/$USER || true
cd /tmp/$USER

if test -f ffte-${ver}.tgz ; then
  echo :FIBER: Already downloaded
else
  echo :FIBER: Downloadig FFTE
  curl -OL http://ffte.jp/ffte-${ver}.tgz
fi

if test -d ffte-${ver} ; then
  echo :FIBER: Already unpacked
else
  echo :FIBER: Unpacking FFTE
  tar -xzf ffte-${ver}.tgz
fi

if test -d fiber ; then
  echo :FIBER: Already cloned
else
  echo :FIBER: Clonining FIBER branch
  git clone https://github.com/luszczek/fiber
fi
# must be on backend/ffte branch
cd fiber ; git pull --rebase ; git checkout backend/ffte ; cd ..

#
# Setup and build FFTE
#
cp fiber/external/ffte/tools/CMakeLists.txt ffte-${ver} # FFTE doesnt use CMake
cp fiber/external/ffte/tools/fftever.f ffte-${ver} # FFTE version
cp fiber/external/ffte/tools/pfftever.f ffte-${ver}/mpi # parallel FFTE version
cp ffte-${ver}/param.h ffte-${ver}/mpi # FFTE parameters for scratch buffers

if test -d ${bffte} ; then
  echo :FIBER: Build dir already created, cleaning
  rm -r ${bffte}/*
else
  mkdir ${bffte}
fi
echo :FIBER: Building FFTE
cd ${bffte} ; cmake -D CMAKE_INSTALL_PREFIX=/tmp/$USER/usr ../ffte-${ver} ; make -j ; make install ; cd ..

#
# Setup and build FFTE
#
if test -d ${bfiber} ; then
  echo :FIBER: Build dir already created, cleaning
  rm -r ${bfiber}/*
else
  mkdir ${bfiber}
fi
echo :FIBER: Building FIBER
cd ${bfiber} ; cmake -D CMAKE_INSTALL_PREFIX=/tmp/$USER/usr -D FIBER_FFT_LIB_DIRS=/tmp/$USER/usr/lib -D FIBER_ENABLE_FFTE=ON ../fiber ; make -j ; make install ; cd ..

echo :FIBER: Benchmarks installed in /tmp/$USER/${bfiber}/benchmarks

date

exit 1
