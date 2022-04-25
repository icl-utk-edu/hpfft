project(install_2decomp_iface)

cmake_minimum_required(VERSION 3.22)

include(ExternalProject)

# FFTW has to be set in the env for 2decomp Makefile.inc to work properly
if (NOT ENV{FFTW_PATH})
  if (NOT DEFINED FFTW_PATH)
    message( FATAL_ERROR 
      "FFTW_PATH cannot be found but is required for 2decomp src/Makefile.inc file."
      " Please provide it through env or by defining it -DFFTW_PATH." )
  endif ()
endif ()
message( STATUS "FFTW path: ${FFTW_PATH}" )

# TODO use findpackage if possible
if (NOT DEFINED ENV{DECOMPFFT_ROOT})
  message( FATAL_ERROR "Cannot find DECOMPFFT_ROOT in the env" )
endif ()

set( DECOMPFFT_ROOT $ENV{DECOMPFFT_ROOT} )
set( DECOMPFFT_IFACE_LIB ${CMAKE_CURRENT_BINARY_DIR} )
set( DECOMPFFT_IFACE_INC ${CMAKE_CURRENT_BINARY_DIR} )
message( STATUS "2DECOMP&FFT root: ${DECOMPFFT_ROOT}" )
message( STATUS "2DECOMP&FFT interface lib: ${DECOMPFFT_IFACE_LIB}" )
message( STATUS "2DECOMP&FFT interface include: ${DECOMPFFT_IFACE_INC}" )

ExternalProject_Add(
    2decomp_c_iface 
    PREFIX ../internal_cmake_dl_files
    SOURCE_DIR 2decomp_fft_build
    GIT_REPOSITORY https://github.com/cayrols/2decompFFT_c_iface.git 
    GIT_TAG master
    CONFIGURE_COMMAND cd <SOURCE_DIR> && cp make.inc.example make.inc
    BUILD_COMMAND make -C <SOURCE_DIR> lib
    INSTALL_COMMAND make -C <SOURCE_DIR> DECOMP_LIB=${DECOMPFFT_IFACE_LIB} DECOMP_INC=${DECOMPFFT_IFACE_INC} install
)
