cmake_minimum_required(VERSION 3.10)

# Most project code is in C. SWFFT needs C++. 2DECOMP&FFT and FFTE need Fortran
project("HPFFT" VERSION 0.1 LANGUAGES C CXX Fortran)

#################
# Default Options
#################
option(HPFFT_ENABLE_HEFFTE        "Enable the HEFFTE backend"       OFF)
option(HPFFT_ENABLE_FFTMPI        "Enable the FFTMPI backend"       OFF)
option(HPFFT_ENABLE_ACCFFT        "Enable the ACCFFT backend"       OFF)
option(HPFFT_ENABLE_P3DFFT        "Enable the P3DFFT backend"       OFF)
option(HPFFT_ENABLE_FFTE          "Enable the FFTE backend"         OFF)
option(HPFFT_ENABLE_SWFFT         "Enable the SWFFT backend"        OFF)
option(HPFFT_ENABLE_2DECOMP       "Enable the 2DECOMP backend"      OFF)
option(HPFFT_ENABLE_NB3DFFT       "Enable the NB3DFFT backend"      OFF)
option(HPFFT_ENABLE_FFTW          "Enable the FFTW backend"         OFF)
option(HPFFT_ENABLE_FFTWPP        "Enable the FFTWPP backend"       OFF)
option(HPFFT_ENABLE_FFTADVMPI     "Enable the FFTADVAMPI backend"   OFF)

option(BUILD_SHARED_LIBS      "Builds shared libraries using CMake conventions"  ON)

set(HPFFT_FFT_INCLUDE_DIRS "" CACHE FILEPATH "Path to FFT library headers")
set(HPFFT_FFT_LIB_DIRS "" CACHE FILEPATH "Path to FFT library files")

set(HPFFT_SWFFT_REPO "$ENV{HOME}/SWFFT" CACHE FILEPATH "Path to SWFFT Git repository")

# Set default install path to build folder
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "default install path" FORCE )
endif()

set(CMAKE_C_FLAGS "-lm")
if(NOT CMAKE_CXX_CREATE_SHARED_LIBRARY)
  set(CMAKE_CXX_CREATE_SHARED_LIBRARY
      "<CMAKE_CXX_COMPILER> <CMAKE_SHARED_LIBRARY_CXX_FLAGS> <LANGUAGE_COMPILE_FLAGS> <LINK_FLAGS> <CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS> <SONAME_FLAG><TARGET_SONAME> -o <TARGET> <OBJECTS> <LINK_LIBRARIES>")
endif()


# HPFFT uses C99 features (in-loop definition of for loop variables)
if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
        set (CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}")
    endif ( CMAKE_C_COMPILER_ID )
else ()
    set ( CMAKE_C_STANDARD 99 )
endif ()

# To add backend flags and paths use CMake flags such as:
# CMAKE_EXE_LINKER_FLAGS

#####################
# Library definitions
#####################
if (HPFFT_ENABLE_HEFFTE)
    add_definitions(-DHPFFT_ENABLE_HEFFTE)
endif()

if (HPFFT_ENABLE_FFTMPI)
    add_definitions(-DHPFFT_ENABLE_FFTMPI)
endif()

if (HPFFT_ENABLE_ACCFFT)
    add_definitions(-DHPFFT_ENABLE_ACCFFT)
endif()

if (HPFFT_ENABLE_P3DFFT)
    add_definitions(-DHPFFT_ENABLE_P3DFFT)
endif()

if (HPFFT_ENABLE_FFTE)
    add_definitions(-DHPFFT_ENABLE_FFTE)
    include(FortranCInterface)
    FortranCInterface_HEADER(${CMAKE_CURRENT_LIST_DIR}/include/hpfft_fortran.h MACRO_NAMESPACE HPFFT_FORTRAN_ SYMBOL_NAMESPACE hpfft_ SYMBOLS pzfft3d pdzfft3d)
endif()

if (HPFFT_ENABLE_SWFFT)
    add_definitions(-DHPFFT_ENABLE_SWFFT)
    add_library(hpfft_swfft benchmarks/hpfft_swfft.cpp ${HPFFT_SWFFT_REPO}/DistributionC.cpp ${HPFFT_SWFFT_REPO}/DfftC.cpp ${HPFFT_SWFFT_REPO}/distribution.c)
endif()

if (HPFFT_ENABLE_2DECOMP)
    add_definitions(-DHPFFT_ENABLE_2DECOMP)
endif()

if (HPFFT_ENABLE_NB3DFFT)
    add_definitions(-DHPFFT_ENABLE_NB3DFFT)
endif()

if (HPFFT_ENABLE_FFTW)
    add_definitions(-DHPFFT_ENABLE_FFTW)
endif()

if (HPFFT_ENABLE_FFTWPP)
    add_definitions(-DHPFFT_ENABLE_FFTWPP)
endif()

if (HPFFT_ENABLE_FFTADVMPI)
    add_definitions(-DHPFFT_ENABLE_FFTADVMPI)
endif()

find_package(MPI REQUIRED)

# Location of the FFT library
include_directories(${HPFFT_FFT_INCLUDE_DIRS})
link_directories(${HPFFT_FFT_LIB_DIRS})

# Location of harness headers
include_directories(include)

# Location of MPI
include_directories(${MPI_C_INCLUDE_DIRS})
link_directories(${MPI_C_LIBRARIES})

#######################
# libhpfft source files
#######################
set(hpfft_common_sources
include/hpfft_backend_accfft.h   
include/hpfft_backend_fftw.h     
include/hpfft_backend_fftwpp.h
include/hpfft_backend_swfft.h
include/hpfft_backend_decomp2d.h 
include/hpfft_backend_heffte.h   
include/hpfft_backends.h
include/hpfft_backend_ffte.h     
include/hpfft_backend_nb3dfft.h  
include/hpfft_utils.h
include/hpfft_backend_fftmpi.h   
include/hpfft_backend_p3dfft.h
include/hpfft_backend_fftadvmpi.h
benchmarks/test3D_CPU_R2C.c
)

# Build libhpfft, with or without GPU support
add_library(hpfft ${hpfft_common_sources})

if (HPFFT_ENABLE_CUDA)
    find_package(CUDA REQUIRED)
    list(APPEND CUDA_NVCC_FLAGS "-std=c++11")

    cuda_add_cufft_to_target(hpfft)
    target_include_directories(hpfft PUBLIC $<INSTALL_INTERFACE:${CUDA_INCLUDE_DIRS}>)
endif()

#########
# install
#########
install(TARGETS hpfft EXPORT hpfft_Targets DESTINATION lib)
install(EXPORT hpfft_Targets FILE hpfftTargets.cmake NAMESPACE hpfft:: DESTINATION lib/cmake/hpfft)

install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h"
                                                              PATTERN "**~" EXCLUDE
                                                              PATTERN "*.cu*" EXCLUDE
                                                              )

##############
# Benchmarks #
##############
add_subdirectory(benchmarks)

###############
# Documentation
###############
if (HPFFT_ENABLE_DOXYGEN)
    add_subdirectory(doxygen)
endif()

# Print summary of the CMake options
if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/print_summary.cmake)
endif()
