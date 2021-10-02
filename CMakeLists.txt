cmake_minimum_required(VERSION 3.10)

project("FIBER" VERSION 0.1 LANGUAGES C)

option(FIBER_ENABLE_HEFFTE    "Enable the HEFFTE backend"       OFF)
option(FIBER_ENABLE_FFTMPI    "Enable the FFTMPI backend"       OFF)
option(FIBER_ENABLE_ACCFFT    "Enable the ACCFFT backend"       OFF)
option(FIBER_ENABLE_P3DFFT    "Enable the P3DFFT backend"       OFF)
option(FIBER_ENABLE_FFTE      "Enable the FFTE backend"         OFF)
option(FIBER_ENABLE_SWFFT     "Enable the SWFFT backend"        OFF)
option(FIBER_ENABLE_2DECOMP   "Enable the 2DECOMP backend"      OFF)
option(FIBER_ENABLE_NB3DFFT   "Enable the NB3DFFT backend"      OFF)
option(FIBER_ENABLE_FFTW      "Enable the FFTW backend"         OFF)

option(BUILD_SHARED_LIBS      "Builds shared libraries using CMake conventions"  ON)

set(FIBER_FFT_INCLUDE_DIRS "" CACHE FILEPATH "Path to FFT library headers")
set(FIBER_FFT_LIB_DIRS "" CACHE FILEPATH "Path to FFT library files")

# Set default install path to build folder
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "default install path" FORCE )
endif()

set(CMAKE_C_FLAGS "-lm")

# FIBER uses C99 features (in-loop definition of for loop variables)
if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
        set (CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}")
    endif ( CMAKE_C_COMPILER_ID )
else ()
    set ( CMAKE_C_STANDARD 99 )
endif ()

if (FIBER_ENABLE_HEFFTE)
    add_definitions(-DFIBER_ENABLE_HEFFTE)
endif()

if (FIBER_ENABLE_FFTMPI)
    add_definitions(-DFIBER_ENABLE_FFTMPI)
endif()

if (FIBER_ENABLE_ACCFFT)
    add_definitions(-DFIBER_ENABLE_ACCFFT)
endif()

if (FIBER_ENABLE_P3DFFT)
    add_definitions(-DFIBER_ENABLE_P3DFFT)
endif()

if (FIBER_ENABLE_FFTE)
    add_definitions(-DFIBER_ENABLE_FFTE)
endif()

if (FIBER_ENABLE_SWFFT)
    add_definitions(-DFIBER_ENABLE_SWFFT)
endif()

if (FIBER_ENABLE_2DECOMP)
    add_definitions(-DFIBER_ENABLE_2DECOMP)
endif()

if (FIBER_ENABLE_NB3DFFT)
    add_definitions(-DFIBER_ENABLE_NB3DFFT)
endif()

if (FIBER_ENABLE_FFTW)
    add_definitions(-DFIBER_ENABLE_FFTW)
endif()

find_package(MPI REQUIRED)

# Location of the FFT library
include_directories(${FIBER_FFT_INCLUDE_DIRS})
link_directories(${FIBER_FFT_LIB_DIRS})

# Location of harness headers
include_directories(include)

# Location of MPI
include_directories(${MPI_DIR}/include)
link_directories(${MPI_DIR}/lib)

# libfiber source files
set(fiber_common_sources
include/fiber_backend_accfft.h   
include/fiber_backend_fftw.h     
include/fiber_backend_swfft.h
include/fiber_backend_decomp2d.h 
include/fiber_backend_heffte.h   
include/fiber_backends.h
include/fiber_backend_ffte.h     
include/fiber_backend_nb3dfft.h  
include/fiber_utils.h
include/fiber_backend_fftmpi.h   
include/fiber_backend_p3dfft.h
# benchmarks/Bench_3D_R2C_CPU.c
benchmarks/test3D_CPU_R2C.c
# benchmarks/test3D_CPU_C2C.c
)

# Build libfiber, with or without GPU support
if (FIBER_ENABLE_CUDA)
    find_package(CUDA REQUIRED)
    cuda_add_library(fiber ${fiber_common_sources}
                            benchmarks/test3D_GPU_R2C.c)

    cuda_add_cufft_to_target(fiber)
else()
    add_library(fiber ${fiber_common_sources})
endif()

##############
# Benchmarks #
##############
add_subdirectory(benchmarks)

###############
# Documentation
###############
if (FIBER_ENABLE_DOXYGEN)
    add_subdirectory(doxygen)
endif()

# Print summary of the CMake options
if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/print_summary.cmake)
endif()
