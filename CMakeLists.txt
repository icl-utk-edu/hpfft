cmake_minimum_required(VERSION 3.10)

# Most project code is in C. SWFFT needs C++. 2DECOMP&FFT and FFTE need Fortran
project("FIBER" VERSION 0.1 LANGUAGES C CXX Fortran)

#################
# Default Options
#################
option(FIBER_ENABLE_HEFFTE        "Enable the HEFFTE backend"       OFF)
option(FIBER_ENABLE_FFTMPI        "Enable the FFTMPI backend"       OFF)
option(FIBER_ENABLE_ACCFFT        "Enable the ACCFFT backend"       OFF)
option(FIBER_ENABLE_P3DFFT        "Enable the P3DFFT backend"       OFF)
option(FIBER_ENABLE_FFTE          "Enable the FFTE backend"         OFF)
option(FIBER_ENABLE_SWFFT         "Enable the SWFFT backend"        OFF)
option(FIBER_ENABLE_2DECOMP       "Enable the 2DECOMP backend"      OFF)
option(FIBER_ENABLE_NB3DFFT       "Enable the NB3DFFT backend"      OFF)
option(FIBER_ENABLE_FFTW          "Enable the FFTW backend"         OFF)
option(FIBER_ENABLE_FFTWPP        "Enable the FFTWPP backend"       OFF)
option(FIBER_ENABLE_FFTADVMPI     "Enable the FFTADVAMPI backend"   OFF)

option(BUILD_SHARED_LIBS      "Builds shared libraries using CMake conventions"  ON)

set(FIBER_FFT_INCLUDE_DIRS "" CACHE FILEPATH "Path to FFT library headers")
set(FIBER_FFT_LIB_DIRS "" CACHE FILEPATH "Path to FFT library files")

set(FIBER_SWFFT_REPO "$ENV{HOME}/SWFFT" CACHE FILEPATH "Path to SWFFT Git repository")

# Set default install path to build folder
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "default install path" FORCE )
endif()

set(CMAKE_C_FLAGS "-lm")
if(NOT CMAKE_CXX_CREATE_SHARED_LIBRARY)
  set(CMAKE_CXX_CREATE_SHARED_LIBRARY
      "<CMAKE_CXX_COMPILER> <CMAKE_SHARED_LIBRARY_CXX_FLAGS> <LANGUAGE_COMPILE_FLAGS> <LINK_FLAGS> <CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS> <SONAME_FLAG><TARGET_SONAME> -o <TARGET> <OBJECTS> <LINK_LIBRARIES>")
endif()


# FIBER uses C99 features (in-loop definition of for loop variables)
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
    include(FortranCInterface)
    FortranCInterface_HEADER(${CMAKE_CURRENT_LIST_DIR}/include/fiber_fortran.h MACRO_NAMESPACE FIBER_FORTRAN_ SYMBOL_NAMESPACE fiber_ SYMBOLS pzfft3d pdzfft3d)
endif()

if (FIBER_ENABLE_SWFFT)
    add_definitions(-DFIBER_ENABLE_SWFFT)
    add_library(fiber_swfft benchmarks/fiber_swfft.cpp ${FIBER_SWFFT_REPO}/DistributionC.cpp ${FIBER_SWFFT_REPO}/DfftC.cpp ${FIBER_SWFFT_REPO}/distribution.c)
endif()

if (FIBER_ENABLE_2DECOMP)
    add_definitions(-DFIBER_ENABLE_2DECOMP)
    add_subdirectory(external/2decompfft)
    include_directories(PUBLIC ${CMAKE_CURRENT_LIST_DIR}/external/2decompfft/include)
    link_directories(${CMAKE_BINARY_DIR}/external/2decompfft/lib)
endif()

if (FIBER_ENABLE_NB3DFFT)
    add_definitions(-DFIBER_ENABLE_NB3DFFT)
endif()

if (FIBER_ENABLE_FFTW)
    add_definitions(-DFIBER_ENABLE_FFTW)
endif()

if (FIBER_ENABLE_FFTWPP)
    add_definitions(-DFIBER_ENABLE_FFTWPP)
endif()

if (FIBER_ENABLE_FFTADVMPI)
    add_definitions(-DFIBER_ENABLE_FFTADVMPI)
endif()

find_package(MPI REQUIRED)

# Location of the FFT library
include_directories(${FIBER_FFT_INCLUDE_DIRS})
link_directories(${FIBER_FFT_LIB_DIRS})

# Location of harness headers
include_directories(include)

# Location of MPI
include_directories(${MPI_C_INCLUDE_DIRS})
link_directories(${MPI_C_LIBRARIES})

#######################
# libfiber source files
#######################
set(fiber_common_sources
include/fiber_backend_accfft.h   
include/fiber_backend_fftw.h     
include/fiber_backend_fftwpp.h
include/fiber_backend_swfft.h
include/fiber_backend_decomp2d.h 
include/fiber_backend_heffte.h   
include/fiber_backends.h
include/fiber_backend_ffte.h     
include/fiber_backend_nb3dfft.h  
include/fiber_utils.h
include/fiber_backend_fftmpi.h   
include/fiber_backend_p3dfft.h
include/fiber_backend_fftadvmpi.h
benchmarks/test3D_CPU_R2C.c
)

# Build libfiber, with or without GPU support
add_library(fiber ${fiber_common_sources})

if (FIBER_ENABLE_CUDA)
    find_package(CUDA REQUIRED)
    list(APPEND CUDA_NVCC_FLAGS "-std=c++11")

    cuda_add_cufft_to_target(fiber)
    target_include_directories(fiber PUBLIC $<INSTALL_INTERFACE:${CUDA_INCLUDE_DIRS}>)
endif()

#########
# install
#########
install(TARGETS fiber EXPORT fiber_Targets DESTINATION lib)
install(EXPORT fiber_Targets FILE fiberTargets.cmake NAMESPACE fiber:: DESTINATION lib/cmake/fiber)

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
if (FIBER_ENABLE_DOXYGEN)
    add_subdirectory(doxygen)
endif()

# Print summary of the CMake options
if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/print_summary.cmake)
endif()
