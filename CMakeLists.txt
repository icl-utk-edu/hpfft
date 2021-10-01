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

set(FIBER_FFT_INCLUDE_DIRS "" CACHE FILEPATH "Path to FFT library headers")
set(FIBER_FFT_LIB_DIRS "" CACHE FILEPATH "Path to FFT library files")

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

if (FIBER_ENABLE_FFTW)
add_definitions(-DFIBER_ENABLE_FFTW)
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

find_package(MPI REQUIRED)

# Location of the FFT library
include_directories(${FIBER_FFT_INCLUDE_DIRS})
link_directories(${FIBER_FFT_LIB_DIRS})

# Location of FFTW
include_directories(${FFTW_DIR}/include)
link_directories(${FFTW_DIR}/lib)

# Location of harness headers
include_directories(include)

# Location of MPI
include_directories(${MPI_DIR}/include)
link_directories(${MPI_DIR}/lib)

if (FIBER_ENABLE_CUDA)
    find_package(CUDA REQUIRED)
    list(APPEND CUDA_NVCC_FLAGS "-std=c++11")
    cuda_add_cufft_to_target(fiber)
    target_include_directories(fiber PUBLIC $<INSTALL_INTERFACE:${CUDA_INCLUDE_DIRS}>)
endif()

macro(fiber_add_benchmark fiber_benchmark)
    
    add_executable(${fiber_benchmark} "${fiber_benchmark}.c")
    target_link_libraries(${fiber_benchmark} mpi)
    target_include_directories(${fiber_benchmark} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include/>)

    if (FIBER_ENABLE_HEFFTE)
        target_link_libraries(${fiber_benchmark} heffte)
    endif()
    if (FIBER_ENABLE_FFTMPI)
    target_link_libraries(${fiber_benchmark} fft3dmpi)
    endif()
    if (FIBER_ENABLE_ACCFFT)
    target_link_libraries(${fiber_benchmark} accfft)
    endif()
    if (FIBER_ENABLE_P3DFFT)
    target_link_libraries(${fiber_benchmark} p3dfft)
    endif()
    if (FIBER_ENABLE_FFTE)
    target_link_libraries(${fiber_benchmark} ffte)
    endif()
    if (FIBER_ENABLE_SWFFT)
    target_link_libraries(${fiber_benchmark} swfft)
    endif()
    if (FIBER_ENABLE_2DECOMP)
    target_link_libraries(${fiber_benchmark} 2decomp_fft)
    endif()
    if (FIBER_ENABLE_NB3DFFT)
    target_link_libraries(${fiber_benchmark} nb3dfft)
    endif()
    if (FIBER_ENABLE_FFTW)
    target_link_libraries(${fiber_benchmark} fftw3)
    endif()

    endmacro()

fiber_add_benchmark(test3D_CPU_R2C)
fiber_add_benchmark(test3D_CPU_C2C)

if (FIBER_ENABLE_CUDA)
fiber_add_benchmark(test3D_GPU_R2C)
endif()

# Print summary of the CMake options
if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/print_summary.cmake)
endif()
