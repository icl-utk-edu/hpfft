cmake_minimum_required(VERSION 3.10)

project("FIBER" VERSION 0.1 LANGUAGES C)

set(CMAKE_C_FLAGS "-lm")

find_package(MPI REQUIRED) # always a dependency

# Add include paths
include_directories(/home/aayalaob/heffte/build/include)
include_directories(/home/aayalaob/Benchmarks_FFT/fftmpi/src)
include_directories(/home/aayalaob/Benchmarks_FFT/accfft/build/include)
include_directories(/home/aayalaob/Benchmarks_FFT/p3dfft.3/include())
include_directories(/home/aayalaob/Benchmarks_FFT/FFTE-C-1.1/include)
include_directories(/home/aayalaob/Benchmarks_FFT/SWFFT/)
include_directories(/home/aayalaob/Benchmarks_FFT/2decomp_fft/src)
include_directories(/home/aayalaob/Benchmarks_FFT/nb3dfft/src)
include_directories(/sw/spack/2021-02-12/opt/spack/linux-scientific7-x86_64/gcc-7.3.0/fftw-3.3.9-xvbumufiowb3nqnifdffjdqga4gqsz2l/include)
include_directories(/sw/openmpi/4.0.0/include)

include_directories(/home/aayalaob/fiber/include)

# Add library paths
link_directories(/home/aayalaob/heffte/build)
link_directories(/home/aayalaob/Benchmarks_FFT/accfft/build)
link_directories(/home/aayalaob/Benchmarks_FFT/fftmpi/src)
link_directories(/home/aayalaob/Benchmarks_FFT/p3dfft.3/build)
link_directories(/home/aayalaob/Benchmarks_FFT/FFTE-C-1.1/src)
link_directories(/home/aayalaob/Benchmarks_FFT/SWFFT/build)
link_directories(/home/aayalaob/Benchmarks_FFT/2decomp_fft/build)
link_directories(/home/aayalaob/Benchmarks_FFT/nb3dfft/src)
link_directories(/sw/spack/2021-02-12/opt/spack/linux-scientific7-x86_64/gcc-7.3.0/fftw-3.3.9-xvbumufiowb3nqnifdffjdqga4gqsz2l/lib)
link_directories(/sw/openmpi/4.0.0/lib)


macro(fiber_add_benchmark fiber_benchmark)
    add_executable(${fiber_benchmark} "${fiber_benchmark}.c")
    target_link_libraries(${fiber_benchmark} heffte)
    target_link_libraries(${fiber_benchmark} fft3dmpi)
    target_link_libraries(${fiber_benchmark} accfft)
    target_link_libraries(${fiber_benchmark} p3dfft.3)
    target_link_libraries(${fiber_benchmark} FFTE-i-mpi-C)
    target_link_libraries(${fiber_benchmark} fftw3)
    target_link_libraries(${fiber_benchmark} mpi)

# target_include_directories(${fiber_benchmark} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include/>)
endmacro()

fiber_add_benchmark(test3D_CPU_R2C)
fiber_add_benchmark(test3D_CPU_C2C)
