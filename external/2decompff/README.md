# Purpose
This repo provides a C interface to the 2decomp&FFT library originally written in Fortran.

# Use

This interface should be integrated into the original library.
An extra library will be created called *2decomp_fft_iface*.

## Simple Installation
In order to install the interface, do:
```
make PREFIX=<2decomp_root_path> copy
```

Or, if this repo is placed at the root of the 2decomp&FFT library, just do:
```
make copy
```

Note: you also can create a make.inc that sets some variables more permanently.
If needed, do `make` to get a help.

## Compilation
To create the interface library, there are two ways:
- Classical way: copy the content of this repo into the folder of 2decomp&FFT
For that, we copy the interface by doing:
```make copy```
Or we provide the path like this:
```make PREFIX=<2decomp_src_path> copy```
- Advanced way: create a make.inc file that sets the following variables:
After the creation of the file, we can:
`make lib` to create the lib in the build folder
`make install` to copy the header and the created lib in the include and lib folder of the library

NOTE: The file make.inc must be located in the same repository as the root Makefile
