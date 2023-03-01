include make.inc

VERBOSE=0
CP_FLAGS="-ir"
DEFAULT_PREFIX = ../
PREFIX ?= $(DEFAULT_PREFIX)

DECOMP_INC ?= $(DECOMP_ROOT)/include
DECOMP_LIB ?= $(DECOMP_ROOT)/lib

ifeq ($(VERBOSE), 1)
CP_FLAGS += "-v"
endif

.PHONY: all copy help install lib

all: help

lib:
	make -C interfaces/C/ DECOMP_ROOT=$(DECOMP_ROOT) all

clean:
	make -C interfaces/C/ DECOMP_ROOT=$(DECOMP_ROOT) clean

install:
	make -C interfaces/C/ DECOMP_ROOT=$(DECOMP_ROOT) DECOMP_LIB=$(DECOMP_LIB) install
	cp $(CP_FLAGS) include/* $(DECOMP_INC)

copy:
	@if [ -d $(PREFIX) ]; then \
		echo "Installation in $(PREFIX)"; \
		cp $(CP_FLAGS) examples $(PREFIX); \
		cp $(CP_FLAGS) interfaces $(PREFIX); \
		cp $(CP_FLAGS) include $(PREFIX); \
	else \
		echo "Cannot install in $(PREFIX): does not exist"; \
	fi

help:
	@echo -e 	"[HELP]\n"\
						"To create the interface library, there are two ways:\n"\
						"1/ Classical way: copy the content of this repo into the folder of 2decomp&FFT\n"\
						"For that, we copy the interface in $(DEFAULT_PREFIX) [by default] by doing:\n"\
						"\t make copy\n"\
						"Or we provide the path like this:\n"\
						"\t make PREFIX=<2decomp_src_path> copy\n"\
						"2/ Advanced way: create a make.inc file that sets the following variables:\n"
	@cat make.inc.example
	@echo -e  "\nAfter the creation of the file, we can:\n"\
						"\t make lib\t\t to create the lib in the build folder\n" \
						"\t make install\t\t to copy the header and the created lib in the include and lib folder of the library\n"
	@echo -e  "NOTE: The file make.inc must be located in the same repository as this Makefile"
