# We assume that in the worse case, the interface is located at the root of
# the library
DECOMP_ROOT ?= ../..
BUILD ?= .
DECOMP_LIB ?= $(DECOMP_ROOT)/lib

include $(DECOMP_ROOT)/src/Makefile.inc

LIBNAME=lib2decomp_fft_iface.a
IFACE_BUILD=$(BUILD)/build
IFACE_LIB=$(BUILD)/lib

SRC=decomp_2d_ciface.F90

OBJ=$(addprefix $(IFACE_BUILD)/, $(SRC:.F90=.o))

.PHONY: all lib

all: $(OBJ) lib

lib: $(OBJ)
	@if [ ! -d $(IFACE_LIB) ]; then mkdir -pv $(IFACE_LIB); fi
	ar vr $(IFACE_LIB)/$(LIBNAME) $^
	ranlib $(IFACE_LIB)/$(LIBNAME)

$(IFACE_BUILD)/%.o: %.F90
	@echo "Wants to compile $<"
	@if [ ! -d $(IFACE_BUILD) ]; then mkdir -pv $(IFACE_BUILD); fi
	$(F90) $(OPTIONS) $(F90FLAGS) -o $@ -I $(DECOMP_ROOT)/include -c $<

clean:
	-${RM} $(OBJ)
	-${RM} $(IFACE_LIB)/$(LIBNAME)

install:
	cp -v $(IFACE_LIB)/$(LIBNAME) $(DECOMP_LIB)
