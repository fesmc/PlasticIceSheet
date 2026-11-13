.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# The template of the Makefile: `configme config PlasticIceSheet` writes ./Makefile from it, with the compiler and
# machine fragments of configme in place of the placeholder below (no machine files in this repository).
# Usage (of the generated Makefile):
#   make static [debug=1] [openmp=1] [pic=1] [lis=0]   (lis=1 default: the Lis of the linked fesm-utils; lis=0 the built-in solver)
#   make check                     build and run the test programs
#   make clean                     (between compilers: .mod files are compiler-specific)

debug  ?= 0
openmp ?= 0
pic    ?= 0
lis    ?= 1

srcdir  = src
objdir  = libplastic/include
bindir  = libplastic/bin
testdir = tests

# `make` alone prints the targets (the first rule would otherwise be an object rule).
.DEFAULT_GOAL := usage

## COMPILER CONFIGURATION ##
# Assembled by configme: the compiler fragment (FC, FFLAGS, FFLAGS_OPENMP, CPPFLAGS_PP, DFLAGS_*), the machine
# fragment (host-specific overrides) and netCDF (for make compare only).

<COMPILER_CONFIGURATION>

## DEPENDENCY WIRING ##
# After the fragments, since it references what they set.

include config/common.mk

DFLAGS = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS = $(DFLAGS_DEBUG)
endif
ifeq ($(debug), 2)
	DFLAGS = $(DFLAGS_PROFILE)
endif

###############################################
## Library sources (no external dependencies)
###############################################

plastic_objs = $(objdir)/plastic_kinds.o \
               $(objdir)/plastic_eikonal.o \
               $(objdir)/plastic_linsolve.o \
               $(objdir)/plastic_solver.o \
               $(objdir)/plastic_flow.o \
               $(objdir)/plastic_steady.o \
               $(objdir)/plastic_inversion.o \
               $(objdir)/plastic_reconstruct.o \
               $(objdir)/plastic_ice.o

$(objdir)/plastic_kinds.o: $(srcdir)/plastic_kinds.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/plastic_eikonal.o: $(srcdir)/plastic_eikonal.f90 $(objdir)/plastic_kinds.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/plastic_linsolve.o: $(srcdir)/plastic_linsolve.F90 $(objdir)/plastic_kinds.o
	$(FC) $(DFLAGS) $(FFLAGS) $(CPPFLAGS_PP) $(DEFS_LIS) $(INC_LIS) -c -o $@ $<

$(objdir)/plastic_solver.o: $(srcdir)/plastic_solver.f90 $(objdir)/plastic_eikonal.o $(objdir)/plastic_linsolve.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/plastic_flow.o: $(srcdir)/plastic_flow.f90 $(objdir)/plastic_solver.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/plastic_steady.o: $(srcdir)/plastic_steady.f90 $(objdir)/plastic_flow.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/plastic_inversion.o: $(srcdir)/plastic_inversion.f90 $(objdir)/plastic_steady.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/plastic_reconstruct.o: $(srcdir)/plastic_reconstruct.f90 $(objdir)/plastic_steady.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/plastic_ice.o: $(srcdir)/plastic_ice.f90 $(objdir)/plastic_reconstruct.o $(objdir)/plastic_inversion.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

###############################################
## Targets
###############################################

.PHONY: usage all static test check compare clean

usage:
	@echo "make static [openmp=1] [lis=0] [debug=1] [pic=1]   -> $(objdir)/libplastic.a + .mod files"
	@echo "make check                                        -> build and run the tests"
	@echo "make compare                                      -> comparison driver (fesm-utils ncio)"
	@echo "make clean"
	@echo "configured for $(MACHINE), $(COMPILER) (FC=$(FC))"

all: static

static: $(objdir)/libplastic.a

$(objdir)/libplastic.a: $(plastic_objs)
	ar rc $@ $(plastic_objs)
	ranlib $@
	@echo " "
	@echo "    $(objdir)/libplastic.a is ready (FC=$(FC))."
	@echo "    link with: -I$(CURDIR)/$(objdir) -L$(CURDIR)/$(objdir) -lplastic $(LIB_LIS)"
	@echo " "

# Test programs: no dependencies beyond the library.
check_progs = test_core test_closures test_surface test_steady test_inversion

$(objdir)/plastic_test_utils.o: $(testdir)/test_utils.f90 $(objdir)/libplastic.a
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(bindir)/%.x: $(testdir)/%.f90 $(objdir)/plastic_test_utils.o $(objdir)/libplastic.a
	$(FC) $(DFLAGS) $(FFLAGS) -o $@ $< $(objdir)/plastic_test_utils.o -L$(objdir) -lplastic $(LIB_LIS) $(LFLAGS)

test: $(addprefix $(bindir)/, $(addsuffix .x, $(check_progs)))

# Build and run the tests; fails on the first failing program.
check: test
	@for t in $(check_progs); do \
		$(bindir)/$$t.x || { echo "FAILED: $$t"; exit 1; }; \
	done
	@echo "make check: all tests passed (FC=$(FC))"

# Comparison driver with PlasticIceSheet.jl (compare/run.sh); needs fesm-utils ncio.
compare: $(objdir)/libplastic.a
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/plastic_compare.x compare/plastic_compare.f90 \
		-L$(objdir) -lplastic $(LIB_LIS) $(LIB_FESMUTILS) $(LIB_NC) $(LFLAGS)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/smb_diag.x compare/smb_diag.f90 \
		-L$(objdir) -lplastic $(LIB_LIS) $(LIB_FESMUTILS) $(LIB_NC) $(LFLAGS)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/smb_sia_proto.x compare/smb_sia_proto.f90 \
		-L$(objdir) -lplastic $(LIB_LIS) $(LIB_FESMUTILS) $(LIB_NC) $(LFLAGS)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/bench_threads.x compare/bench_threads.f90 \
		-L$(objdir) -lplastic $(LIB_LIS) $(LIB_FESMUTILS) $(LIB_NC) $(LFLAGS)
	$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/plastic_native.x compare/native/plastic_native.f90 \
		-L$(objdir) -lplastic $(LIB_LIS) $(LFLAGS)

clean:
	rm -f $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.smod $(bindir)/*.x
