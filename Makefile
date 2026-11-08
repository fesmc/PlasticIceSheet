.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# Usage:
#   make static [config=dkrz_levante_ifx] [debug=1] [openmp=1] [pic=1]
#   make check  [config=...]       build and run the test programs
#   make clean                     (between compilers: .mod files are compiler-specific)
#
# The compiler/machine fragment is config/$(config); config/common.mk follows it.

config ?= dkrz_levante_ifx
debug  ?= 0
openmp ?= 0
pic    ?= 0

srcdir  = src
objdir  = libplastic/include
bindir  = libplastic/bin
testdir = tests

# `make` alone prints the targets (the first rule would otherwise be an object rule).
.DEFAULT_GOAL := usage

include config/$(config)
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
               $(objdir)/plastic_solver.o \
               $(objdir)/plastic_ice.o

$(objdir)/plastic_kinds.o: $(srcdir)/plastic_kinds.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/plastic_eikonal.o: $(srcdir)/plastic_eikonal.f90 $(objdir)/plastic_kinds.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/plastic_solver.o: $(srcdir)/plastic_solver.f90 $(objdir)/plastic_eikonal.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/plastic_ice.o: $(srcdir)/plastic_ice.f90 $(objdir)/plastic_solver.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

###############################################
## Targets
###############################################

.PHONY: usage all static clean

usage:
	@echo "make static [config=<machine>_<compiler>]   -> $(objdir)/libplastic.a + .mod files"
	@echo "make check  [config=...]                    -> build and run the tests"
	@echo "make clean"
	@echo "configs: $(notdir $(filter-out config/common.mk,$(wildcard config/*)))"

all: static

static: $(objdir)/libplastic.a

$(objdir)/libplastic.a: $(plastic_objs)
	ar rc $@ $(plastic_objs)
	ranlib $@
	@echo " "
	@echo "    $(objdir)/libplastic.a is ready (FC=$(FC))."
	@echo "    link with: -I$(CURDIR)/$(objdir) -L$(CURDIR)/$(objdir) -lplastic"
	@echo " "

clean:
	rm -f $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(bindir)/*.x
