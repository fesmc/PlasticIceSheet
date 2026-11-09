# Shared build configuration, loaded after the config/<machine>_<compiler> fragment.
#
# The library has no dependencies. Only the comparison driver (make compare) links
# fesm-utils (ncio) and netCDF; the fragment or the command line sets FESMUTILSROOT
# and LIB_NC for it.

# OpenMP build (make openmp=1). The library has no OpenMP regions; this only matters
# to match a host that builds everything with the flag.
ifeq ($(openmp), 1)
    FFLAGS := $(FFLAGS) $(FFLAGS_OPENMP)
endif

# Position-independent code (make pic=1), to link into a shared library.
ifeq ($(pic), 1)
    FFLAGS := $(FFLAGS) -fPIC
endif

# Lis (make lis=1): BiCGStab + ILU(0) of the steady SMB-driven solver through Lis
# instead of the built-in solver. LISROOT (fragment or command line) is the Lis
# install matching openmp= (fesm-utils' lis-omp / lis-serial, ifx-built).
DEFS_LIS =
INC_LIS  =
LIB_LIS  =
ifeq ($(lis), 1)
    ifeq ($(openmp), 1)
        LISROOT ?= $(LISROOT_OMP)
    else
        LISROOT ?= $(LISROOT_SERIAL)
    endif
    DEFS_LIS = -DPLASTIC_LIS
    INC_LIS  = -I$(LISROOT)/include
    LIB_LIS  = -L$(LISROOT)/lib -llis
endif

LFLAGS_EXTRA ?=
LFLAGS = $(LFLAGS_EXTRA)

FESMUTILSROOT ?= ../fesm-utils
INC_FESMUTILS = -I$(FESMUTILSROOT)/include-serial
LIB_FESMUTILS = -L$(FESMUTILSROOT)/include-serial -lfesmutils
