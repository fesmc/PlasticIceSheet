# Shared build configuration (dependency wiring), included by the generated Makefile after configme's compiler and
# machine fragments.
#
# The library has no dependencies. Optional: Lis (lis=1, the steady SMB-driven solve) and, for the comparison driver
# (make compare) only, fesm-utils (ncio) and netCDF. Both come from fesm-utils, the checkout or link `fesm-utils` here
# (configme links it; in a host such as kryos, the host's fesm-utils), the serial or OpenMP variant as openmp=.

# Floating point: the library tests NaN (RegularizedCoulomb at its cap) and must keep IEEE semantics, while the
# machine fragments optimise with -Ofast. The value-safe flag of the compiler is appended (the last one wins).
FPSAFE_ifx      = -fp-model precise
FPSAFE_ifort    = -fp-model precise
FPSAFE_gfortran = -fno-fast-math
DFLAGS_NODEBUG += $(FPSAFE_$(COMPILER))

# OpenMP build (make openmp=1). The library has no OpenMP regions except through Lis; this matters to match a host
# that builds everything with the flag.
FESMUTILSROOT ?= fesm-utils
ifeq ($(openmp), 1)
    FFLAGS := $(FFLAGS) $(FFLAGS_OPENMP)
    FESMUTILS_VARIANT = omp
else
    FESMUTILS_VARIANT = serial
endif

# Position-independent code (make pic=1), to link into a shared library.
ifeq ($(pic), 1)
    FFLAGS := $(FFLAGS) -fPIC
endif

# Lis (make lis=1): BiCGStab + ILU(0) of the steady SMB-driven solver through Lis instead of the built-in solver;
# fesm-utils' Lis of the same variant (make -C fesm-utils lis), or LISROOT on the command line.
DEFS_LIS =
INC_LIS  =
LIB_LIS  =
ifeq ($(lis), 1)
    LISROOT ?= $(FESMUTILSROOT)/lis/lis-$(FESMUTILS_VARIANT)
    DEFS_LIS = -DPLASTIC_LIS
    INC_LIS  = -I$(LISROOT)/include
    LIB_LIS  = -L$(LISROOT)/lib -llis
endif

# make compare: fesm-utils' ncio and netCDF (LIB_NC from configme)
INC_FESMUTILS = -I$(FESMUTILSROOT)/include-$(FESMUTILS_VARIANT)
LIB_FESMUTILS = -L$(FESMUTILSROOT)/include-$(FESMUTILS_VARIANT) -lfesmutils

# Extra link flags (a machine fragment may set them; macOS ld takes none)
LFLAGS_EXTRA ?=
LFLAGS = $(LFLAGS_EXTRA)
