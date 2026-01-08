# ===========================================================================
# User-configurable options
# ARCH options = ser | mpi | ampi
# requires gnu compiler collection and openmpi/mpich for parallel version
# ===========================================================================
DEBUG         ?= false
TEST_CASE     ?= climate
PARAM         ?= param_J5
ARCH          ?= mpi
MPIF90        ?= mpif90
PREFIX        ?= .
SRC_DIR       ?= src
BIN_DIR       ?= bin
BUILD_DIR     ?= build
MODDIR        ?= $(BUILD_DIR)/mod

# Libraries
LAPACK        ?= -llapack
NETCDF        ?= -lnetcdff

# =========================
# Derived toggles
# =========================
TOPO    := false
PHYSICS := false
ifeq ($(TEST_CASE),make_NCAR_topo)
TOPO := true
endif
ifeq ($(TEST_CASE),$(filter $(TEST_CASE),climate spherical_harmonics))
PHYSICS := true
endif

# =========================
# Toolchain selection
# =========================
ifeq ($(ARCH),ser)
FC   := gfortran
PROC := ser
else ifeq ($(ARCH),mpi)
FC   := $(MPIF90)
PROC := mpi
else ifeq ($(ARCH),ampi)
CHARM_DIR   ?= $(HOME)/charm
CHARM_BUILD ?= multicore-linux-x86_64
FC          := $(CHARM_DIR)/$(CHARM_BUILD)/bin/mpif90.ampi
PROC        := mpi
else
$(error Unknown ARCH='$(ARCH)'; use ser|mpi|ampi)
endif

# =========================
# Require GNU compiler collection
# =========================
ifeq ($(shell command -v gfortran >/dev/null 2>&1 || echo no),no)
  $(error GNU Fortran compiler (gfortran) not found. Please load/install GCC.)
endif
$(info Using GNU Fortran compiler: $(shell gfortran --version | head -n 1))

LD := $(FC)

# =========================
# Flags (clean separation)
# =========================
CPPFLAGS += -cpp
ifeq ($(ARCH),ser)
# nothing
else
CPPFLAGS += -DMPI
endif
ifeq ($(ARCH),ampi)
CPPFLAGS += -DAMPI
endif
ifeq ($(PHYSICS),true)
CPPFLAGS += -DPHYSICS
endif

# Identify shell
UNAME_S := $(shell uname -s)

# Debug vs release
ifeq ($(DEBUG),true)
ifeq ($(UNAME_S),Darwin)
$(error DEBUG=true is not supported on macOS (ASan/UBSan + mpif90 + -flat_namespace). Please use DEBUG=true on Linux/HPC)
endif

# Run with export ASAN_OPTIONS=abort_on_error=1:detect_leaks=0 when using mpi and > 1 tasks
SANFLAGS := -fsanitize=address,undefined -fno-omit-frame-pointer
FFLAGS   += -O0 -g2 -Wall -Wextra -Wno-unused-dummy-argument -Wno-trampolines -fcheck=all $(SANFLAGS)
LDFLAGS  += $(SANFLAGS) -Wl,-undefined,error

else

FFLAGS   += -O2 -g -Wall -Wextra -Wno-unused-dummy-argument -Wno-trampolines \
-ffpe-trap=invalid,zero,overflow -fbacktrace -ftrapping-math -march=native -funroll-loops

endif


ifeq ($(UNAME_S),Linux)
  # Force non-executable stack (ELF/GNU ld); fixes ".note.GNU-stack is executable" warnings
  LDFLAGS += -Wl,-z,noexecstack
endif

# Modules/includes
FFLAGS   += -J$(MODDIR) -I$(MODDIR)

# Link libraries
LDLIBS   += $(LAPACK)

# =========================
# Machine / platform tweaks
# =========================
ifeq ($(UNAME_S),Darwin)
MACHINE := mac
# Homebrew lapack path (adjust if you want)
LDLIBS  += -L/opt/homebrew/opt/lapack/lib

ifeq ($(TOPO),true)
NETCDF_DIR ?= /opt/homebrew/Cellar/netcdf-fortran/4.6.2
CPPFLAGS   += -I$(NETCDF_DIR)/include
LDFLAGS    += -L$(NETCDF_DIR)/lib
endif
else
MACHINE := $(shell scontrol show config 2>/dev/null | \
awk -F= '/^ClusterName/{sub(/^[[:space:]]+/,"",$$2); sub(/[[:space:]]+$$/,"",$$2); print tolower($$2); exit}')

$(info Machine: $(MACHINE))

# Compute Canada: prefer flexiblas if you load it
ifeq ($(MACHINE),$(filter $(MACHINE),nibi trillium narval))
LAPACK  := -lflexiblas
LDLIBS  += $(LAPACK)
endif

# bbcluster2
ifeq ($(MACHINE),bbcluster2)
ifdef NETLIB_LAPACK_ROOT
CPPFLAGS += -I$(NETLIB_LAPACK_ROOT)/include
LDFLAGS  += -L$(NETLIB_LAPACK_ROOT)/lib64
endif
ifeq ($(TOPO),true)
ifdef NETCDF_FORTRAN_ROOT
CPPFLAGS += -I$(NETCDF_FORTRAN_ROOT)/include
LDFLAGS  += -L$(NETCDF_FORTRAN_ROOT)/lib -Wl,-rpath,$(NETCDF_FORTRAN_ROOT)/lib
endif
endif
endif
endif

ifeq ($(TOPO),true)
LDLIBS += $(NETCDF)
endif

# ==================================================
# Optional SHTOOLS / FFTW for spherical_harmonics
# ==================================================
ifeq ($(TEST_CASE),spherical_harmonics)
ifeq ($(MACHINE),mac)
CPPFLAGS += -I/opt/homebrew/include
LDFLAGS  += -L/opt/homebrew/lib
else ifeq ($(MACHINE),bbcluster2)
CPPFLAGS += -I$(SHTOOLS_ROOT)/include
LDFLAGS  += -L$(SHTOOLS_ROOT)/lib -Wl,-rpath,$(SHTOOLS_ROOT)/lib \
            -L$(FFTW_ROOT)/lib    -Wl,-rpath,$(FFTW_ROOT)/lib
else
CPPFLAGS += -I/usr/local/include
LDFLAGS  += -L/usr/local/lib
endif
LDLIBS   += -lSHTOOLS -lfftw3 -lm
FFLAGS   += -fPIC
endif

# ==================================================
# Test case source links
# (done as real build steps, not parse-time)
# ==================================================
TESTMOD_SRC := $(SRC_DIR)/test_case_module.f90
TEST_SRC    := $(SRC_DIR)/test.f90

.DEFAULT_GOAL := all

all: $(EXE)

$(TESTMOD_SRC):
	ln -nsf ../test/$(TEST_CASE)/test_case_module.f90 $@

$(TEST_SRC):
	ln -nsf ../test/$(TEST_CASE)/$(TEST_CASE).f90 $@

# =========================
# Sources / objects
# =========================
vpath %.f90 $(SRC_DIR)

SRC = kind.f90 $(PARAM).f90 shared.f90 coord_arithmetic.f90 calendar.f90 geom.f90 patch.f90 dyn_array.f90 \
      base_$(PROC).f90 spline.f90 domain.f90 domain_ops.f90 init.f90 comm.f90 comm_$(PROC).f90 utils.f90 \
      projection.f90 equation_of_state.f90 wavelet.f90 lnorms.f90 mask.f90 refine_patch.f90 coarse_grid.f90 ops.f90 \
      multi_level.f90 adapt.f90 lin_solve.f90 barotropic_2d.f90 time_integr.f90 io.f90 vert_diffusion.f90 io_vtk.f90 \
      remap.f90 std_atm_profile.f90 sso.f90

ifeq ($(TOPO),true)
SRC += topo_grid_descriptor.f90
endif

ifeq ($(PHYSICS),true)
SIMPLEPHYSMODPATH := $(SRC_DIR)/physics/simple_physics/phyparam/include
PHYSLIB_PATH      := $(SRC_DIR)/physics/simple_physics/phyparam/driver
CPPFLAGS          += -I$(SIMPLEPHYSMODPATH)
LDFLAGS           += -L$(PHYSLIB_PATH)
LDLIBS            += -lphyparam
-include $(SRC_DIR)/physics/Makefile.inc
endif

# Ensure the link-generated test sources are built before compiling objects that USE them
SRC += main.f90 test_case_module.f90 test.f90

OBJ := $(addprefix $(BUILD_DIR)/,$(SRC:.f90=.o))

# =========================
# Targets
# =========================
EXE := $(PREFIX)/$(BIN_DIR)/$(TEST_CASE)

.PHONY: all clean phys_package topography dirs
all: $(EXE)

dirs: $(BUILD_DIR) $(BIN_DIR) $(MODDIR)

$(BUILD_DIR):
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

$(MODDIR): | $(BUILD_DIR)
	mkdir -p $@

# Link
$(EXE): $(OBJ) | dirs
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LDLIBS)

# Compile (ensure mod dir exists; ensure test symlinks exist)
$(BUILD_DIR)/%.o: %.f90 | dirs
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $< -o $@

# Force creation of the symlinked test sources before they compile
$(BUILD_DIR)/test_case_module.o: $(TESTMOD_SRC)
$(BUILD_DIR)/test.o:            $(TEST_SRC)

# Optional helper targets
phys_package:
	$(MAKE) -C $(SRC_DIR)/physics/simple_physics/phyparam F90=$(MPIF90)

topography:
	$(MAKE) -C topo
	$(MAKE) -C topo clean

clean:
	rm -rf $(BUILD_DIR) $(SRC_DIR)/test_case_module.f90 $(SRC_DIR)/test.f90
ifeq ($(PHYSICS),true)
	$(MAKE) -C $(SRC_DIR)/physics/simple_physics/phyparam clean
endif

