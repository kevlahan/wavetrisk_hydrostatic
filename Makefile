# ===========================================================================
# User-configurable options
#
# Build modes:
#   DEBUG=false  → optimized
#   DEBUG=asan   → ASan/UBSan (Linux/HPC only)
#   DEBUG=check  → Fortran runtime checks and floating-point traps
# ===========================================================================
DEBUG         ?= false
TEST_CASE     ?= climate
PARAM         ?= param_J5
MPIF90        ?= mpif90
PREFIX        ?= .
SRC_DIR       ?= src
BIN_DIR       ?= bin
BUILD_DIR     ?= build
MODDIR        ?= $(BUILD_DIR)/mod

# Default LAPACK implementation
LAPACK        ?= -llapack

# =========================
# Derived toggles
# =========================
PHYSICS := false

ifeq ($(TEST_CASE),$(filter $(TEST_CASE),climate spherical_harmonics))
PHYSICS := true
endif

# =========================
# Toolchain selection
# =========================
FC   := $(MPIF90)
PROC := mpi

LD := $(FC)

# =========================
# Require GNU Fortran
# =========================
ifeq ($(shell command -v gfortran >/dev/null 2>&1 || echo no),no)
$(error GNU Fortran compiler gfortran not found)
endif

$(info Using GNU Fortran compiler: $(shell gfortran --version | head -n 1))

# Operating system
UNAME_S := $(shell uname -s)

# =========================
# Preprocessor flags
# =========================
CPPFLAGS += -cpp

ifeq ($(PHYSICS),true)
CPPFLAGS += -DPHYSICS
endif

# =========================
# Build flags
# =========================
ifeq ($(DEBUG),asan)

ifeq ($(UNAME_S),Darwin)
$(error DEBUG=asan is not supported on macOS; use DEBUG=check)
endif

# For MPI runs with more than one task:
# export ASAN_OPTIONS=abort_on_error=1:detect_leaks=0
SANFLAGS := -fsanitize=address,undefined -fno-omit-frame-pointer

FFLAGS += -O0 -g2 \
          -Wall \
          -Wextra \
          -Wno-unused-dummy-argument \
          -Wno-trampolines \
          -fcheck=all \
          $(SANFLAGS)

LDFLAGS += $(SANFLAGS) -Wl,-undefined,error

else ifeq ($(DEBUG),check)

FFLAGS += -O0 -g \
          -Wall \
          -Wextra \
          -Wno-unused-dummy-argument \
          -Wimplicit-interface \
          -Werror \
          -fmax-errors=10 \
          -fcheck=all \
          -finit-real=snan \
          -finit-integer=-999999 \
          -finit-logical=true \
          -ffpe-trap=invalid,zero,overflow \
          -ftrapping-math \
          -fbacktrace

else ifeq ($(DEBUG),false)

FFLAGS += -O2 -g \
              -Wall -Wextra \
              -Wno-unused-dummy-argument \
              -Wno-trampolines \
              -fbacktrace \
              -march=native

else

$(error Unknown DEBUG='$(DEBUG)'; use false, asan, or check)

endif

# Force a non-executable stack on Linux
ifeq ($(UNAME_S),Linux)
LDFLAGS += -Wl,-z,noexecstack
endif

# Fortran module output and search directory
FFLAGS += -J$(MODDIR) -I$(MODDIR)

# =========================
# NetCDF
# =========================
NF_CONFIG ?= nf-config
NC_CONFIG ?= nc-config

ifeq ($(shell command -v $(NF_CONFIG) >/dev/null 2>&1 || echo no),no)
$(error nf-config not found; load or install NetCDF-Fortran)
endif

ifeq ($(shell command -v $(NC_CONFIG) >/dev/null 2>&1 || echo no),no)
$(error nc-config not found; load or install NetCDF-C)
endif

# Compilation and module include flags
FFLAGS += $(shell $(NF_CONFIG) --fflags)

# nf-config --flibs supplies the libraries in the correct order:
#     -lnetcdff -lnetcdf
#
# On Homebrew, it may omit the separate NetCDF-C -L directory.
NETCDF_LIBS     := $(shell $(NF_CONFIG) --flibs)
NETCDF_C_LIBDIR := $(shell $(NC_CONFIG) --libdir)

# Extract any -L directories supplied by nf-config
NETCDF_F_LIBDIRS := \
	$(sort $(patsubst -L%,%,$(filter -L%,$(NETCDF_LIBS))))

# Add the NetCDF-C search directory without adding -lnetcdf again
LDFLAGS += -L$(NETCDF_C_LIBDIR)

# Runtime search paths
LDFLAGS += $(foreach d,$(NETCDF_F_LIBDIRS),-Wl,-rpath,$(d))
LDFLAGS += -Wl,-rpath,$(NETCDF_C_LIBDIR)

# =========================
# Machine/platform tweaks
# =========================
ifeq ($(UNAME_S),Darwin)

MACHINE := mac

# Homebrew LAPACK
LDFLAGS += -L/opt/homebrew/opt/lapack/lib

else

MACHINE := $(shell scontrol show config 2>/dev/null | \
	awk -F= '/^ClusterName/ { \
		sub(/^[[:space:]]+/,"",$$2); \
		sub(/[[:space:]]+$$/,"",$$2); \
		print tolower($$2); \
		exit \
	}')

$(info Machine: $(MACHINE))

# Alliance Canada systems
ifeq ($(MACHINE),$(filter $(MACHINE),nibi trillium narval))
LAPACK := -lflexiblas
endif

# bbcluster2
ifeq ($(MACHINE),bbcluster2)

ifdef NETLIB_LAPACK_ROOT
CPPFLAGS += -I$(NETLIB_LAPACK_ROOT)/include
LDFLAGS  += -L$(NETLIB_LAPACK_ROOT)/lib64
endif

endif
endif

# Libraries added exactly once
LDLIBS += $(LAPACK)
LDLIBS += $(NETCDF_LIBS)

# ==================================================
# Optional SHTOOLS/FFTW for spherical_harmonics
# ==================================================
ifeq ($(TEST_CASE),spherical_harmonics)

ifeq ($(MACHINE),mac)

CPPFLAGS += -I/opt/homebrew/include
LDFLAGS  += -L/opt/homebrew/lib

else ifeq ($(MACHINE),bbcluster2)

CPPFLAGS += -I$(SHTOOLS_ROOT)/include

LDFLAGS += -L$(SHTOOLS_ROOT)/lib \
           -Wl,-rpath,$(SHTOOLS_ROOT)/lib \
           -L$(FFTW_ROOT)/lib \
           -Wl,-rpath,$(FFTW_ROOT)/lib

else

CPPFLAGS += -I/usr/local/include
LDFLAGS  += -L/usr/local/lib

endif

LDLIBS += -lSHTOOLS -lfftw3 -lm
FFLAGS += -fPIC

endif

# ==================================================
# Test-case source links
# ==================================================
TESTMOD_SRC := $(SRC_DIR)/test_case_module.f90
TEST_SRC    := $(SRC_DIR)/test.f90

$(TESTMOD_SRC):
	ln -nsf ../test/$(TEST_CASE)/test_case_module.f90 $@

$(TEST_SRC):
	ln -nsf ../test/$(TEST_CASE)/$(TEST_CASE).f90 $@

# =========================
# Sources and objects
# =========================
vpath %.f90 $(SRC_DIR)

SRC = kind.f90 \
      $(PARAM).f90 \
      shared.f90 \
      coord_arithmetic.f90 \
      calendar.f90 \
      geom.f90 \
      patch.f90 \
      dyn_array.f90 \
      arch.f90 \
      spline.f90 \
      domain.f90 \
      domain_ops.f90 \
      init.f90 \
      comm.f90 \
      comm_mpi.f90 \
      utils.f90 \
      projection.f90 \
      equation_of_state.f90 \
      wavelet.f90 \
      lnorms.f90 \
      mask.f90 \
      refine_patch.f90 \
      coarse_grid.f90 \
      ops.f90 \
      multi_level.f90 \
      adapt.f90 \
      lin_solve.f90 \
      barotropic_2d.f90 \
      time_integr.f90 \
      checkpoint.f90 \
      NCAR_topo.f90 \
      vert_diffusion.f90 \
      io_vtk.f90 \
      remap.f90 \
      std_atm_profile.f90 \
      sso.f90 \
      topo_grid_descriptor.f90

# =========================
# Optional physics package
# =========================
ifeq ($(PHYSICS),true)

SIMPLEPHYSMODPATH := $(SRC_DIR)/physics/simple_physics/phyparam/include
PHYSLIB_PATH      := $(SRC_DIR)/physics/simple_physics/phyparam/driver

CPPFLAGS += -I$(SIMPLEPHYSMODPATH)
LDFLAGS  += -L$(PHYSLIB_PATH)
LDLIBS   += -lphyparam

-include $(SRC_DIR)/physics/Makefile.inc

endif

SRC += main.f90 \
       test_case_module.f90 \
       test.f90

OBJ := $(addprefix $(BUILD_DIR)/,$(SRC:.f90=.o))

# =========================
# Main target
# =========================
EXE := $(PREFIX)/$(BIN_DIR)/$(TEST_CASE)

.DEFAULT_GOAL := all

.PHONY: all clean phys_package topography dirs

all: $(EXE)

# =========================
# Directory targets
# =========================
dirs: $(BUILD_DIR) $(BIN_DIR) $(MODDIR)

$(BUILD_DIR):
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

$(MODDIR): | $(BUILD_DIR)
	mkdir -p $@

# =========================
# Link
# =========================
$(EXE): $(OBJ) | dirs
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LDLIBS)

# =========================
# Compile
# =========================
$(BUILD_DIR)/%.o: %.f90 | dirs
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $< -o $@

# Ensure symlinked test sources exist before compilation
$(BUILD_DIR)/test_case_module.o: $(TESTMOD_SRC)
$(BUILD_DIR)/test.o: $(TEST_SRC)

# =========================
# Optional helper targets
# =========================
phys_package:
	$(MAKE) -C $(SRC_DIR)/physics/simple_physics/phyparam \
		F90=$(MPIF90)

clean:
	rm -rf $(BUILD_DIR) \
	       $(SRC_DIR)/test_case_module.f90 \
	       $(SRC_DIR)/test.f90

ifeq ($(PHYSICS),true)
	$(MAKE) -C $(SRC_DIR)/physics/simple_physics/phyparam clean
endif
