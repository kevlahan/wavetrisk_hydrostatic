# ===========================================================================
# User-configurable options
#
# Build modes:
#   DEBUG=false  → optimized
#   DEBUG=asan   → ASan/UBSan (Linux/HPC only)
#   DEBUG=check  → Fortran runtime checks and floating-point traps
# ===========================================================================
DEBUG         ?= false
# Disable contraction for reproducible checked/optimized comparisons.
# FP_CONTRACT=default restores compiler behavior. Use fresh BUILD_DIR and
# BIN_DIR values when changing the arithmetic setting.
FP_CONTRACT   ?= off
# Validated batching and resident inverse workspace (Stage229).
# Use fresh BUILD_DIR and BIN_DIR values when changing these options.
# Stage226 comparison: STAGE202=on STAGE203=off; original path: both off.
STAGE202      ?= on
STAGE203      ?= on
# Bounded interstage pressure-owner integration; not yet the full tendency path.
BLOCK_COLUMNS ?= off
BLOCK_PRIMITIVES ?= off
BLOCK_DIFFUSION ?= off
BLOCK_MASS_OWNER ?= off
BLOCK_VELOCITY_OWNER ?= off
BLOCK_TENDENCY_INITIAL ?= off
BLOCK_RK ?= off
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
FFLAGS += -ffree-line-length-132

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
          -Wline-truncation \
          -Wno-trampolines \
          -fcheck=all \
          $(SANFLAGS)

LDFLAGS += $(SANFLAGS) -Wl,-undefined,error

else ifeq ($(DEBUG),check)

FFLAGS += -O0 -g \
          -Wall \
          -Wextra \
          -Wno-unused-dummy-argument \
          -Wline-truncation \
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
              -Wall \
              -Wextra \
              -Wline-truncation \
              -Wno-unused-dummy-argument \
              -Wno-trampolines \
              -fbacktrace \
              -march=native

else

$(error Unknown DEBUG='$(DEBUG)'; use false, asan, or check)

endif

# Optional arithmetic control for solver code; external physics keeps its flags.
ifneq ($(filter $(BLOCK_COLUMNS),on check),)
ifneq ($(STAGE203),on)
$(error BLOCK_COLUMNS=on requires STAGE203=on)
endif
FFLAGS += -DBLOCK_COLUMNS
endif
ifeq ($(BLOCK_COLUMNS),check)
FFLAGS += -DBLOCK_COLUMNS_ORACLE
endif
ifneq ($(filter $(BLOCK_PRIMITIVES),on check),)
ifeq ($(BLOCK_COLUMNS),off)
$(error BLOCK_PRIMITIVES requires BLOCK_COLUMNS=on or check)
endif
FFLAGS += -DBLOCK_PRIMITIVES
endif
ifeq ($(BLOCK_PRIMITIVES),check)
FFLAGS += -DBLOCK_PRIMITIVES_ORACLE
endif
ifeq ($(filter $(BLOCK_PRIMITIVES),off on check),)
$(error BLOCK_PRIMITIVES must be off, on or check)
endif
ifeq ($(filter $(BLOCK_COLUMNS),off on check),)
$(error BLOCK_COLUMNS must be off, on or check)
endif
ifeq ($(filter $(FP_CONTRACT),default off fast),)
$(error Unknown FP_CONTRACT='$(FP_CONTRACT)'; use default, off, or fast)
endif
ifeq ($(filter $(STAGE203),off fields inverse on),)
$(error Unknown STAGE203='$(STAGE203)'; use off, fields, inverse, or on)
endif
ifneq ($(STAGE203),off)
ifneq ($(STAGE202),on)
$(error STAGE203 requires STAGE202=on)
endif
endif
ifneq ($(filter $(STAGE203),fields on),)
FFLAGS += -DSTAGE203_FIELDS
endif
ifneq ($(filter $(STAGE203),inverse on),)
FFLAGS += -DSTAGE203_INVERSE
endif
ifeq ($(filter $(STAGE202),off records inverse on),)
$(error Unknown STAGE202='$(STAGE202)'; use off, records, inverse, or on)
endif
ifneq ($(filter $(STAGE202),records on),)
FFLAGS += -DSTAGE202_RECORDS
endif
ifneq ($(filter $(STAGE202),inverse on),)
FFLAGS += -DSTAGE202_INVERSE
endif

ifneq ($(FP_CONTRACT),default)
FFLAGS += -ffp-contract=$(FP_CONTRACT)
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
      domain.f90 \
      parallel_block_primitive_walk.f90 \
      parallel_block_resident_operands.f90 \
      parallel_block_primitive_registry.f90 \
      parallel_block_primitives.f90 \
      parallel_block_pentagons.f90 \
      parallel_block_gradients.f90 \
      parallel_block_primitive_kernels.f90 \
      parallel_block_columns.f90 \
      parallel_block_profile.f90 \
      parallel_block_scalar_storage.f90 \
      parallel_block.f90 \
      parallel_block_velocity.f90 \
      parallel_block_mass.f90 \
      parallel_block_program.f90 \
      parallel_block_inverse.f90 \
      parallel_block_mpi.f90 \
      spline.f90 \
      domain_ops.f90 \
      init.f90 \
      comm.f90 \
      comm_mpi.f90 \
      utils.f90 \
      parallel_block_diffusion.f90 \
      integrate.f90 \
      diagnostics.f90 \
      projection.f90 \
      equation_of_state.f90 \
      wavelet.f90 \
      lnorms.f90 \
      mask.f90 \
      refine_patch.f90 \
      coarse_grid.f90 \
      ops.f90 \
      parallel_block_build.f90 \
      parallel_block_definition_exchange.f90 \
      parallel_block_velocity_owner.f90 \
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

# Explicit dependencies for the parallel-block modules.  These are
# required in addition to source-list ordering when make is allowed to
# update more than one target concurrently.
$(BUILD_DIR)/parallel_block.o: \
	$(BUILD_DIR)/parallel_block_profile.o \
	$(BUILD_DIR)/kind.o \
	$(BUILD_DIR)/shared.o \
	$(BUILD_DIR)/patch.o

$(BUILD_DIR)/parallel_block_inverse.o: \
	$(BUILD_DIR)/kind.o \
	$(BUILD_DIR)/shared.o \
	$(BUILD_DIR)/arch.o \
	$(BUILD_DIR)/domain.o \
	$(BUILD_DIR)/parallel_block.o

$(BUILD_DIR)/parallel_block_mpi.o: $(BUILD_DIR)/parallel_block_inverse.o $(BUILD_DIR)/parallel_block_velocity.o $(BUILD_DIR)/parallel_block_mass.o $(BUILD_DIR)/parallel_block_scalar_storage.o
$(BUILD_DIR)/parallel_block_scalar_storage.o: $(BUILD_DIR)/kind.o

$(BUILD_DIR)/parallel_block_velocity.o: $(BUILD_DIR)/kind.o
$(BUILD_DIR)/comm_mpi.o: $(BUILD_DIR)/parallel_block_profile.o
$(BUILD_DIR)/adapt.o $(BUILD_DIR)/time_integr.o: $(BUILD_DIR)/parallel_block_profile.o
$(BUILD_DIR)/parallel_block_mass.o: $(BUILD_DIR)/kind.o

$(BUILD_DIR)/multi_level.o: $(BUILD_DIR)/parallel_block_velocity.o

$(BUILD_DIR)/ops.o: \
	$(BUILD_DIR)/parallel_block_mpi.o \
	$(BUILD_DIR)/kind.o \
	$(BUILD_DIR)/shared.o \
	$(BUILD_DIR)/diagnostics.o \
	$(BUILD_DIR)/domain_ops.o \
	$(BUILD_DIR)/patch.o \
	$(BUILD_DIR)/utils.o \
	$(BUILD_DIR)/init.o \
	$(BUILD_DIR)/domain.o

$(BUILD_DIR)/parallel_block_build.o: \
	$(BUILD_DIR)/kind.o \
	$(BUILD_DIR)/shared.o \
	$(BUILD_DIR)/ops.o \
	$(BUILD_DIR)/patch.o \
	$(BUILD_DIR)/arch.o \
	$(BUILD_DIR)/domain.o \
	$(BUILD_DIR)/parallel_block.o

$(BUILD_DIR)/multi_level.o: \
	$(BUILD_DIR)/parallel_block_mpi.o

$(BUILD_DIR)/adapt.o $(BUILD_DIR)/remap.o: $(BUILD_DIR)/parallel_block_mpi.o

$(BUILD_DIR)/time_integr.o: \
	$(BUILD_DIR)/parallel_block_mpi.o \
	$(BUILD_DIR)/adapt.o \
	$(BUILD_DIR)/multi_level.o

# main.f90 is appended to SRC below the optional physics include.  Keep
# that separate placement, but make its new module prerequisites explicit.
$(BUILD_DIR)/main.o: \
	$(BUILD_DIR)/parallel_block_build.o \
	$(BUILD_DIR)/parallel_block_mpi.o \
	$(BUILD_DIR)/time_integr.o \
	$(BUILD_DIR)/remap.o

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

$(BUILD_DIR)/parallel_block_primitive_kernels.o: $(BUILD_DIR)/kind.o
$(BUILD_DIR)/parallel_block_columns.o: $(BUILD_DIR)/parallel_block_inverse.o $(BUILD_DIR)/parallel_block_primitive_kernels.o $(BUILD_DIR)/init.o
$(BUILD_DIR)/parallel_block_inverse.o: $(BUILD_DIR)/domain_ops.o
$(BUILD_DIR)/diagnostics.o $(BUILD_DIR)/time_integr.o $(BUILD_DIR)/multi_level.o: $(BUILD_DIR)/parallel_block_columns.o

$(BUILD_DIR)/parallel_block_primitive_walk.o: $(BUILD_DIR)/domain.o
$(BUILD_DIR)/parallel_block_resident_operands.o: $(BUILD_DIR)/parallel_block_inverse.o
$(BUILD_DIR)/parallel_block_primitives.o: $(BUILD_DIR)/parallel_block_resident_operands.o $(BUILD_DIR)/parallel_block_columns.o $(BUILD_DIR)/parallel_block_primitive_walk.o
$(BUILD_DIR)/parallel_block_inverse.o: $(BUILD_DIR)/parallel_block_primitive_walk.o
$(BUILD_DIR)/ops.o $(BUILD_DIR)/multi_level.o: $(BUILD_DIR)/parallel_block_primitives.o

$(BUILD_DIR)/parallel_block_pentagons.o: $(BUILD_DIR)/parallel_block_primitives.o $(BUILD_DIR)/domain_ops.o
$(BUILD_DIR)/ops.o $(BUILD_DIR)/multi_level.o: $(BUILD_DIR)/parallel_block_pentagons.o

$(BUILD_DIR)/parallel_block_gradients.o: $(BUILD_DIR)/parallel_block_primitives.o $(BUILD_DIR)/parallel_block_velocity.o
$(BUILD_DIR)/multi_level.o: $(BUILD_DIR)/parallel_block_gradients.o

$(BUILD_DIR)/parallel_block_mpi.o $(BUILD_DIR)/time_integr.o: $(BUILD_DIR)/parallel_block_gradients.o

$(BUILD_DIR)/parallel_block_program.o: $(BUILD_DIR)/parallel_block_mass.o $(BUILD_DIR)/parallel_block_resident_operands.o

ifneq ($(filter $(BLOCK_DIFFUSION),on check),)
ifneq ($(TEST_CASE),climate)
$(error BLOCK_DIFFUSION currently requires TEST_CASE=climate)
endif
ifeq ($(filter $(BLOCK_PRIMITIVES),on check),)
$(error BLOCK_DIFFUSION requires BLOCK_PRIMITIVES=on/check)
endif
FFLAGS += -DBLOCK_DIFFUSION
endif
ifeq ($(BLOCK_DIFFUSION),check)
FFLAGS += -DBLOCK_DIFFUSION_ORACLE
endif
$(BUILD_DIR)/parallel_block_diffusion.o: $(BUILD_DIR)/parallel_block_program.o $(BUILD_DIR)/parallel_block_columns.o $(BUILD_DIR)/parallel_block_primitive_walk.o
$(BUILD_DIR)/multi_level.o: $(BUILD_DIR)/parallel_block_diffusion.o

$(BUILD_DIR)/parallel_block_program.o: $(BUILD_DIR)/parallel_block_primitive_kernels.o
$(BUILD_DIR)/parallel_block_diffusion.o: $(BUILD_DIR)/utils.o
$(BUILD_DIR)/ops.o: $(BUILD_DIR)/parallel_block_diffusion.o

ifneq ($(filter $(BLOCK_MASS_OWNER),on check),)
ifeq ($(filter $(BLOCK_DIFFUSION),on check),)
$(error BLOCK_MASS_OWNER requires BLOCK_DIFFUSION=on/check)
endif
FFLAGS += -DBLOCK_MASS_OWNER
endif
ifeq ($(BLOCK_MASS_OWNER),check)
FFLAGS += -DBLOCK_MASS_OWNER_ORACLE
endif

$(BUILD_DIR)/parallel_block_primitive_registry.o: $(BUILD_DIR)/domain.o
$(BUILD_DIR)/parallel_block_primitives.o $(BUILD_DIR)/parallel_block_pentagons.o: $(BUILD_DIR)/parallel_block_primitive_registry.o

$(BUILD_DIR)/parallel_block_program.o: $(BUILD_DIR)/parallel_block_velocity.o

ifneq ($(filter $(BLOCK_VELOCITY_OWNER),on check),)
ifeq ($(filter $(BLOCK_MASS_OWNER),on check),)
$(error BLOCK_VELOCITY_OWNER requires BLOCK_MASS_OWNER=on/check)
endif
FFLAGS += -DBLOCK_VELOCITY_OWNER
endif
ifeq ($(BLOCK_VELOCITY_OWNER),check)
FFLAGS += -DBLOCK_VELOCITY_OWNER_ORACLE
endif
$(BUILD_DIR)/parallel_block_definition_exchange.o: $(BUILD_DIR)/domain.o
$(BUILD_DIR)/parallel_block_velocity_owner.o: $(BUILD_DIR)/parallel_block_definition_exchange.o $(BUILD_DIR)/parallel_block_diffusion.o $(BUILD_DIR)/parallel_block_pentagons.o
$(BUILD_DIR)/multi_level.o: $(BUILD_DIR)/parallel_block_velocity_owner.o

$(BUILD_DIR)/parallel_block_gradients.o: $(BUILD_DIR)/parallel_block_velocity_owner.o
$(BUILD_DIR)/utils.o: $(BUILD_DIR)/comm_mpi.o
$(BUILD_DIR)/comm_mpi.o: $(BUILD_DIR)/comm.o $(BUILD_DIR)/domain_ops.o
$(BUILD_DIR)/comm.o: $(BUILD_DIR)/domain.o

ifneq ($(filter $(BLOCK_TENDENCY_INITIAL),on check),)
ifeq ($(filter $(BLOCK_VELOCITY_OWNER),on check),)
$(error BLOCK_TENDENCY_INITIAL requires BLOCK_VELOCITY_OWNER=on/check)
endif
FFLAGS += -DBLOCK_TENDENCY_INITIAL
endif

ifneq ($(filter $(BLOCK_RK),on check),)
ifeq ($(filter $(BLOCK_TENDENCY_INITIAL),on check),)
$(error BLOCK_RK requires BLOCK_TENDENCY_INITIAL=on/check)
endif
FFLAGS += -DBLOCK_RK
endif

$(BUILD_DIR)/parallel_block_mpi.o: $(BUILD_DIR)/parallel_block_diffusion.o

# Simple-physics hydrostatic scratch migration; keep the legacy route available.
PHYSICS_COLUMN_LOCAL ?= off
ifneq ($(filter $(PHYSICS_COLUMN_LOCAL),on check),)
FFLAGS += -DPHYSICS_COLUMN_LOCAL
endif
ifeq ($(PHYSICS_COLUMN_LOCAL),check)
FFLAGS += -DPHYSICS_COLUMN_ORACLE
endif
$(BUILD_DIR)/physics_simple.o: $(BUILD_DIR)/parallel_block_primitive_kernels.o

# Explicit velocity read versions for the serial Simple-physics column walk.
PHYSICS_VELOCITY_PLAN ?= off
ifneq ($(filter $(PHYSICS_VELOCITY_PLAN),on check),)
ifeq ($(filter $(PHYSICS_COLUMN_LOCAL),on check),)
$(error PHYSICS_VELOCITY_PLAN requires PHYSICS_COLUMN_LOCAL=on/check)
endif
FFLAGS += -DPHYSICS_VELOCITY_PLAN
endif
ifeq ($(PHYSICS_VELOCITY_PLAN),check)
FFLAGS += -DPHYSICS_VELOCITY_ORACLE
endif

# Conservative ordered package execution on final block owners.
PHYSICS_OWNER ?= off
ifneq ($(filter $(PHYSICS_OWNER),on check),)
ifeq ($(filter $(PHYSICS_VELOCITY_PLAN),on check),)
$(error PHYSICS_OWNER requires PHYSICS_VELOCITY_PLAN=on/check)
endif
FFLAGS += -DPHYSICS_OWNER
endif
ifeq ($(PHYSICS_OWNER),check)
FFLAGS += -DPHYSICS_OWNER_ORACLE
endif

# Owner-resident physics columns with explicit velocity dependency waves.
PHYSICS_OWNER_BATCH ?= off
ifneq ($(filter $(PHYSICS_OWNER_BATCH),on check),)
ifeq ($(filter $(PHYSICS_OWNER),on check),)
$(error PHYSICS_OWNER_BATCH requires PHYSICS_OWNER=on/check)
endif
FFLAGS += -DPHYSICS_OWNER_BATCH
endif
ifeq ($(PHYSICS_OWNER_BATCH),check)
FFLAGS += -DPHYSICS_OWNER_BATCH_ORACLE
endif

PHYSICS_NATIVE_INPUT ?= off
ifneq ($(filter $(PHYSICS_NATIVE_INPUT),on check),)
ifeq ($(filter $(PHYSICS_OWNER_BATCH),on check),)
$(error PHYSICS_NATIVE_INPUT requires PHYSICS_OWNER_BATCH=on/check)
endif
ifneq ($(BLOCK_RK),on)
$(error PHYSICS_NATIVE_INPUT requires BLOCK_RK=on)
endif
FFLAGS += -DPHYSICS_NATIVE_INPUT
endif
ifeq ($(PHYSICS_NATIVE_INPUT),check)
FFLAGS += -DPHYSICS_NATIVE_INPUT_ORACLE
endif

PHYSICS_NATIVE_STATE ?= off
ifneq ($(filter $(PHYSICS_NATIVE_STATE),on check),)
ifeq ($(filter $(PHYSICS_NATIVE_INPUT),on check),)
$(error PHYSICS_NATIVE_STATE requires PHYSICS_NATIVE_INPUT=on/check)
endif
FFLAGS += -DPHYSICS_NATIVE_STATE
endif
ifeq ($(PHYSICS_NATIVE_STATE),check)
FFLAGS += -DPHYSICS_NATIVE_STATE_ORACLE
endif

SOIL_NATIVE ?= off
ifneq ($(filter $(SOIL_NATIVE),on check),)
ifeq ($(filter $(PHYSICS_NATIVE_STATE),on check),)
$(error SOIL_NATIVE requires PHYSICS_NATIVE_STATE=on/check)
endif
ifneq ($(STAGE203),on)
$(error SOIL_NATIVE requires STAGE203=on)
endif
FFLAGS += -DSOIL_NATIVE
endif
ifeq ($(SOIL_NATIVE),check)
FFLAGS += -DSOIL_NATIVE_ORACLE
endif

# Native numerical topology transition; old/off paths remain available.
ADAPT_NATIVE ?= off
ifneq ($(filter $(ADAPT_NATIVE),on check),)
ifeq ($(filter $(SOIL_NATIVE),on check),)
$(error ADAPT_NATIVE requires SOIL_NATIVE=on/check)
endif
FFLAGS += -DADAPT_NATIVE
endif
ifeq ($(ADAPT_NATIVE),check)
FFLAGS += -DADAPT_NATIVE_ORACLE
endif
$(BUILD_DIR)/mask.o $(BUILD_DIR)/parallel_block_build.o: $(BUILD_DIR)/parallel_block_inverse.o

$(BUILD_DIR)/lnorms.o: $(BUILD_DIR)/parallel_block_inverse.o

# Native between-step numerical lifecycle and explicit I/O adapters.
LIFECYCLE_NATIVE ?= off
ifneq ($(filter $(LIFECYCLE_NATIVE),on check),)
ifeq ($(filter $(ADAPT_NATIVE),on check),)
$(error LIFECYCLE_NATIVE requires ADAPT_NATIVE=on/check)
endif
FFLAGS += -DLIFECYCLE_NATIVE
endif
ifeq ($(LIFECYCLE_NATIVE),check)
FFLAGS += -DLIFECYCLE_NATIVE_ORACLE
endif
$(BUILD_DIR)/remap.o: $(BUILD_DIR)/parallel_block_inverse.o

COEFFICIENT_NATIVE ?= off
ifneq ($(filter $(COEFFICIENT_NATIVE),on check),)
ifeq ($(filter $(LIFECYCLE_NATIVE),on check),)
$(error COEFFICIENT_NATIVE requires LIFECYCLE_NATIVE=on/check)
endif
ifneq ($(BLOCK_RK),on)
$(error COEFFICIENT_NATIVE requires BLOCK_RK=on)
endif
FFLAGS += -DCOEFFICIENT_NATIVE
endif

MIRROR_FREE ?= off
ifeq ($(MIRROR_FREE),on)
ifneq ($(COEFFICIENT_NATIVE),on)
$(error MIRROR_FREE requires COEFFICIENT_NATIVE=on)
endif
ifneq ($(filter check,$(BLOCK_COLUMNS) $(BLOCK_PRIMITIVES) $(BLOCK_DIFFUSION) $(BLOCK_MASS_OWNER) $(BLOCK_VELOCITY_OWNER) $(PHYSICS_COLUMN_LOCAL) $(PHYSICS_VELOCITY_PLAN) $(PHYSICS_OWNER) $(PHYSICS_OWNER_BATCH) $(PHYSICS_NATIVE_INPUT) $(PHYSICS_NATIVE_STATE) $(SOIL_NATIVE) $(ADAPT_NATIVE) $(LIFECYCLE_NATIVE)),)
$(error MIRROR_FREE requires separate oracle builds with legacy views)
endif
FFLAGS += -DMIRROR_FREE
endif
