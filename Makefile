# Default general options
TEST_CASE = jet
PARAM     = param_J5
ARCH      = mpi
OPTIM     = 2
F90       = gfortran
MPIF90    = mpif90
BIN_DIR   = bin
BUILD_DIR  = build
LAPACK     = -llapack # link to lapack library
TOPO       = false    # use NCAR topography
PREFIX     = .

# Topgraphy options
TOPO_DIR   = ~/Topo
NETCDF_DIR = /opt/homebrew/Cellar/netcdf-fortran/4.6.1

# AMPI options
CHARM_DIR   = ~/charm
CHARM_BUILD = ucx-linux-x86_64-openpmix-smp
#CHARM_BUILD  = multicore-linux-x86_64 
AMPIF90     = $(CHARM_DIR)/$(CHARM_BUILD)/bin/mpif90.ampi

# Link to test case module file	
$(shell ln -nsf ../test/$(TEST_CASE)/test_case_module.f90 src)
$(shell ln -nsf ../test/$(TEST_CASE)/$(TEST_CASE).f90 src/test.f90)

# Make directories
$(shell mkdir -p $(PREFIX)/$(BUILD_DIR))
$(shell mkdir -p $(PREFIX)/$(BIN_DIR))

vpath %.f90 src

SYSTEM = $(shell uname -a | cut -c 1-6 -)
ifeq ($(SYSTEM),Darwin)
  MACHINE = mac
  LIB    += -L/opt/homebrew/opt/lapack/lib
else
  MACHINE = $(shell uname -n | sed -e "s/[^a-z].*//")

  ifeq ($(MACHINE),$(filter $(MACHINE),orc bul gra nia))
    LAPACK = -lopenblas # use module load openblas
  endif
endif

ifeq ($(F90),ifort)
  FLAGS_COMP = -O$(OPTIM) -c -Isrc/ppr -cpp -diag-disable 8291
  FLAGS_LINK = -O$(OPTIM)
  ifeq ($(ARCH),mpi) # problem with -module when using AMPI
    FLAGS_COMP += -module $(BUILD_DIR)
    FLAGS_LINK += -module $(BUILD_DIR)
  endif
else
  FLAGS_COMP = -O$(OPTIM) -c -J$(BUILD_DIR) -cpp
  FLAGS_LINK = -O$(OPTIM) -J $(BUILD_DIR)
endif

ifeq ($(OPTIM),0)
  ifeq ($(F90),ifort)
    FLAGS_COMP += -g -traceback
  else
    FLAGS_COMP += -g -fbacktrace -fcheck=all
  endif
endif

ifeq ($(ARCH),ser)
  COMPILER = $(F90)
  PROC = ser
else
   PROC = mpi
   FLAGS_COMP += -DMPI
   FLAGS_LINK += -DMPI
  ifeq ($(ARCH),mpi)
    COMPILER = $(MPIF90)
  else
    ARCH        = mpi
    F90         = $(AMPIF90)
    COMPILER    = $(AMPIF90)
    FLAGS_COMP += -DAMPI -pieglobals
    FLAGS_LINK += -DAMPI -pieglobals
  endif
endif
LINKER = $(COMPILER)

LIBS   = $(LAPACK)

ifeq ($(TEST_CASE), spherical_harmonics) # add shtools and supporting libraries (MUST use gfortran/openmpi)
  ifeq ($(MACHINE),$(filter $(MACHINE),orc bul gra nia))
    # module load gcc openmpi fftw mkl
    SHTOOLSLIBPATH = $(HOME)/SHTOOLS-4.7.1/lib
    SHTOOLSMODPATH = $(HOME)/SHTOOLS-4.7.1/include
    LAPACK = -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core
  else ifeq ($(MACHINE), mac)
    SHTOOLSMODPATH = /opt/homebrew/include
    SHTOOLSLIBPATH = /opt/homebrew/lib
    LAPACK = -llapack
  else
    SHTOOLSMODPATH = /usr/local/include
    SHTOOLSLIBPATH = /usr/local/lib
    LAPACK = -llapack
  endif
  LIBS += -L$(SHTOOLSLIBPATH) -lSHTOOLS -lfftw3 -lm $(LAPACK)
  FLAGS_COMP += -I$(SHTOOLSMODPATH) -m64 -fPIC
endif

ifeq ($(TOPO), true)
  FLAGS_COMP += -I$(NETCDF_DIR)/include
  LIBS       += -L$(NETCDF_DIR)/lib -lnetcdff
endif

SRC = $(PARAM).f90 shared.f90 coord_arithmetic.f90 sphere.f90  patch.f90 dyn_array.f90 \
base_$(PROC).f90 spline.f90 domain.f90 init.f90 comm.f90 comm_$(PROC).f90 utils.f90 \
projection.f90 equation_of_state.f90 wavelet.f90 lnorms.f90 mask.f90 refine_patch.f90 smooth.f90 ops.f90 \
multi_level.f90 adapt.f90 lin_solve.f90 barotropic_2d.f90 time_integr.f90 vert_diffusion.f90 io.f90 \
remap.f90 main.f90 std_atm_profile.f90

ifeq ($(TOPO), true)
  SRC += topo_grid_descriptor.f90
endif

SRC += test_case_module.f90 test.f90

OBJ = $(patsubst %.f90,$(BUILD_DIR)/%.o,$(SRC))

$(PREFIX)/$(BIN_DIR)/$(TEST_CASE): $(OBJ)
	$(LINKER) $(FLAGS_LINK) -o $@ $^ $(LIBS) 

$(BUILD_DIR)/%.o: %.f90 shared.f90 $(PARAM).f90
	$(COMPILER) $(FLAGS_COMP) $< -o $@ 

clean:
	\rm -f $(BUILD_DIR)/* src/test_case_module.f90 src/test.f90
