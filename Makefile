# Default general options
TEST_CASE     = climate
PARAM         = param_J6
ARCH          = mpi
OPTIM         = 2
COMPILER_TYPE = gnu
MPIF90        = mpif90
BIN_DIR       = bin
BUILD_DIR     = build
LAPACK        = -llapack
NETCDF        = -lnetcdff
TOPO          = false
ifeq ($(TEST_CASE), make_NCAR_topo)
 TOPO = true
endif

LIBS = $(LAPACK)

PHYSICS = false
ifeq ($(TEST_CASE), climate)
 PHYSICS = true
endif
ifeq ($(TEST_CASE), spherical_harmonics)
 PHYSICS = true
endif

PREFIX     = .

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
 LIBS   += -L/opt/homebrew/opt/lapack/lib
 ifeq ($(TOPO), true)
  NETCDF_DIR  = /opt/homebrew/Cellar/netcdf-fortran/4.6.2
  FLAGS_COMP += -I$(NETCDF_DIR)/include
  LIBS       += -L$(NETCDF_DIR)/lib
 endif
else
 MACHINE = $(shell uname -n | sed -e "s/[^a-z].*//")

 ifeq ($(MACHINE),$(filter $(MACHINE), orc bul gra nia narval)) # module load StdEnv netcdf 
  LAPACK = -lflexiblas  # module load flexiblas
 endif

 ifeq ($(MACHINE),$(filter $(MACHINE), bbserv))
  FLAGS_COMP += -I$(NETLIB_LAPACK_ROOT)/include
  LIBS       += -L$(NETLIB_LAPACK_ROOT)/lib64
  ifeq ($(TOPO), true)
   FLAGS_COMP += -I$(NETCDF_FORTRAN_ROOT)/include
   LIBS       += -L$(NETCDF_FORTRAN_ROOT)/lib
  endif
 endif
endif

ifeq ($(TOPO), true)
  LIBS += $(NETCDF)  # bbserv: module load netcdf netcdf-fortran
endif

ifeq ($(COMPILER_TYPE),gnu)
 F90 = gfortran
 FLAGS_COMP += -O$(OPTIM) -c -J$(BUILD_DIR) -cpp -fallow-argument-mismatch 
else ifeq ($(COMPILER_TYPE),amd)
 F90 = flang
 FLAGS_COMP += -O$(OPTIM) -c -module $(BUILD_DIR) -cpp
else ifeq ($(COMPILER_TYPE),intel)
 F90 = ifort
 FLAGS_COMP += -O$(OPTIM) -c -Isrc/ppr -cpp -diag-disable 8291
endif
FLAGS_LINK += -O$(OPTIM)

ifeq ($(OPTIM),0)
 ifeq ($(COMPILER_TYPE),intel)
   FLAGS_COMP += -g -traceback
 else
   FLAGS_COMP += -g -fbacktrace -fcheck=all
 endif
endif

ifeq ($(ARCH),ser)
 COMPILER = $(F90)
 PROC     = ser
else
 PROC        = mpi
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

ifeq ($(PHYSICS), true)
 FLAGS_COMP += -DPHYSICS
 FLAGS_LINK += -DPHYSICS
endif

LINKER += $(COMPILER)
LIBS   += $(LAPACK)

ifeq ($(TEST_CASE), spherical_harmonics) # add shtools and supporting libraries (MUST use gfortran/openmpi)
 ifeq ($(MACHINE),$(filter $(MACHINE),orc bul gra nia))
  # module load fftw
  SHTOOLSLIBPATH = $(HOME)/SHTOOLS-4.7.1/lib
  SHTOOLSMODPATH = $(HOME)/SHTOOLS-4.7.1/include
 else ifeq ($(MACHINE), mac)
  SHTOOLSMODPATH = /opt/homebrew/include
  SHTOOLSLIBPATH = /opt/homebrew/lib
 else
  SHTOOLSMODPATH = /usr/local/include
  SHTOOLSLIBPATH = /usr/local/lib
 endif
 LIBS       += -L$(SHTOOLSLIBPATH) -lSHTOOLS -lfftw3 -lm 
 FLAGS_COMP += -I$(SHTOOLSMODPATH) -m64 -fPIC
endif

SRC = kind.f90 $(PARAM).f90 shared.f90 coord_arithmetic.f90 calendar.f90 sphere.f90  patch.f90 dyn_array.f90 \
base_$(PROC).f90 spline.f90 domain.f90 domain_ops.f90 init.f90 comm.f90 comm_$(PROC).f90 utils.f90 \
projection.f90 equation_of_state.f90 wavelet.f90 lnorms.f90 mask.f90 refine_patch.f90 smooth.f90 ops.f90 \
multi_level.f90 adapt.f90 lin_solve.f90 barotropic_2d.f90 time_integr.f90 io.f90 vert_diffusion.f90 io_vtk.f90 \
remap.f90 std_atm_profile.f90 sso.f90

ifeq ($(TOPO), true)
 SRC += topo_grid_descriptor.f90
endif

ifeq ($(PHYSICS), true)
 SIMPLEPHYSMODPATH = src/physics/simple_physics/phyparam/include
 FLAGS_COMP       += -I$(SIMPLEPHYSMODPATH) # mod 
 PHYSLIB_PATH      = src/physics/simple_physics/phyparam/driver
 LIBS             += -L$(PHYSLIB_PATH) -lphyparam
 -include src/physics/Makefile.inc
endif

SRC += main.f90 test_case_module.f90 test.f90

OBJ = $(patsubst %.f90,$(BUILD_DIR)/%.o,$(SRC))

$(PREFIX)/$(BIN_DIR)/$(TEST_CASE): $(OBJ)
	$(LINKER) $(FLAGS_LINK) -o $@ $^ $(LIBS) 

$(BUILD_DIR)/%.o: %.f90 shared.f90 $(PARAM).f90
	$(COMPILER) $(FLAGS_COMP) $< -o $@ 

phys_package:
	@echo "Compiling Physics Package" 
	@$(MAKE) -C src/physics/simple_physics/phyparam F90=mpif90

topography:
	@echo "Compiling NCAR topography package"
	@$(MAKE) -C topo 
	@$(MAKE) -C topo clean

clean:
	\rm -f $(BUILD_DIR)/* src/test_case_module.f90 src/test.f90

ifeq ($(PHYSICS), true)
	$(MAKE) -C src/physics/simple_physics/phyparam clean
endif
