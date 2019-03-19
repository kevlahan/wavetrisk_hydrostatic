# options: ser mpi mpi-lb
TEST_CASE  = DCMIP2012c4
OPTIM      = -O2
ARCH       = mpi-lb
PARAM      = param_J4
GEOM       = sphere
ARRAYS     = dyn_array
BUILD_DIR  = build

PREFIX = .

# Remove files associated with previous test case
$(shell \rm $(BUILD_DIR)/test_case_module.o $(BUILD_DIR)/test_case_mod.mod)

# Link to test case module file	
$(shell ln -nsf ../test/$(TEST_CASE)/test_case_module.f90 src/.)

vpath %.f90 src
SRC = $(PARAM).f90 shared.f90 $(GEOM).f90 patch.f90 $(ARRAYS).f90 \
      base_$(ARCH).f90 domain.f90 init.f90 comm.f90 comm_$(ARCH).f90 \
      wavelet.f90 lnorms.f90 mask.f90 refine_patch.f90 smooth.f90 test_case_module.f90 ops.f90 \
      multi_level.f90 adapt.f90 io.f90 time_integr.f90 remap.f90 main.f90
OBJ = $(patsubst %.f90,$(BUILD_DIR)/%.o,$(SRC))

SYSTEM = $(shell uname -a | cut -c 1-6 -)
ifeq ($(SYSTEM),Darwin)
  MACHINE = mac
else
  MACHINE = $(shell uname -n | sed -e "s/[^a-z].*//")
endif

# lapack/mkl used only to solve 6x6 linear system
ifeq ($(MACHINE),if)
  F90    = gfortran	
  MPIF90 = mpif90
  LIBS   = -llapack
  FLAGS  = $(OPTIM) -J$(BUILD_DIR) -fbacktrace -fcheck=all
else ifeq ($(MACHINE),$(filter $(MACHINE),orc bul gra, nia))
  # Need module load intel; module load openmpi
  F90    = ifort	
  MPIF90 = mpif90
  LIBS   = -mkl
  FLAGS  = -traceback -module $(BUILD_DIR) -diag-disable 8291
# else ifeq ($(MACHINE),$(filter $(MACHINE),nia))
#   # Need to load: module load gcc/7.3.0 openmpi/3.1.0 mkl/2018.2                                                                                                        
#   F90    = gfortran
#   MPIF90 = mpif90
#   LIBS   = -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
#   FLAGS  = $(OPTIM) -J$(BUILD_DIR) -fbacktrace -fcheck=all -I${MKLROOT}/include
else ifeq ($(MACHINE),$(filter $(MACHINE),login)) # occigen
  # Need to load: module load openmpi/gnu/2.0.2 mkl/18.1                                                                                                           
  F90    = gfortran
  MPIF90 = mpif90
  LIBS   = -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
  FLAGS  = $(OPTIM) -J$(BUILD_DIR) -fbacktrace -fcheck=all -B /opt/software/common/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64
else ifeq ($(MACHINE),mac)
  F90    = gfortran
  LIBS   = -llapack	
  MPIF90 = mpif90
  FLAGS  = $(OPTIM) -J$(BUILD_DIR) -fbacktrace -fcheck=all
else ifeq ($(MACHINE),froggy)
  F90    = gfortran
  MPIF90 = mpif90
  FLAGS  = $(OPTIM) -J$(BUILD_DIR) -fbacktrace -fcheck=all
else # try gfortran as default
  F90    = gfortran
  MPIF90 = mpif90
  FLAGS  = $(OPTIM) -J$(BUILD_DIR) -fbacktrace -fcheck=all
endif

ifeq ($(ARCH),ser)
  COMPILER = $(F90)
else
  COMPILER = $(MPIF90)
endif
LINKER = $(COMPILER)

$(PREFIX)/bin/$(TEST_CASE): $(OBJ) test/$(TEST_CASE)/$(TEST_CASE).f90
	mkdir -p $(PREFIX)/bin
	$(LINKER) $(FLAGS) -o $@ $^ $(LIBS)

$(BUILD_DIR)/%.o: %.f90 shared.f90 $(PARAM).f90
	$(COMPILER) -c $< -o $@ $(FLAGS) 

clean:
	rm -f $(BUILD_DIR)/*
