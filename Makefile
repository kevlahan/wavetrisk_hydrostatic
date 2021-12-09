# Default options
TEST_CASE = jet
ARCH      = mpi-lb
PARAM     = param_J5
GEOM      = sphere
ARRAYS    = dyn_array
BUILD_DIR = build
MPIF90    = mpi
F90       = gfortran
OPTIM     = 2
DEBUG     = 
AMPIF90   = ~/charm/bin/mpif90.ampi
LIBS      = 
PREFIX    = .

# Link to test case module file	
$(shell ln -nsf ../test/$(TEST_CASE)/test_case_module.f90 src)
$(shell ln -nsf ../test/$(TEST_CASE)/$(TEST_CASE).f90 src/test.f90)

# Make directories
$(shell mkdir -p $(PREFIX)/$(BUILD_DIR))
$(shell mkdir -p $(PREFIX)/bin)

vpath %.f90 src

SYSTEM = $(shell uname -a | cut -c 1-6 -)
ifeq ($(SYSTEM),Darwin)
  MACHINE = mac
else
  MACHINE = $(shell uname -n | sed -e "s/[^a-z].*//")
endif

ifeq ($(MACHINE),$(filter $(MACHINE),orc bul gra nia))
# Need: module load NiaEnv/.2021a cmake intel intelmpi ucx
  F90 = ifort
endif

ifeq ($(OPTIM),0)
  ifeq ($(F90),ifort)
    DEBUG = -g -traceback 
  else
    DEBUG = -ggdb3 -fbacktrace -fcheck=all
  endif
endif

ifeq ($(F90),ifort)
  FLAGS_COMP = -O$(OPTIM) $(DEBUG) -c -Isrc/ppr -cpp -diag-disable 8291
  FLAGS_LINK = -O$(OPTIM)
  ifeq ($(MPIF90),mpi) # problem with -module when using AMPI
    FLAGS_COMP += -module $(BUILD_DIR)
    FLAGS_LINK += -module $(BUILD_DIR)
  endif
else
# Use -fallow-argument-mismatch to deal with mpi argument mismatch bug in gcc 10.1
  FLAGS_COMP = -O$(OPTIM) $(DEBUG) -c -J$(BUILD_DIR) -cpp #-fallow-argument-mismatch
  FLAGS_LINK = -O$(OPTIM) -J$(BUILD_DIR)
endif

SRC = $(PARAM).f90 shared.f90 $(GEOM).f90 patch.f90 $(ARRAYS).f90 \
base_$(ARCH).f90 dgesv.f90 dgtsv.f90 spline.f90 domain.f90 init.f90 utils.f90 comm.f90 comm_$(ARCH).f90  projection.f90 equation_of_state.f90 \
wavelet.f90 lnorms.f90 mask.f90 refine_patch.f90 smooth.f90 ops.f90 multi_level.f90 adapt.f90 lin_solve.f90  \
barotropic_2d.f90 time_integr.f90 vert_diffusion.f90 lateral_diffusion.f90 io.f90 remap.f90 main.f90 test_case_module.f90 test.f90 

OBJ = $(patsubst %.f90,$(BUILD_DIR)/%.o,$(SRC))

ifeq ($(TEST_CASE), spherical_harmonics) # add shtools and supporting libraries (MUST use gfortran/openmpi)
  ifeq ($(MACHINE),$(filter $(MACHINE),orc bul gra nia))
    # need to do:
    # module load gcc openmpi fftw mkl
    SHTOOLSLIBPATH = /home/k/kevlahan/kevlahan/SHTOOLS-4.7.1/lib
    SHTOOLSMODPATH = /home/k/kevlahan/kevlahan/SHTOOLS-4.7.1/include
    BLAS   =-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lmkl_blas95_lp64
    FLAGS  = -O$(OPTIM) -J$(BUILD_DIR) -cpp -fbacktrace -fcheck=all -std=gnu -ffast-math -I$(SHTOOLSMODPATH) -m64 -fPIC
    LIBS   =-L$(SHTOOLSLIBPATH) -lSHTOOLS -lfftw3 -lm $(BLAS) $(LAPACK)
  else 
    SHTOOLSMODPATH = /usr/local/include
    SHTOOLSLIBPATH = /usr/local/lib
    FLAGS += -I$(SHTOOLSMODPATH) -m64 -fPIC 
    LIBS  += -L$(SHTOOLSLIBPATH) -lSHTOOLS -lfftw3 -lm -lblas -llapack
  endif
endif

ifeq ($(ARCH),ser)
  COMPILER = $(F90)	
else
   FLAGS_COMP += -DMPI
   FLAGS_LINK += -DMPI
  ifeq ($(MPIF90),mpi)
    COMPILER = mpif90
  else
    ARCH        = mpi
    F90         = $(AMPIF90)
    COMPILER    = $(AMPIF90)
    FLAGS_COMP += -DAMPI -pieglobals
    FLAGS_LINK += -DAMPI -pieglobals
  endif
endif
LINKER = $(COMPILER)

$(PREFIX)/bin/$(TEST_CASE): $(OBJ) 
	$(LINKER) $(FLAGS_LINK) -o $@ $^ $(LIBS) 

$(BUILD_DIR)/%.o: %.f90 shared.f90 $(PARAM).f90
	$(COMPILER) $(FLAGS_COMP) $< -o $@ 

clean:
	\rm -f $(BUILD_DIR)/* src/test_case_module.f90 src/test.f90
