# Default options
TEST_CASE  = jet
ARCH       = mpi-lb
PARAM      = param_J5
GEOM       = sphere
ARRAYS     = dyn_array
BUILD_DIR  = build
MPIF90     = mpi
AMPIF90    = ~/charm/bin/mpif90.ampi
F90        = gfortran
OPTIM      = -O2
FLAGS      = $(OPTIM) -J$(BUILD_DIR) -cpp -fbacktrace -fcheck=all
LIBS       = 

PREFIX = .

ifeq ($(ARCH),ser)
  COMPILER = $(F90)
else
  ifeq ($(MPIF90),mpi)
    COMPILER = mpif90
  else
    ARCH     = mpi
    F90      = $(AMPIF90)
    COMPILER = $(AMPIF90)
    FLAGS    = $(OPTIM) -J$(BUILD_DIR) -cpp -DAMPI -pieglobals 
  endif
endif
LINKER = $(COMPILER)

# Remove files associated with previous test case
$(shell \rm $(BUILD_DIR)/test_case_module.o $(BUILD_DIR)/test_case_mod.mod)

# Link to test case module file	
$(shell ln -nsf ../test/$(TEST_CASE)/test_case_module.f90 src/.)

vpath %.f90 src

SRC = $(PARAM).f90 shared.f90 $(GEOM).f90 patch.f90 $(ARRAYS).f90 \
base_$(ARCH).f90 dgesv.f90 dgtsv.f90 spline.f90 domain.f90 init.f90 utils.f90 comm.f90 comm_$(ARCH).f90  projection.f90 equation_of_state.f90 \
wavelet.f90 lnorms.f90 mask.f90 refine_patch.f90 smooth.f90 ops.f90 multi_level.f90 adapt.f90 lin_solve.f90  \
barotropic_2d.f90 time_integr.f90 vert_diffusion.f90 lateral_diffusion.f90 io.f90 remap.f90 main.f90 test_case_module.f90

OBJ = $(patsubst %.f90,$(BUILD_DIR)/%.o,$(SRC))

SYSTEM = $(shell uname -a | cut -c 1-6 -)
ifeq ($(SYSTEM),Darwin)
  MACHINE = mac
else
  MACHINE = $(shell uname -n | sed -e "s/[^a-z].*//")
endif 

ifeq ($(MACHINE),$(filter $(MACHINE),orc bul gra nia))
  # Need module load intel; module load intelmpi
  F90    = ifort	
  FLAGS  = $(OPTIM) -traceback -module $(BUILD_DIR) -Isrc/ppr -cpp -diag-disable 8291
endif

ifeq ($(TEST_CASE), spherical_harmonics) # add shtools and supporting libraries (MUST use gfortran/openmpi)
  ifeq ($(MACHINE),$(filter $(MACHINE),orc bul gra nia))
    # need to do:
    # module load gcc openmpi fftw mkl
    SHTOOLSLIBPATH = /home/k/kevlahan/kevlahan/SHTOOLS-4.7.1/lib
    SHTOOLSMODPATH = /home/k/kevlahan/kevlahan/SHTOOLS-4.7.1/include
    BLAS   =-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lmkl_blas95_lp64
    FLAGS  = $(OPTIM) -J$(BUILD_DIR) -cpp -fbacktrace -fcheck=all -std=gnu -ffast-math -I$(SHTOOLSMODPATH) -m64 -fPIC
    LIBS   =-L$(SHTOOLSLIBPATH) -lSHTOOLS -lfftw3 -lm $(BLAS) $(LAPACK)
  else 
    SHTOOLSMODPATH = /usr/local/include
    SHTOOLSLIBPATH = /usr/local/lib
    FLAGS += -I$(SHTOOLSMODPATH) -m64 -fPIC 
    LIBS  += -L$(SHTOOLSLIBPATH) -lSHTOOLS -lfftw3 -lm -lblas -llapack
  endif
endif

$(PREFIX)/bin/$(TEST_CASE): $(OBJ) test/$(TEST_CASE)/$(TEST_CASE).f90
	mkdir -p $(PREFIX)/bin
	$(LINKER) $(FLAGS) -o $@ $^ $(LIBS)

$(BUILD_DIR)/%.o: %.f90 shared.f90 $(PARAM).f90
	$(COMPILER) $(FLAGS) -c $< -o $@ 

clean:
	\rm -f $(BUILD_DIR)/*
