# options: ser mpi mpi-lb
TEST_CASE  = DCMIP2012c4
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
      base_$(ARCH).f90 dgesv.f90 domain.f90 init.f90 comm.f90 comm_$(ARCH).f90 \
      wavelet.f90 lnorms.f90 mask.f90 refine_patch.f90 smooth.f90 test_case_module.f90 ops.f90 multi_level.f90   \
      adapt.f90 lin_solve.f90 barotropic_2d.f90 io.f90 time_integr.f90 remap.f90 main.f90 
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
  MPIF90 = mpif90
  OPTIM  = -O2
  #FLAGS  =  $(OPTIM) -g -trace -profile=vtmc -module $(BUILD_DIR) -Isrc/ppr -cpp -diag-disable 8291	
  FLAGS  =  $(OPTIM) -traceback -module $(BUILD_DIR) -Isrc/ppr -cpp -diag-disable 8291
  LIBS   = 
else # gfortran as default
  F90    = gfortran
  MPIF90 = mpif90
  OPTIM  = -O2
  FLAGS  = $(OPTIM) -J$(BUILD_DIR) -cpp -fbacktrace -fcheck=all
  LIBS   = 
endif

ifeq ($(TEST_CASE), spherical_harmonics) # add shtools and supporting libraries (must use gfortran)
  F90    = gfortran
  MPIF90 = mpif90
  OPTIM  = -O2
  ifeq ($(MACHINE),$(filter $(MACHINE),orc bul gra nia))
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
