TEST_CASE = tsunami
# options: ser mpi mpi-lb
ARCH = ser
PARAM = param
GEOM = sphere
ARRAYS = dyn_array

BUILD_DIR = build
OPTIM_FLAGS = -O3 
MPIF90 = mpif90
PREFIX = .

vpath %.f90 src
SRC = $(PARAM).f90 shared.f90 $(GEOM).f90 patch.f90 $(ARRAYS).f90 \
      base_$(ARCH).f90 domain.f90 init.f90 comm.f90 comm_$(ARCH).f90 \
      wavelet.f90 mask.f90 refine_patch.f90 viscous.f90 ops.f90 \
      multi_level.f90 adapt.f90 smooth.f90 io.f90 time_integr.f90 main.f90
OBJ = $(patsubst %.f90,$(BUILD_DIR)/%.o,$(SRC))

INTEL_FLAGS = -module $(BUILD_DIR) #-traceback -check bounds
GNU_FLAGS = -J$(BUILD_DIR)

MACHINE = $(shell uname -n | sed -e "s/[^a-z].*//")

# lapack/mkl used only to solve 6x6 eq. sys.
ifeq ($(MACHINE),if)
  VENDOR = gnu
  LIBS = -llapack
else
 ifeq ($(MACHINE),orc)
   VENDOR = intel
   LIBS = -llapack
 else
  ifeq ($(MACHINE),req)
    VENDOR = path
    LIBS = -llapack
  else
   ifeq ($(MACHINE),gpc)
     VENDOR = intel
     LIBS = -mkl
   else # try gfortran and liblapack as default
     VENDOR = gnu
     LIBS = -llapack
   endif
  endif
 endif
endif

ifeq ($(VENDOR),gnu)
    FLAGS = $(OPTIM_FLAGS) $(GNU_FLAGS)
    F90 = gfortran
else
 ifeq ($(VENDOR),intel)
    FLAGS = $(OPTIM_FLAGS) $(INTEL_FLAGS)
    F90 = ifort
 else # default puts all mod files in current directory
    FLAGS = $(OPTIM_FLAGS)
 endif
endif

LIBMETIS = -lmetis
ifeq ($(ARCH),ser)
  COMPILER = $(F90)
else
  COMPILER = $(MPIF90)
  ifeq ($(ARCH),mpi-metis)
    LIBS += $(LIBMETIS)
  endif
endif
LINKER = $(COMPILER)

$(PREFIX)/bin/$(TEST_CASE): $(OBJ) test/$(TEST_CASE)/$(TEST_CASE).f90
	mkdir -p $(PREFIX)/bin
	$(LINKER) $(FLAGS) -o $@ $^ $(LIBS)

$(BUILD_DIR)/%.o: %.f90 shared.f90 $(PARAM).f90
	$(COMPILER) -c $< -o $@ $(FLAGS)

clean:
	rm -f $(BUILD_DIR)/*
