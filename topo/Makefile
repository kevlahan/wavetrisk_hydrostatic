FC      = $(shell nf-config --fc)
LDFLAGS = -lnetcdff 
FFLAGS  = -c $(shell nf-config --fflags) -O 
EXEDIR  = ../bin
EXENAME = cube_to_target

# Identify shell
UNAME_S := $(shell uname -s)

# =========================
# Machine / platform tweaks
# =========================
ifeq ($(UNAME_S),Darwin)
MACHINE := mac
# Homebrew lapack path (adjust if you want)
LDLIBS     += -L/opt/homebrew/opt/lapack/lib
NETCDF_DIR ?= /opt/homebrew/Cellar/netcdf-fortran/4.6.2
LDFLAGS    += -L$(NETCDF_DIR)/lib
else
MACHINE := $(shell scontrol show config 2>/dev/null | \
awk -F= '/^ClusterName/{sub(/^[[:space:]]+/,"",$$2); sub(/[[:space:]]+$$/,"",$$2); print tolower($$2); exit}')

$(info Machine: $(MACHINE))

# bbcluster2
ifeq ($(MACHINE),bbcluster2)
ifdef NETLIB_LAPACK_ROOT
LDFLAGS  += -L$(NETLIB_LAPACK_ROOT)/lib64
endif
ifdef NETCDF_FORTRAN_ROOT
LDFLAGS  += -L$(NETCDF_FORTRAN_ROOT)/lib -Wl,-rpath,$(NETCDF_FORTRAN_ROOT)/lib
endif
endif
endif

.SUFFIXES: .F90 .o

.F90.o:
	$(FC) $(FFLAGS) $< $(LDFLAGS)

OBJS := reconstruct.o remap.o shr_kind_mod.o shared_vars.o rot.o smooth_topo_cube.o ridge_utils.o ridge_ana.o cube_to_target.o f90getopt.o

$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	\rm -f $(OBJS)  *.mod 

cube_to_target.o: shr_kind_mod.o remap.o reconstruct.o shared_vars.o f90getopt.o
reconstruct.o: remap.o 
remap.o      : shr_kind_mod.o
shared_vars.o: shr_kind_mod.o
ridge_ana.o: ridge_utils.o
