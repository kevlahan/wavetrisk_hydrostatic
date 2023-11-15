FC      = $(shell nf-config --fc)
LDFLAGS = -L/opt/homebrew/Cellar/netcdf-fortran/4.6.1/lib -lnetcdff 
FFLAGS  = -c $(shell nf-config --fflags) -O -fallow-argument-mismatch
EXEDIR  = ../bin
EXENAME = cube_to_target

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
