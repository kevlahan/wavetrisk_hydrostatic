FC = gfortran
FFLAGS = -O2 -Wall -std=f2008

# Source directories
SRCDIR = ../src
CURDIR = .

SRC = \
    ../src/kind.f90 \
    ../src/param_J5.f90 \
    ../src/shared.f90 \
    ../src/coord_arithmetic.f90 \
    ../src/geom.f90 \
    ./vertices2wavetrisk.f90 \

# Object files (automatic replacement of .f90 with .o)
OBJ = $(SRC:.f90=.o)

# Executable name
EXEC = vertices2wavetrisk

# Default rule
all: $(EXEC)

# Linking rule
$(EXEC): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

# Compilation rule
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@ 

# Clean rule
clean:
	rm -f $(CURDIR)/*.o $(CURDIR)/*.mod $(SRCDIR)/*.o $(SRCDIR)/*.mod


