# known to compile with gfortran version 7, ifort version 17, pgfortran version 18
# Examples
# F90=ifort FFLAGS=-fast make clean all
# F90=gfortran F90FLAGS=-ffree-line-length-none make clean all
# F90=pgfortran F90FLAGS="-fast -Minfo=loop" make clean all

# if XCodeML F_Front is installed, set FFRONT=F_Front to parse code into XML
FFRONT ?= true

F90 ?= gfortran
F90FLAGS ?= -ffree-line-length-none -g -fbacktrace -fcheck=all -ffpe-trap=zero,overflow,invalid -std=f2003 -freal-4-real-8 #added real-real(8) flag

NAMES = $(basename $(notdir $(wildcard physics/*.F90)))
OBJECTS = $(addprefix obj/, $(addsuffix .o, $(NAMES)))
#EXE = driver/driver.exe #Compiling on windows
EXE = driver/driver #Compiling on unix

all : obj $(OBJECTS) lib/libphyparam.so $(EXE)

test: all
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/lib driver/driver.exe

py : clean all
	bash/extract_cython.sh
	cd python ; rm -rf build *.c ; python setup.py build_ext --inplace
obj :
	@rm -rf obj include lib xml tmp *~ */*~
	@mkdir obj include lib xml tmp outputs

clean : nice
	@rm -rf obj include lib xml tmp outputs *~ */*~
	@mkdir obj include lib xml tmp outputs
	@bash/concatenate_all_code.sh $(NAMES) > tmp/all_code
	@bash/makedeps.sh tmp/all_code >  Makefile.inc

nice:
	bash/beautify.sh physics/*.F90 DYNAMICO/*.F90 driver/driver.F90

doc:	all
	cd tmp ; doxygen ../bash/Doxyfile

%.so : $(OBJECTS)
	$(F90) -shared $^ -o $@ 

#%.exe : %.F90 lib/libphyparam.so #If compiling on windows
$(EXE) : $(addsuffix .F90,$(EXE)) lib/libphyparam.so # Compiling on unix
	ar cr driver/libphyparam.a obj/*.o
	$(F90) -I include -L driver -fPIC $(F90FLAGS) $< -lphyparam $(LDFLAGS) -o $@
	mv *.mod driver/

SRC = $(basename $(notdir $@))
obj/%.o: physics/%.F90
	$(F90) -I include -fPIC $(F90FLAGS) -c $<
	@mv $(SRC).mod include
	@mv $(SRC).o obj
	@cpp -DXCODEML $< > tmp/$(SRC).f90
	@$(FFRONT) -M tmp tmp/$(SRC).f90 > xml/$(SRC).xml

# dependencies between modules
$(OBJECTS) : obj/logging.o
-include Makefile.inc
