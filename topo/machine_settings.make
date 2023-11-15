DEBUG=FALSE

#FC = $(shell nc-config --fc)
#LDFLAGS = $(shell nc-config --flibs)
#FFLAGS = -c $(shell nc-config --fflags)

FC      = $(shell nf-config --fc)
LDFLAGS = -L/opt/homebrew/Cellar/netcdf-fortran/4.6.1/lib -lnetcdff 
FFLAGS  = -c $(shell nf-config --fflags)

#------------------------------------------------------------------------
# Gfortran
#------------------------------------------------------------------------
ifeq ($(findstring gfortran, $(FC)),gfortran)
    ifneq ($(findstring pgfortran, $(FC)),pgfortran)
      ifeq ($(DEBUG),TRUE)
        FFLAGS += -Wall -fbacktrace -fbounds-check -fno-range-check 
      else
        FFLAGS += -O -fallow-argument-mismatch
    endif
  endif
endif

#------------------------------------------------------------------------
# NAG
#------------------------------------------------------------------------
ifeq ($(findstring nagfor, $(FC)),nagfor)
  ifeq ($(DEBUG),TRUE)
    FFLAGS += -g -C
  else
    FFLAGS += -O
  endif

endif

#------------------------------------------------------------------------
# PGF95
#------------------------------------------------------------------------
ifeq ($(findstring pgf95, $(FC)),pgf95)
  FFLAGS += -Mlarge_arrays

  ifeq ($(DEBUG),TRUE)
    FFLAGS += -g -Mbounds -traceback -Mchkfpstk
  else
    FFLAGS += -O
  endif
endif

export LDFLAGS
export FFLAGS
export FC
