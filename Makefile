# makefile for OpTool
#

# The compiler and linker
FC	  = gfortran
LINKER	  = gfortran

ifeq ($(ifort),true)
  FC	  = ifort
  LINKER  = ifort
endif

# Multicore support
ifeq ($(multi),true)
    ifeq ($(ifort),true)
	MULTICORE = -openmp -fp-model strict
    else
	MULTICORE = -fopenmp
    endif
endif

# Debugging flags
ifeq ($(debug),true)
  ifeq ($(ifort),true)
    DEBUGGING = -check all -traceback -check bounds -O0 -g -check -fpe1
  else	
    DEBUGGING = -fbounds-check -fbacktrace
  endif
endif

# CFITSIO support
ifeq ($(fitsio),true)
  FLAG_FITS		= -DUSE_FITSIO
  LIBS_FITS		= -lcfitsio -L/usr/local/lib/
endif

# Platform specific compilation options
ifeq ($(ifort),true)
  FLAG_ALL      = -O3 -g -extend-source -zero -prec-div $(DEBUGGING) $(MULTICORE)
  FLAG_LINUX    = -xHOST -fpp
  FLAG_MAC      = -xHOST -opt-prefetch -static-intel -fpp
else
  FLAG_ALL      = -O3 -g -fdefault-double-8 -fdefault-real-8 $(DEBUGGING) $(MULTICORE)
  FLAG_LINUX    = -cpp
  FLAG_MAC      = -m64 -cpp
endif

# Combine the flags
FFLAGS  = $(FLAG_ALL) $(FLAG_LINUX) $(FLAG_FITS)
LDFLAGS = $(FLAG_ALL) $(FLAG_LINUX) $(FLAG_FITS)
LIBS    = $(LIBS_FITS)

# Object files to link
OBJS	= optool.o \
	  optool_extra.o \
	  ref_ind.o \
	  dmilay_f95.o

# Program name and install location
PROGRAM       = optool
DEST	      = ${HOME}/bin

# make actions 
all:		$(PROGRAM)
cleanoutput:;   rm -rf dustkap*
clean:;		rm -f $(OBJS) $(PROGRAM) *.mod *.i
		make cleanoutput
		rm -rf *~ \#* *.tex *.log auto notes.pdf notes.html optool.dSYM
install:	$(PROGRAM)
		mv $(PROGRAM) $(DEST)
test:; 		echo Computing size-integrated opacities ...
		make cleanoutput
		make
		./optool -s
		ipython -i optool.py
testdiv:;	echo computing size-dependant opacities ...
		make cleanoutput
		make
		./optool -d 3 -s
		ipython -i optool.py
quicktest:;	echo Computing size-integrated opacities ...
		make cleanoutput
		make
		./optool -na 10 -nl 30 -s -t
		ipython -i optool.py
quicktestdiv:;	echo computing size-dependant opacities ...
		make cleanoutput
		make
		./optool -na 10 -nl 30 -d 3 -s
		ipython -i optool.py

# how to compile program 
.SUFFIXES : .o .f .f90

.f.o:
	$(FC) $(LDFLAGS) -c $<

.f90.o:
	$(FC) $(LDFLAGS) -c $<

$(PROGRAM):     $(OBJS)
		$(LINKER) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)

ingest:
	./misc/ingestlnk.pl lnk_data/*.lnk > ref_ind.f90

