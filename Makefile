# makefile for OpTool
#

FC	  = gfortran
LINKER	  = gfortran

ifeq ($(ifort),true)
  FC	  = ifort
  LINKER  = ifort
endif

# array checks for debugging
ifeq ($(debug),true)
  ifeq ($(ifort),true)
    DEBUGGING = -check all -traceback -check bounds -O0 -g -check -fpe1
  else	
    DEBUGGING = -fbounds-check -fbacktrace
  endif
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

ifeq ($(fitsio),true)
  FLAG_FITS		= -DUSE_FITSIO
  LIBS_FITS		= -lcfitsio -L/usr/local/lib/
endif

FFLAGS  = $(FLAG_ALL) $(FLAG_LINUX) $(FLAG_FITS)
LDFLAGS = $(FLAG_ALL) $(FLAG_LINUX) $(FLAG_FITS)
LIBS    = $(LIBS_FITS)

# Object files to link
OBJS	= optool.o \
	  optool_extra.o \
	  ref_ind.o \
	  dmilay_f95.o \
          meanopacity.o \


# program name and install location
PROGRAM       = optool
DEST	      = ${HOME}/bin

# make actions 
all:		$(PROGRAM)
cleanoutput:;   rm -rf dustkap*
clean:;		rm -f $(OBJS) $(PROGRAM) *.mod *.i
		make cleanoutput
		rm -rf *~ \#* *.tex *.log auto notes.pdf notes.html 
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
		./optool -na 10 -nl 30 -s
		ipython -i optool.py
quicktestdiv:;	echo computing size-dependant opacities ...
		make cleanoutput
		make
		./optool -na 10 -nl 30 -d 3 -s
		ipython -i optool.py

# how to compile program 
.SUFFIXES : .o .f .f90 .F

.f.o:
	$(FC) $(LDFLAGS) -c $<

.f90.o:
	$(FC) $(LDFLAGS) -c $<

.F.o:
	$(FC) $(LDFLAGS) -c $<

$(PROGRAM):     $(OBJS)
		$(LINKER) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)

ingest:
	./misc/ingestlnk.pl lnk_data/*.lnk > ref_ind.f90

README.org: notes.org
	perl -ne 'if (/^\* OpTool User Guide/../^\* Materials/) {s/@@.*?@@//g; print unless /^\* Materials/}' notes.org > README.org
