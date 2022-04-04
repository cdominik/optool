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
	MULTICORE = -qopenmp -fp-model strict
    else
	MULTICORE = -fopenmp
    endif
endif

# Debugging flags
ifeq ($(debug),true)
  ifeq ($(ifort),true)
    DEBUGGING = -check all -traceback -check bounds -O0 -g -check -fpe1
  else	
    DEBUGGING = -fbounds-check -fbacktrace -O0 -ffpe-trap=invalid,zero,overflow -fcheck=all
  endif
endif

# CFITSIO support
ifeq ($(fits),true)
  FLAG_FITS		= -DUSE_FITSIO
  LIBS_FITS		= -lcfitsio -L/usr/local/lib/ -L.
endif
ifeq ($(fitsio),true)
  FLAG_FITS		= -DUSE_FITSIO
  LIBS_FITS		= -lcfitsio -L/usr/local/lib/ -L.
endif

# Platform specific compilation options
ifeq ($(ifort),true)
  FLAG_ALL      = -O3 -g -extend-source -zero -prec-div $(DEBUGGING) $(MULTICORE)
  FLAG_ALL_NM   = -O3 -g -extend-source -zero -prec-div $(DEBUGGING)
  FLAG_LINUX    = -xHOST -fpp
  FLAG_MAC      = -xHOST -opt-prefetch -static-intel -fpp
else
  FLAG_ALL      = -O3 -g -fdefault-double-8 -fdefault-real-8 $(DEBUGGING) $(MULTICORE)
  FLAG_ALL_NM   = -O3 -g -fdefault-double-8 -fdefault-real-8 $(DEBUGGING)
  FLAG_LINUX    = -cpp
  FLAG_MAC      = -m64 -cpp
endif

# Combine the flags
FFLAGS  = $(FLAG_ALL) $(FLAG_LINUX) $(FLAG_FITS)
LDFLAGS = $(FLAG_ALL) $(FLAG_LINUX) $(FLAG_FITS)
LIBS    = $(LIBS_FITS)

# Object files to link
OBJS	= optool.o optool_guts.o optool_manual.o optool_refind.o optool_fractal.o optool_geofractal.o

# Program name and install location
PROGRAM       = optool
bindir	      = ${HOME}/bin
BINRELEASE    = ~/Dropbox/Websites/uva.nl/WWW/optool


# make actions 
all:		$(PROGRAM)
install:	$(PROGRAM)
		cp optool optool2tex optool-complete $(bindir)

full:;		make multi=true fits=true
it:;		make clean
		make full
		make clean1
cleanoutput:;   rm -rf dustkap*.dat dustkap*.inp blended.lnk optool_sd.dat optool_lam.dat optool_tmp_output_dir_*
cleanbin:;	rm -f bin/optool* bin.zip
cleanlatex:;	rm -rf *.tex *.aux *.log *.dvi *.blg *.bbl auto optool.pdf
cleanpython:;	rm -rf optool.dSYM tmp.py __pycache__
cclean:;	rm -f $(OBJS) $(PROGRAM)
clean:;		make clean1
		rm -rf $(PROGRAM) optool.egg-info
clean1:;	rm -f $(OBJS) *.mod *.i *.html
		make cleanoutput cleanbin cleanlatex cleanpython
		rm -rf *~ \#*  selftest_optool

manual:;        /Applications/Emacs.app/Contents/MacOS/Emacs UserGuide.org \
		     --batch -f org-ascii-export-to-ascii --kill
		maint/bake_manual.pl > optool_manual.f90
		rm UserGuide.txt
pdf:;		/Applications/Emacs.app/Contents/MacOS/Emacs -l maint/bake_manual.el \
		     UserGuide.org --batch -f org-latex-export-to-pdf --kill
ingest:;	echo Compiling in datasets in lnk_data...
		./maint/ingestlnk.pl lnk_data/*.lnk > optool_refind.f90

release:;	make clean
		make pdf
		make manual
		make ingest
		git add *
		git commit -m "Release $(version)"
		git tag release_$(version)
		git push
		git push origin release_$(version)
		make binmac
		make binmv

binmac:;	make cclean
		make
		mv optool bin/optool-mac
		make cclean
		make multi=true
		mv optool bin/optool-mac-OpenMP
		make cclean
		make fits=true
		mv optool bin/optool-mac-fits
		make cclean
		make multi=true fits=true
		mv optool bin/optool-mac-fits-OpenMP
binlinux:;	make cclean
		make
		mv optool bin/optool-linux
		make cclean
		make multi=true
		mv optool bin/optool-linux-OpenMP
		make cclean
		make fits=true
		mv optool bin/optool-linux-fits
		make cclean
		make multi=true fits=true
		mv optool bin/optool-linux-fits-OpenMP
binzip:;	rm -f bin.zip
		(cd bin; zip ../bin.zip optool*)
binmv:;		mv bin/optool* ~/Dropbox/Websites/uva.nl/WWW/optool/

.SUFFIXES : .o .f .f90

.f.o:
	$(FC) $(LDFLAGS) -c $<

.f90.o:
	$(FC) $(LDFLAGS) -c $<

$(PROGRAM):     $(OBJS)
		$(LINKER) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)


