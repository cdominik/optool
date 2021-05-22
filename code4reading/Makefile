

all:
	make readfitsfortran
	make readasciifortran

readfitsfortran: readfits.f90
	gfortran readfits.f90 -lcfitsio -L/usr/local/lib/ -o readfitsfortran

readasciifortran: readascii.f90
	gfortran readascii.f90 -o readasciifortran

clean:
	rm -rf readfitsfortran readasciifortran *~
