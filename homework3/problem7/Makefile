#CC/codes/fortran/newton/Makefile
OBJECTS = functions.o newton.o test.o
OBJECTS2 = functions.o newton.o testquartic.o
MODULES = functions.mod newton.mod

FFLAGS = -g

.PHONY: test clean testquartic

test: test.exe
		./test.exe

test.exe: $(MODULES) $(OBJECTS)
		gfortran $(FFLAGS) $(OBJECTS) -o test.exe

testquartic: testquartic.exe
		./testquartic.exe
#output.txt: intersections.exe
#		./intersections.exe > output.txt

testquartic.exe: $(MODULES) $(OBJECTS2)
		gfortran $(FFLAGS) $(OBJECTS2) -o testquartic.exe

%.o : %.f90
		gfortran $(FFLAGS) -c  $<

%.mod: %.f90
		gfortran $(FFLAGS) -c $<

clean:
		rm -f *.o *.exe *.mod
