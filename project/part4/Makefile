
OBJECTS = random_util.o fun_u.o mc_walk.o Laplace.o
FFLAGS = 
NUM_PROCS ?= 4
.PHONY: Laplace plot clean clobber

Laplace: Laplace.exe
	mpiexec -n $(NUM_PROCS) Laplace.exe

Laplace.exe: $(OBJECTS)
	mpif90 $(FFLAGS) $(OBJECTS) -o Laplace.exe


%.o : %.f90
	mpif90 $(FFLAGS) -c  $< 

mc_Laplace_error.txt: Laplace.exe
	mpiexec -n $(NUM_PROCS) Laplace.exe

plot: mc_Laplace_error.txt
	python plot_mc_laplace_error.py

clean:
	rm -f *.o *.exe *.mod

clobber: clean
	rm -f mc_laplace_error.txt mc_laplace_error.png
