all : exec

# Options de débogage et de compilation
# DEBUG = -Wall -fcheck=all -ffpe-trap=invalid -pedantic -O0

# Options de débogage et de compilation
DEBUG = -Wall -fcheck=all -ffpe-trap=invalid -pedantic -O0 
        # -finit-integer=zero # -finit-real=snan 


# Utilisation de mpif90 pour gérer MPI correctement
FC = gfortran
FCFLAGS = $(DEBUG)

exec : const.o visualisation.o functions.o main.o
	$(FC) $(FCFLAGS) -o exec const.o visualisation.o functions.o main.o

functions.o : functions.f90 const.o visualisation.o
	$(FC) $(FCFLAGS) -c functions.f90

visualisation.o : visualisation.f90 const.o
	$(FC) $(FCFLAGS) -c visualisation.f90

const.o : const.f90
	$(FC) $(FCFLAGS) -c const.f90

# algebre.o : const.o algebre.f90
# 	$(FC) $(FCFLAGS) -c algebre.f90

# solveur.o : const.o algebre.o solveur.f90
# 	$(FC) $(FCFLAGS) -c solveur.f90

# fonc.o : const.o algebre.o solveur.o fonc.f90
# 	$(FC) $(FCFLAGS) -c fonc.f90

# sch_temps.o : const.o algebre.o solveur.o fonc.o sch_temps.f90
# 	$(FC) $(FCFLAGS) -c sch_temps.f90

main.o : const.o visualisation.o functions.o main.f90
	$(FC) $(FCFLAGS) -c main.f90

clean :
	rm -f *.o *.mod exec
	rm -f *.dat