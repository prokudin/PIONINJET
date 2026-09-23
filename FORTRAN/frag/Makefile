OS := $(shell uname)
#on Linux it gives Linux
ifeq ($(OS),Linux)
FC = gfortran
else
FC = gfortran 
endif

FFLAGS = -g

OBJECTS = fDSS.o pkhff.o kkp.o polin2.o polint.o grille_had_charged.o locate.o dlib.o akk.o hknsff07.o

toy.run: toy.o $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $+

clean:
	rm *.o; rm toy.run
