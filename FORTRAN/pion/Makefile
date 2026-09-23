OS := $(shell uname)
#on Linux it gives Linux
ifeq ($(OS),Linux)
FC = gfortran
else
FC = gfortran 
endif

FFLAGS = -g

OBJECTS = grvpilo.o GRVFGP.o GRVFQBP.o GRVFVP.o mrsspi.o

toy.run: toy.o $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $+

clean:
	rm *.o; rm toy.run
