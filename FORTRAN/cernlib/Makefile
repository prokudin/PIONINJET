OS := $(shell uname)
#on Linux it gives Linux
ifeq ($(OS),Linux)
FC = gfortran
else
FC = gfortran 
endif

FFLAGS = -g

OBJECTS = dadapt.o dgs56p.o


toy.run: toy.o $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $+

clean:
	rm *.o; rm toy.run
