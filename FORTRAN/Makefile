# Here we will create a library out of all object files I do not
# want to change.

# Library of FF parametrizations updated to be posted on               
# http://www.pv.infn.it/~radici/FFdatabase/                            
# site maintained by Marco Radici (radici@pv.infn.it
FRAGMENT = ./frag/pkhff.o ./frag/kkp.o ./frag/polin2.o ./frag/polint.o \
	 ./frag/grille_had_charged.o ./frag/locate.o ./frag/dlib.o \
	 ./frag/fDSS.o	./frag/hknsff07.o  ./frag/akk.o

CERNLIB = ./cernlib/dadapt.o ./cernlib/dgs56p.o

PION = ./pion/GRVFGP.o ./pion/GRVFQBP.o ./pion/GRVFVP.o ./pion/grvpilo.f ./pion/mrsspi.o

OBJALL = $(FRAGMENT) $(CERNLIB) $(PION)

all:
	cd ./frag && make
	cd ./cernlib && make
	cd ./pion && make
	rm libMYLIB.a ; ar rvs libMYLIB.a $(OBJALL) ;

clean:
	cd ./frag && make clean
	cd ./cernlib && make clean
	cd ./pion && make clean
	rm libMYLIB.a
