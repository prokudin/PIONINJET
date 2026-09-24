# Here we will create a library out of all object files I do not
# want to change.
# Under MAC
OS := $(shell uname)
#on Linux it gives Linux
ifeq ($(OS),Linux)
FC = gfortran
else
FC = gfortran
endif

FFLAGS = -g


PDF = aac03e.o grsv2000pdf_g1.o DSSV-2008.o


OBJALL = $(PDF)

sample1.exe: dssv-uncertainties-example.o $(PDF)
	$(FC) $(FFLAGS) -o $@ $+

sample.exe: sample.o $(PDF)
	$(FC) $(FFLAGS) -o $@ $+


clean:	
	rm *.o; rm sample.exe;

all:
	rm libPDFPOLARIZED.a ; ar rvs libPDFPOLARIZED.a $(OBJALL) ;
