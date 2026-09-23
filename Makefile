#OS := $(shell uname)
OS := $(shell whoami)
#on Linux it gives Linux
ifeq ($(OS),prokudin) # This is what Tanya's macbook gives
CPP = g++
FF = gfortran
else # Mac
CPP = g++
FF = gfortran
endif


# Root is the same for both systems...
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs) -lMinuit
ROOTGLIBS := $(shell root-config --glibs)

LDFLAGS = -O
CFLAGS += $(ROOTCFLAGS)

#LHAPDF
LHAPDFLIB  =  -L/usr/local/lib -lLHAPDF
#header can be found in
#ls /usr/local/include/LHAPDF/
LHAPDFINCLUDE := $(shell lhapdf-config --incdir)
# for BOOST libraries -I/usr/local/Cellar/boost/1.59.0/include/

ifeq ($(OS),prokudin) # This is what Tanya's macbook gives
OPT = -c -O3 -arch x86_64 -m64 -Wall -I./ -I./hoppet-1.1.5/src/ -I$(LHAPDFINCLUDE) \
	-I$(LHAPDFINCLUDE)/LHAPDF -I/usr/local/Cellar/boost/1.60.0_2/include/ -I/Users/prokudin/ScienceSoftware/Cuba-4.2/
else # Mac
OPT = -c -O3 -arch x86_64 -m64 -Wall -I./ -I./hoppet-1.1.5/src/ -I$(LHAPDFINCLUDE) \
	-I$(LHAPDFINCLUDE)/LHAPDF -I/usr/local/Cellar/boost/1.59.0/include/ -I/Users/avp5627/ScienceSoftware/Cuba-4.2 \
	-I/usr/local/Cellar/boost/1.64.0_1/include/
endif

# Hoppet implements DGLAP
HOPPETLIB =  -L./hoppet-1.1.5/src  -lhoppet_v1


# In order for gfortran to work I need to include this
FLIBS =  -L/usr/local/gfortran/lib/ -lgfortran

# My fortran lib
FORTRAN_LIB = -L./FORTRAN -lMYLIB

#PDF_POLARISED
PDF_POLARIZED_LIB = -L./pdf_polarized   -lPDFPOLARIZED

#CUBA_LIB
ifeq ($(OS),prokudin) # This is what Tanya's macbook gives
CUBA_LIB = -L/Users/prokudin/ScienceSoftware/Cuba-4.2/  -lcuba
else # Mac
CUBA_LIB = -L/Users/avp5627/ScienceSoftware/Cuba-4.2  -lcuba
endif

# All libraries:
LIB = $(ROOTLIBS) $(LHAPDFLIB) $(CUBA_LIB) \
	$(PDF_POLARIZED_LIB)  $(FORTRAN_LIB) $(FLIBS) $(HOPPETLIB)



# This is transversity coming from Werner Vogelsang
TRANS_LO = ./transv_pdf.o

OBJECTS_UNIVERSAL = particle.o hadron.o lepton.o parton.o parameters.o $(TRANS_LO) \
	evolution_common.o experimental_data.o gauss_2d_sphere.o gauss_legendre.o

OBJECTS_SIDIS = sidis.o cuba_sidis.o \
	setexperiment_sidis.o \
	reader_sidis.o \
	draw_result_sidis.o \
	twist3.o \
	approximate_collins_pp.o \
	pionjet_pp.o \
	approximate_sivers_pp.o




#drawAN
drawan_pp.exe: $(OBJECTS_UNIVERSAL) $(OBJECTS_SIDIS) drawan_pp.o
	$(CPP)  -o drawan_pp.exe drawan_pp.o $(OBJECTS_UNIVERSAL)  \
	$(OBJECTS_SIDIS) \
	$(LIB)

#draw pion jet
drawcollins_pionjet_pp.exe: $(OBJECTS_UNIVERSAL) $(OBJECTS_SIDIS) drawcollins_pionjet_pp.o
	$(CPP)  -o drawcollins_pionjet_pp.exe drawcollins_pionjet_pp.o $(OBJECTS_UNIVERSAL)  \
	$(OBJECTS_SIDIS) \
	$(LIB)

#drawsivers JET!
drawsivers_pp_jet_sidis.exe: $(OBJECTS_UNIVERSAL) $(OBJECTS_SIDIS) drawsivers_pp_jet_sidis.o
		$(CPP)  -o drawsivers_pp_jet_sidis.exe drawsivers_pp_jet_sidis.o $(OBJECTS_UNIVERSAL)  \
		$(OBJECTS_SIDIS) \
		$(LIB)

#drawsivers PHOTON!
#drawsivers_pp_photon.exe: $(OBJECTS_UNIVERSAL) $(OBJECTS_SIDIS) drawsivers_pp_photon.o
#		$(CPP)  -o drawsivers_pp_photon.exe drawsivers_pp_photon.o $(OBJECTS_UNIVERSAL)  \
#		$(OBJECTS_SIDIS) \
#		$(LIB)



#some generic definitions BEGIN
particle.o:	particle.cpp
	$(CPP) $(OPT) $(CFLAGS) -o $@ particle.cpp
	@echo "..................done Particle."

hadron.o:	hadron.cpp
	$(CPP) $(OPT) $(CFLAGS) -o $@ hadron.cpp
	@echo "..................done Hadron."

lepton.o:	lepton.cpp
	$(CPP) $(OPT) $(CFLAGS) -o $@ lepton.cpp
	@echo "..................done Lepton."

parton.o:	parton.cpp
	$(CPP) $(OPT) $(CFLAGS) -o $@ parton.cpp
	@echo "..................done Parton."

parameters.o:	parameters.cpp
	$(CPP) $(OPT) $(CFLAGS) -o $@ parameters.cpp
	@echo "..................done Parameters."


experimental_data.o:	experimental_data.cpp
	$(CPP) $(OPT) $(CFLAGS) -o $@ experimental_data.cpp
	@echo "..................done Experimental_data."

#experiments are set here
setexperiment_sidis.o:	setexperiment_sidis.cpp
	$(CPP) $(OPT) $(CFLAGS) -o $@ setexperiment_sidis.cpp
	@echo "..................done Set Experiment SIDIS."


#reader for ROOT
reader_sidis.o:	reader_sidis.cpp
	$(CPP) $(OPT) $(CFLAGS) -o $@ reader_sidis.cpp
	@echo "..................done Reader SIDIS."

#transversity
transv_pdf.o: transv_pdf.f
	$(FF) -c transv_pdf.f
#some generic definitions END


#this is the MAIN file almost all definitions are here!
sidis.o:	sidis.cpp sidis.h
	$(CPP) $(OPT) $(CFLAGS) -o $@ sidis.cpp
	@echo "..................done SIDIS"

# numerical integration over 2d disc x^2 + y^2 < R^2, Pavel Holoborodko:
gauss_2d_sphere.o:	gauss_2d_sphere.c gauss_2d_sphere.h
	$(CPP) $(OPT) $(CFLAGS) -o $@ gauss_2d_sphere.c
	@echo "..................done gauss_2d_sphere"

#this one is used for twist-3 calculations
gauss_legendre.o: gauss_legendre.c
	$(CPP) $(OPT) $(CFLAGS) -o $@  gauss_legendre.c
	@echo "..................done  Gaus Legendre."


#main subroutine for integration
cuba_sidis.o:	cuba_sidis.cpp cuba_sidis.h
	$(CPP) $(OPT) $(CFLAGS) -o $@ cuba_sidis.cpp
	@echo "..................done Cuba SIDIS."


# evolution common functions
evolution_common.o:	evolution_common.cpp
	$(CPP) $(OPT) $(CFLAGS) -o $@ evolution_common.cpp
	@echo "..................done Evolution Common."



#various calculations are here
approximate_collins_pp.o:	approximate_collins_pp.cpp
	$(CPP) $(OPT) $(CFLAGS) -o $@ approximate_collins_pp.cpp
	@echo "..................done Approximate Collins PP."

approximate_sivers_pp.o:	approximate_sivers_pp.cpp
		$(CPP) $(OPT) $(CFLAGS) -o $@ approximate_sivers_pp.cpp
		@echo "..................done Approximate Sivers PP."

pionjet_pp.o:	pionjet_pp.cpp
	$(CPP) $(OPT) $(CFLAGS) -o $@ pionjet_pp.cpp
	@echo "..................done Pionjet PP."


#twist-3
twist3.o:	twist3.cpp
	$(CPP) $(OPT) $(CFLAGS) -o $@ twist3.cpp
	@echo "..................done Twist3."

#Files for Drawing are here: BEGIN
draw_result_sidis.o:	draw_result_sidis.cpp
	$(CPP) $(OPT) $(CFLAGS) -o $@ draw_result_sidis.cpp
	@echo "..................done Draw Result SIDIS."

draw_sidis.o:	draw_sidis.cpp
	$(CPP) $(OPT) $(CFLAGS) -o $@ draw_sidis.cpp
	@echo "..................done Draw SIDIS."



#drawan_pp object
drawan_pp.o:	drawan_pp.cpp
	$(CPP) $(OPT) $(CFLAGS) -o $@ drawan_pp.cpp
	@echo "..................done Draw AN PP."


#drawan_pion jet object
drawcollins_pionjet_pp.o:	drawcollins_pionjet_pp.cpp
	$(CPP) $(OPT) $(CFLAGS) -o $@ drawcollins_pionjet_pp.cpp
	@echo "..................done Draw drawcollins_pionjet_pp."


drawsivers_pp_jet_sidis.o:	drawsivers_pp_jet_sidis.cpp 
		$(CPP) $(OPT) $(CFLAGS) -o $@ drawsivers_pp_jet_sidis.cpp
		@echo "..................done Draw Sivers PP JET SIDIS."


clean:
	rm *.o
