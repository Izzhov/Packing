CPPFLAGS = -I"../../NumDiff" -I"../../ModPacking" -O3 -Wall -c -fmessage-length=0

libdir = ../../ModPacking/Cluster
libnumdir = ../../NumDiff/Cluster

lib = $(libdir)/libModPacking.a
libnum = $(libnumdir)/libNumDiff.a

LDDIRS = -L$(libdir) -L$(libnumdir)

LDLIBS = -lNumDiff -lModPacking -lgsl -lgslcblas -lboost_program_options

all: 
	@echo "making all."
	cd $(libnumdir) && make
	cd $(libdir) && make
	make ModSph2DSteep.out

ModSph2DSteep.o: ../src/ModSph2DSteep.cpp
	$(CXX) -c $(CPPFLAGS) $< -o $@

ModSph2DSteep.out: ModSph2DSteep.o $(lib) $(libnum)
	$(CXX) -O3 $(LDDIRS) -o $@ $< $(LDLIBS)

$(lib):%: $(libdir)/Makefile
	cd $(libdir) && make

$(libnum):%: $(libnumdir)/Makefile
	cd $(libnumdir) && make

clean:
	rm -f ModSph2DSteep.out ModSph2DSteep.o
	cd $(libdir) && make clean
	cd $(libnumdir) && make clean

clean_local:
	rm -f ModSph2DSteep.out ModSph2DSteep.o
