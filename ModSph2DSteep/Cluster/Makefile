CPPFLAGS = -I"../../ModPacking" -O3 -Wall -c -fmessage-length=0

libdir = ../../ModPacking/Cluster

lib = $(libdir)/libModPacking.a

LDDIRS = -L$(libdir)

LDLIBS = -lModPacking -lgsl -lgslcblas -lboost_program_options

all: 
	@echo "making all."
	cd $(libdir) && make
	make ModSph2DSteep.out

ModSph2DSteep.o: ../src/ModSph2DSteep.cpp
	$(CXX) -c $(CPPFLAGS) $< -o $@

ModSph2DSteep.out: ModSph2DSteep.o $(lib)
	$(CXX) -O3 $(LDDIRS) -o $@ $< $(LDLIBS)

$(lib):%: $(libdir)/Makefile
	cd $(libdir) && make

clean:
	rm -f ModSph2DSteep.out ModSph2DSteep.o
	cd $(libdir) && make clean

clean_local:
	rm -f ModSph2DSteep.out ModSph2DSteep.o
