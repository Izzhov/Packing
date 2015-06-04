CPPFLAGS = -I"../../PackingLib" -O3 -Wall -c -fmessage-length=0

libdir = ../../PackingLib/Cluster

lib = $(libdir)/libPackingLib.a

LDDIRS = -L$(libdir)

LDLIBS = -lPackingLib -lgsl -lgslcblas -lboost_program_options

all: 
	@echo "making all."
	cd $(libdir) && make
	make Sphere2DSteepest.out

Sphere2DSteepest.o: ../src/Sphere2DSteepest.cpp
	$(CXX) -c $(CPPFLAGS) $< -o $@

Sphere2DSteepest.out: Sphere2DSteepest.o $(lib)
	$(CXX) -O3 $(LDDIRS) -o $@ $< $(LDLIBS)

$(lib):%: $(libdir)/Makefile
	cd $(libdir) && make

clean:
	rm -f Sphere2DSteepest.out Sphere2DSteepest.o
	cd $(libdir) && make clean

clean_local:
	rm -f Sphere2DSteepest.out Sphere2DSteepest.o
