CPPFLAGS=-O3 -Wall -c -fmessage-length=0 -lgsl -lgslcblas

objects = mymath.o NbrList2D.o Potential2DNoNbr.o Sphere2D.o SphereConNet2D.o SphereSizes.o Sphero2D.o SpheroPot2D.o SpheroSizes.o SpheroTorus2D.o TimeRNG01.o Torus2D.o

all: libPackingLib.a
	@echo "making all."

%.o: ../%.cpp ../%.h
	$(CXX) -c $(CPPFLAGS) $< -o $@

libPackingLib.a: $(objects)
	ar -r $@ $(objects)

clean:
	rm -f libPackingLib.a $(objects)
