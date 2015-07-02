CPPFLAGS= -I"../../NumDiff" -O3 -Wall -c -fmessage-length=0 -lgsl -lgslcblas

objects = ConNet.o HarmPot.o mm.o MultiMD.o SingleMD.o Sizes.o SphDynMatHarmOsc.o Sphere.o SphereConNet.o Spherocyl.o SteepDesc.o TimeRNG01.o TimeRNGNormal.o Torus.o

all: libModPacking.a
	@echo "making all."

%.o: ../%.cpp ../%.h
	$(CXX) -c $(CPPFLAGS) $< -o $@

libModPacking.a: $(objects)
	ar -r $@ $(objects)

clean:
	rm -f libModPacking.a $(objects)
