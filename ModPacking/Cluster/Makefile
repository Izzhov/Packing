CPPFLAGS=-O3 -Wall -c -fmessage-length=0 -lgsl -lgslcblas

objects = HarmPot.o mm.o Sizes.o Sphere.o SteepDesc.o TimeRNG01.o Torus.o

all: libModPacking.a
	@echo "making all."

%.o: ../%.cpp ../%.h
	$(CXX) -c $(CPPFLAGS) $< -o $@

libModPacking.a: $(objects)
	ar -r $@ $(objects)

clean:
	rm -f libModPacking.a $(objects)
