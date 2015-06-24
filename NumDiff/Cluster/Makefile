CPPFLAGS= -O3 -Wall -c -fmessage-length=0

objecks = Const.o Div.o DPow.o Fcon.o Function.o IntPow.o Minus.o Plus.o Times.o Variable.o

all: libNumDiff.a
	@echo "making all."

%.o: ../%.cpp ../%.h
	$(CXX) -c $(CPPFLAGS) $< -o $@

libNumDiff.a: $(objecks)
	ar -r $@ $(objecks)

clean:
	rm -f liNumDiff.a $(objecks)
