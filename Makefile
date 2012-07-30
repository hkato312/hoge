# testOctave.c

program = testOctave
objs = testOctave.o

CFLAGS = -g -Wall -I/usr/include/octave-3.6.1
LDFLAGS = -L/usr/lib/octave/3.6.1
LIBS = -loctave

$(program):	$(objs)
	$(CXX) -o $(program) $^ $(LDFLAGS) $(LIBS)

testOctave.o:	testOctave.cpp
	$(CXX) $(CFLAGS) -c $<

clean:
	$(RM) $(program) $(objs)
