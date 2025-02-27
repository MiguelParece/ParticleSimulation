CXX=g++
CXXFLAGS=-O2 -I. -fopenmp
LDFLAGS=-lm

all: parsim

parsim: parsim.cpp
	$(CXX) $(CXXFLAGS) -o parsim parsim.cpp $(LDFLAGS)

clean:
	rm -f parsim