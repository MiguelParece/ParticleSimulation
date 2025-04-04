CXX=g++
CXXFLAGS=-O2 -I. -fopenmp
LDFLAGS=-lm

all: parsim-omp

parsim-omp:
	$(CXX) $(CXXFLAGS) -o parsim-omp parsim-omp.cpp $(LDFLAGS)

clean:
	rm -f parsim-omp
