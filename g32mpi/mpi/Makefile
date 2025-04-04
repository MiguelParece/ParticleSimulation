CXX=mpic++
CXXFLAGS=-g -O2 -I. -fopenmp
LDFLAGS=-lm

all: parsim-omp

parsim-omp:
	$(CXX) $(CXXFLAGS) -o parsim-mpi parsim-mpi.cpp $(LDFLAGS)

clean:
	rm -f parsim-mpi