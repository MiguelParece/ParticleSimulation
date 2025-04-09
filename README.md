# MPI Parallel N-Body Particle Simulation

A distributed-memory parallel N-body simulation modeling gravitational interactions between particles with periodic boundary conditions, implemented using MPI for high-performance computing.

## Features

- **MPI Parallelization**: Distributed spatial decomposition across MPI processes.
- **Cell-Based Spatial Partitioning**: Reduces force calculation complexity via grid-based particle grouping.
- **Periodic Boundary Conditions**: Particles wrap around domain edges (vertical/horizontal).
- **Collision Detection**: Particles merge when within `EPSILON` distance.
- **Performance Optimizations**: 
  - Non-blocking MPI communication.
  - Ghost cell exchange for boundary handling.
  - Hybrid OpenMP+MPI support.

## Dependencies

- MPI (OpenMPI/Intel MPI)
- C++17 compiler
- (Optional) `libunwind` and `binutils` for profiling

## Build

```bash
mpic++ -O3 -fopenmp -std=c++17 main.cpp -o particle_sim