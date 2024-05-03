# MPI Projects

This repository contains implementations of parallel algorithms for solving the Poisson equation and the two-dimensional heat equation using MPI, OpenMP, and hybrid MPI+OpenMP approaches.
Also, repo contains some implemented MPI exercises.

## Contents

### Poisson Equation
- `poisson_mpi.c`: MPI implementation of the Poisson equation solver.
- `poisson_openmp.c`: OpenMP implementation of the Poisson equation solver.
- `poisson_hybrid.c`: Hybrid MPI+OpenMP implementation of the Poisson equation solver.
- `run_mpi.sh`: Bash script for running the MPI implementation.
- `run_openmp.sh`: Bash script for running the OpenMP implementation.
- `run_hybrid.sh`: Bash script for running the hybrid MPI+OpenMP implementation.
- `command.txt`: Information about the compilation commands for each implementation.

### Heat Equation
- `2D_heat_sequential_program.cpp`: Sequential program for solving the two-dimensional heat equation.
- `2Dcoefficients1Dec.h`: Header file for the sweep method (forward step and backward) in the Thomas algorithm.
- `2DparamUnkn1Dec.h`: Header file for the sweep method in the second step of the Thomas algorithm.
- `heat_parallel_1D_decomposition.cpp`: Parallel program implementing the Thomas algorithm using 1D decomposition.
- `heat_parallel_block_decomposition.cpp`: Parallel program implementing the Thomas algorithm using 2D (block) decomposition.
- `run.sh`, `run1024.sh`, `run2048.sh`, `run4096.sh`, `run512.sh`: Shell scripts for running the programs with different matrix sizes.

### MPI Practice Exercises
This section contains various MPI practice exercises for learning and understanding MPI concepts and functions.
- `io_mpi_file`: Implementation of MPI file I/O operations.
- `creating_topology.cpp`: Example code demonstrating the creation of custom MPI topologies.
- `ping-pong.cpp`: Implementation of the ping-pong communication pattern using MPI.
- `p-pong_sendrecv.cpp`: Implementation of the ping-pong communication pattern using MPI_Sendrecv.
- `ring.cpp`: Implementation of the ring communication pattern using MPI.
- `simpson.cpp`: Implementation of numerical integration using Simpson's rule with MPI parallelization.
- `trapezoidal_rule.cpp`: Implementation of numerical integration using the trapezoidal rule with MPI parallelization.

## Results and References

Results for the Poisson equation solver can be found in the [article](https://doi.org/10.26577/JMMCS-2018-3-523), while results for the two-dimensional heat equation solver are available in the [paper](https://doi.org/10.26577/JMMCS-2019-3-24). These papers showcase the speedup and efficiency of the parallel algorithms, along with the details of their implementation and analysis.
