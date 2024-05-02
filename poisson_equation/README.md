# MPI, OpenMP, and Hybrid Parallel Algorithm for the Poisson Equation

This repository contains implementations of a hybrid parallel algorithm for solving the Dirichlet problem for the two-dimensional Poisson equation using MPI, OpenMP, and a hybrid combination of both.

## Contents

- `poisson_mpi.c`: Implementation of the Poisson equation solver using MPI.
- `poisson_openmp.c`: Implementation of the Poisson equation solver using OpenMP.
- `poisson_hybrid.c`: Implementation of the Poisson equation solver using a hybrid MPI+OpenMP approach.
- `run_mpi.sh`: Bash script for running the MPI implementation.
- `run_openmp.sh`: Bash script for running the OpenMP implementation.
- `run_hybrid.sh`: Bash script for running the hybrid MPI+OpenMP implementation.
- `command.txt`: Information about the compilation commands for each implementation.

## Running the Codes

For each implementation, refer to the respective `run_*.sh` script for execution on your system. Ensure that you have the necessary MPI and OpenMP compilers installed.

## Cluster Execution

The scripts provided are configured for execution on a cluster. Modify the PBS directives in the `.sh` files according to your cluster's requirements.

## Results

Tables and graphs showcasing the acceleration and efficiency of each parallel algorithm are presented in the [article](https://doi.org/10.26577/JMMCS-2018-3-523). The hybrid MPI+OpenMP algorithm demonstrates notable improvements in solving such problems, achieving a speedup of 1.5-2 times.
