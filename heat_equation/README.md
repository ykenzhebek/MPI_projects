# Parallel Implementation of Thomas Algorithm for 2D Heat Equation

This repository presents a parallel implementation of the Thomas algorithm for solving the two-dimensional heat equation using MPI. The numerical solution of the heat conduction problem is achieved through the two-step iteration process of the Alternating Direction Implicit (ADI) method. [paper](https://doi.org/10.26577/JMMCS-2019-3-24)

## Contents

- `2D_heat_sequential_program.cpp`: Sequential program for solving the two-dimensional heat equation.
- `2Dcoefficients1Dec.h`: Header file for the sweep method (forward step and backward) in the Thomas algorithm.
- `2DparamUnkn1Dec.h`: Header file for the sweep method in the second step of the Thomas algorithm.
- `heat_parallel_1D_decomposition.cpp`: Parallel program implementing the Thomas algorithm using 1D decomposition for the two-dimensional heat equation.
- `heat_parallel_block_decomposition.cpp`: Parallel program implementing the Thomas algorithm using 2D (block) decomposition for the two-dimensional heat equation.
- `run.sh`, `run1024.sh`, `run2048.sh`, `run4096.sh`, `run512.sh`: Shell scripts for running the programs with different matrix sizes.

## Parallelization Approach

The parallelization of the Thomas algorithm is challenging due to dependent data transfers. In this work, the Thomas algorithm was employed for parallelization using both 1D and 2D data decomposition. Particularly, in 2D data decomposition, the Thomas algorithm was applied along each x-axis and y-axis direction.

## Results and Analysis

The speedup and efficiency of parallel programs using 1D and 2D data decomposition are presented in the form of tables and graphs. The algorithms were tested on a cluster at the computing center of Novosibirsk State University for various matrix sizes ranging from 512x512 to 4096x4096. The obtained test results are analyzed to describe the features of the used decompositions.
