# Gdon3d
Three-dimensional Discontinuous Galerkin solver for linear conservation law systems.
### Summary
This solver implements a Discontinuous Galerkin method of order two on unsctructured 
3D tetrahedral meshes. To achieve the time integration, two algorithms are provided, one is an order 3 low-storage Runge-kutta scheme, the second an explicit CFL-less scheme (see: [HAL](https://hal.archives-ouvertes.fr/hal-03218086 "HAL")) derived from a relaxation method and a Crank-Nicolson.

### How to Build and Install
This project, written in C99, uses CMake > 3.13 and the libraries : `SuiteSparse`, `OpenMP` and `Hdf5`. On Ubuntu systems you can install those libs using : 
`sudo apt-get install libsuitesparse-dev libomp-dev libhdf5-dev`

To now build gdon3d on your host system, follow the following steps:
1. `git clone https://github.com/p-gerhard/dg-cfl-less-maxwell-solver.git` -- download the source
2. `mkdir build && cd build` -- create a build directory outside the source tree
3. `cmake ..` -- run CMake to setup the build
4. `make -j` -- compile the code

### Example
One example of how to use the solver is provided in `./examples/solve_example.c`. Before running it, please check that the mesh file (`.msh ASCII-v2`) exists in `./data/mesh/` directory.