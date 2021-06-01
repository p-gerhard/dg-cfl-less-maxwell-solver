# gdon3d
Three-dimensional Discontinuous Galerkin solver for linear conservation law systems.
### Summary
This solver implements a Discontinuous Galerkin method of order two on unsctructured 
3D tetrahedral meshes. To achieve the time integration, two algorithms are provided, one is an order 3 low-storage Runge-kutta scheme, the second an explicit CFL-less scheme (see: [hal](https://hal.archives-ouvertes.fr/hal-03218086 "hal")) derived from relaxation method and a Crank-Nicolson.


