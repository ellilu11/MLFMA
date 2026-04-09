## Multilevel Fast Multipole Method for Helmholtz Kernel

C++ implementation of the Multilevel Fast Multipole Method (MLFMM) for simulating Helmholtz-kernel mediated radiation and scattering phenomena.

### Algorithm

The general FMM algorithm accelerates computation of `N`-particle interactions represented by Laplace/Helmholtz-like kernels,
effecting a speedup from order `O(N^2)` to `O(N log N)`. It recursively divides the computational domain into successively 
finer boxes, aggregates source contributions via their multipole expansions, and transforms these into local expansions 
in each box from which solutions due to distant sources can be quickly evaluated.

For the Helmholtz kernel, the so-called "multipole" expansion is actually a modal plane wave expansion (Weyl expansion) along a set of angular samples. 
To account for the frequency dependence of these plane waves, the algorithm adjusts the sampling rate for the box level (with larger boxes sampling at a higher rate), 
hence the term "Multilevel" FMM. This code uses Lagrange interpolation to pass field coefficients between boxes at different levels while preserving their angular content.
A future enhancement involves accelerating this scheme using FFT-based techniques.

A typical Method of Moments application uses an iterative solver to deduce a surface unknown (e.g. electric or magnetic current) by directly or FMM-wise evaluating interactions 
between source and testing functions (e.g. point dipoles or RWG functions). This code uses the generalized minimal residual method (GMRES) with incomplete LU (ILU) preconditioner 
to solve for the surface currents.

### Build

CMake

##### Requirements

* C++ compiler with C++20 or newer
* Eigen library (included)

### Run

[TODO]