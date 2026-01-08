# Installation and Compilation Instructions

## Prerequisites
- A Fortran compiler (e.g., `gfortran`).
- The basis set file `6-31g.1.dalton` must be present in the directory where you run the executable.

## Compilation

To compile the different exercises, navigate to the main directory and run the following commands:

### Exercise 3.2 (Overlap Norms)
```bash
gfortran src/basis_module.f90 src/math_module.f90 src/main_3_2.f90 -o ex3_2
./ex3_2

#### Exercise 3.3 (MacMurchie-DavidSon Diatomic Overlap)
```bash
gfortran src/basis_module.f90 src/mmd_module.f90 src/main_3_3.f90 -o ex3_3
./ex3_3
