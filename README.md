# ecoSLM

This repository houses software and datasets (experimental and simulation) on the ecoSLM - a fast and approximate slip link model for predicting the linear viscoelasticity of blends of entangled star and linear polymers.

## Code

Code files contain three Fortran 90 programs, which can be compiled as follows:

gfortran -O4 Utils1.f90 Initialize1.f90 -o ini
gfortran -O4 Utils1.f90 Dynamics1.f90 -o dyna

The executable file `ini` can be used to initialize an ensemble, while `dyna` can be used to carry out dynamics.

The input file `inp.dat` is used to specify the simulation ensemble for a binary blend (star/linear, linear/linear, star/star).

> Is fraction linear? [F if fraction is linear, T if it is linear]
> F T
> Number of chains, Np1, Np2 [# chains of component 1 and component 2;]
> 2950 50
> Number of SL, Z1, Z2
> 8.0 116.1
> Random number seed
> -3
> ConstraintRelease
> T

