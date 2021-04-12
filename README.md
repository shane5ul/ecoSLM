# ecoSLM

This repository houses software and datasets (experimental and simulation) on the ecoSLM - a fast and approximate slip link model for predicting the linear viscoelasticity of blends of entangled star and linear polymers.

## References

1. "Fast Slip Link Model for Bidisperse Linear Polymer Melts", *Macromolecules*, **2019**, 52, 3092−3103 [doi: 10.1021/acs.macromol.8b02367].
2. "Mathematical foundations of an ultra coarse-grained slip link model", *J. Chem. Phys.*, **2019**, 151, 044903 [doi: 10.1063/1.5111032].
3. "Temporal Coarse-Graining in a Slip Link Model for Polydisperse Polymer Melts", *Front. Phys.*, **2020**, 8, 579499 [doi: 10.3389/fphy.2020.579499].

## Code

Code files contain three Fortran 90 programs, which can be compiled as follows:

gfortran -O4 Utils1.f90 Initialize1.f90 -o ini
gfortran -O4 Utils1.f90 Dynamics1.f90 -o dyna

The executable file `ini` can be used to initialize an ensemble, while `dyna` can be used to carry out dynamics.

The input file `inp.dat` is used to specify the simulation ensemble for a binary blend (star/linear, linear/linear, star/star).
<pre><code>
Is fraction linear? [F if component is star, T if it is linear] 
F T
Number of chains, Np1, Np2 [# chains of component 1 and component 2; must be even, and Np1 cannot be zero]
2950 50
Number of SL, Z1, Z2 [avg #slip links on component 1 and 2]
8.0 116.1
Random number seed [random number seed for multiple replicas]
-3
ConstraintRelease [mark F if you want to suppress CR]
T
</code></pre>

`ini` reads the input file, creates an ensemble, and equilibrates it by writing a snapshot file `inSnap.dat`.

`dyna` uses `inp.dat` and `inSnap.dat` and writes the final snapshot `outSnap.dat` and a file `relax.dat`. The latter file has three columns: the first column is the time, the second is the normalized stress and the last column is the normalized dielectric relaxation.

**Expert Note**: 

In Dynamics1.f90, calibration for stars and linears is currently hard-coded for PBd at 25C. The following lines in the code may be suitably modified for any other situation.

<pre><code>
gammaLin  = 1.4				
gammaStr  = 3.5

tau0Lin   = 0.12 
tau0Str   = 1.8e-3
</code></pre>

