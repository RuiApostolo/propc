# propc

This program aims to be a modular toolkit for the calculation of physical properties from MD simulation outputs.

Tools:
- Calculates Rg from average distance of all atoms (in the same molecule, or all).
- Calculates Ree.
- Calculates p(q) from the positions of atoms for one or more molecules.
- Calculates p(q) from the positions of atoms for smae-molecule and different-molecule in multi-molecule systems.
- Rebuilds polymer, removing solvent and PBC.

## Compile

To compile use:

```
gfortran propc.f90 -o propc -fopenmp -O3
```

## Usage

To run, use:
```
./propc < parameters.in
```

## Parameters

Starting on line 2:


inputfile_name  
outputfile_prefix (must be string)  
bool_calc_rg?  
bool_calc_ree?  
bool_calc_pq?  
bool_calc_same/diff_pq? (for systems with number of molecules > 1, requires regular pq)  
bool_output_trj?  
Number of columns in inputfile (lammpstrj)  
Total number of timesteps  
Ignore first # timesteps  
Number of target molecules  
Size of target molecule  
Atomtype of first atom in target molecule  
Output information every # timesteps  
Exponent of lower limit of p(q) calculation  
Exponent of higher limit of p(q) calculation  
Number of points used as q in p(q) calculation  

## Changelog

* Version 2.1:
  * Added 'Indivudual Rg' calculation switch: if true calculates Rg for each molecule, if false consideres every molecule part of the same 'cluste' and calculates the Rg for said cluster (new 'Total Rg' subroutine).
  * Added OMP support for the new 'Total Rg' subroutine.
  * Fixed timing modules
