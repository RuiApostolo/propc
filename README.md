# propc

This program aims to be a modular toolkit for the calculation of physical properties from MD simulation outputs.

Tools:
- Calculates Rg from average distance of all atoms (in the same molecule, or all).
- Calculates Ree, if ReeFirstAtom and ReeLastAtom are 0, they will be taken as the first and last atom. Any other number (between 1 and MolSize) will be taken as the index of the atom to be used, starting with the first atom of MolStartType as index 1, and going in the order of the datafile.
- Calculates p(q) from the positions of atoms for one or more molecules.
- Calculates p(q) from the positions of atoms for smae-molecule and different-molecule in multi-molecule systems.
- Rebuilds polymer, removing solvent and PBC.

## Compile

To compile use:

```bash
gfortran propc.f90 -o propc -fopenmp -O3
```

## Usage

To run, use:
```bash
./propc < parameters.in
```

## Parameters

Starting on line 2:

```
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
Index of the first atom for the Ree calculation (0 is the same as index 1, i.e., the first atom of type MolStartType)
Index of the last atom for the Ree calculation (0 is the same as index MolSize)
Output information every # timesteps  
Exponent of lower limit of p(q) calculation  
Exponent of higher limit of p(q) calculation  
Number of points used as q in p(q) calculation  
```

## Changelog

* Version 2.2.1 -- 2021/04/05
  * fixed timing module - ETA wasn't considering skipped timesteps into calculation, so it was underestimating total time in the cases where timesteps were being skipped.

* Version 2.2 -- 2021/03/30:
  * fixed timing modules - ETA had a bug.
  * added ReeFirstAtom and ReeLastAtom to allow calculation of Ree in molcules where the edge atoms where not the first and last in the datafile order.
  * added version number to params.in

* Version 2.1 -- 2021/02/03:
  * Added 'Indivudual Rg' calculation switch: if true calculates Rg for each molecule, if false consideres every molecule part of the same 'cluste' and calculates the Rg for said cluster (new 'Total Rg' subroutine).
  * Added OMP support for the new 'Total Rg' subroutine.
  * Fixed timing modules.
