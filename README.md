# propc

This program aims to be a modular toolkit for the calculation of physical properties from MD simulation outputs.

Tools:
Calculates Rg from average distance of all atoms.
Calculates Ree.
Calculates p(q) from the positions of atoms.
Rebuilds polymer, removing solvent and PBC.

## Compile

To compile use:

```
f95 -o propc -fopenmp -O3 propc
```

## Usage

To run, use:
```
./propc < parameters.in
```

## Parameters

Starting on line 2:
inputfile_name
outputfile_prefix
bool_calc_rg?
bool_calc_ree?
bool_calc_pq?
bool_output_trj?A
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

