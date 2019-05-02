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
f95 csv_file.f90 probc.f08 -fopenmp -o ./probc
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
bool_output_trj?