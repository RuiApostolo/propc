# propc

This program aims to be a modular toolkit for the calculation of physical properties from MD simulation outputs.

Tools:
- Calculates Rg from average distance of all atoms.
- Calculates Ree.
- Calculates p(q) from the positions of atoms.
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

```
inputfile_name
outputfile_prefix (must be string)
bool_calc_rg?
bool_calc_ree?
bool_calc_pq?
bool_calc_individual_pq? (for systems with number of molecules > 1, requires regular pq)
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
```


before array change
real    57m0.304s
user    357m27.428s
sys     0m11.928s


after array change
real    58m15.939s
user    354m20.100s
sys     0m13.420s


OMP dynamic
real    55m17.439s
user    365m4.392s
sys     0m15.424s


FF loop changes
real    52m16.084s
user    366m29.540s
sys     0m14.940s


with ind_pq - wrong
real    12m5.175s
user    71m49.188s
sys     0m1.204s

with ind_pq - fixed?
real    42m22.843s
user    197m6.848s
sys     0m5.020s

added real(x,dp)
real    43m1.547s
user    199m14.588s
sys     0m6.360s
