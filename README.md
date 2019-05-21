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



## Speedup improvements:

+ before changes  
calculating only p(q)  
real    57m0.304s  
user    357m27.428s  
sys     0m11.928s  


+ After Changes  
Collapsed do loop into arrays  
Improved OMP scheduling  
Expanded use of double_precision variables  
Also calculating p_ind(q)  
real    39m31.588s
user    182m7.256s
sys     0m4.876s
