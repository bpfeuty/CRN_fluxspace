# metafluxspace_BJ2024
This repository contains the program, executable and data files necessary to reproduce the results of the scientific publication entitled "Free-energy transduction mechanisms shape the flux space of metabolic networks" by Benjamin Pfeuty, submitted to Biophysical Journal.

The file extensions are:
- ".f90" for fortran files
- ".py" for python files
- ".dat" for data files
- ".a" for fortran library file
- ".gnu" for gnuplot file used to generate figure ?
  
The file names begins respectively by:
- "manifold" for computing the feasible flux space of positive entropy production rate.
  - manifold_fig34.f90 and manifold_fig56.f90
  - no library need to be linked
- "sampling" for computing the sampled flux solutions
  - sampling_fig2_nmX.f90, sampling_fig34_nmX.f90 and sampling_fig56.f90
  - library used: lseulex2.a, lminpack.a, -llapack.a
- "figXy" for data file associated to Fig X, panel Y.
  
The fortran compiler used :
GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0
Copyright (C) 2021 Free Software Foundation, Inc.

Some Fortran programs require the LAPACK, MINPACK, and SEULEX libraries. These libraries can be installed on a Linux system and linked to gfortran using the following commands:
> apt-get install liblapack-dev

> apt-get install minpack-dev

> cp libseulex2.a /usr/lib/

> gfortran code.f90 -L/path/to/lib -llapack -lminpack -lseulex2
 
The files are organized into folders associated to a set of figures:
- FIG_2-S1 (random reaction network with two chemostats)
- FIG_3-4-S2 (random reaction network with three chemostats)
- FIG_5-6-S3 (coarse grained metabolic network) contains:
Each folder contains fortran files and gnuplot to generate data and figures, as well as a subfolder named DATA containing most data files.

