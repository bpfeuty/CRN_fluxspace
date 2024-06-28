# metafluxspace_BJ2024
This repository contains the program, executable and data files necessary to reproduce the results of the scientific publication entitled "Free-energy transduction mechanisms shape the flux space of metabolic networks" by Benjamin Pfeuty, submitted to Biophysical Journal.

The file extensions are:
- ".f90" for fortran files
- ".py" for python files
- ".out" for compiled fortran files
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
- "figXY" for data file associated to Fig X, panel Y.
  
The fortran compiler used :
GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0
Copyright (C) 2021 Free Software Foundation, Inc.

Some Fortran programs require the LAPACK, MINPACK, and SEULEX libraries. These libraries can be installed on a Linux system and linked to gfortran using the following commands:
> apt-get install liblapack-dev

> apt-get install minpack-dev

> cp libseulex2.a /usr/lib/

> gfortran code.f90 -L/path/to/lib -llapack -lminpack -lseulex2
 
The files are organized into folders associated to a set of figures:
- Folder FIG_2-S1 (random reaction network with two chemostats) contains:
  - Code files: sampling_fig2_nm1.f90 and sampling_fig2_nm1.f90
  - Data files (Fig2B): dis_n8_nm1.dat dis_n16_nm1.dat dis_n8_nm2.dat, dis_n8_k001.dat dis_n8_k01.dat dis_n8_k1.dat dis_n8_k10.dat
  - Data files (Fig2C): fig2C-reac_nm1.dat, fig2C-reac_nm2.dat, fig2C-mu_nm1.dat, fig2C-mu_nm2.dat
- Folder FIG_3-4-S2 (random reaction network with three chemostats) contains:
  - manifold_fig34.f90 (generates manifold_fig3.dat and manifold_fig4.dat)
  - manifold_fig3.dat (used in fig3BCG)
  - manifold_fig4.dat (used in fig4BCG)
  - sampling_fig34_nmX.f90 (samples steady state flux solutions and the distribution of objective values)
  - fig3b_nm1.dat
  - fig3c_nm2.dat
  - fig4b_nm1.dat
  - fig4c_nm2.dat 
  - xxx
- Folder FIG_5-6-S3 (coarse grained metabolic network) contains:
  - manifold_fig56.f90 (generates the boundary manifold of the feasible flux space)
  - sampling_fig56.f90  (samples steady state flux solutions)
  - xxxx.dat
