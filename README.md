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

The files are organized into folders associated to a set of figures:
- Folder FIG_2-S1 (random reaction network with two chemostats) contains:
  - sampling_fig2_nmX.f90
  - fig2B.dat, fig2C.dat, fig2D.dat
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
