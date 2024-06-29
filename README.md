# metafluxspace_BJ2024

This repository contains the program, executable and data files necessary to reproduce the results of the scientific publication entitled "Free-energy transduction mechanisms shape the flux space of metabolic networks" by Benjamin Pfeuty, submitted to Biophysical Journal.

**File extensions:**
- ".f90" for fortran files
- ".py" for python files
- ".dat" for data files
- ".a" for fortran library file
- ".gnu" for gnuplot file
- ".pdf" for portable document format
  
**File naming convention:**
- "manifold" for computing the feasible flux space of positive entropy production rate.
  - manifold_fig34.f90 and manifold_fig56.f90
- "sampling" for computing the sampled flux solutions
  - sampling_fig2_nmX.f90, sampling_fig34_nmX.f90 and sampling_fig56.f90
- "figX" for figure, gnuplot and data files associated to Figure X.
  
**Fortran compiler used :**
GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0
Copyright (C) 2021 Free Software Foundation, Inc.

**Required libraries:**
Some Fortran programs require the LAPACK, MINPACK, and SEULEX libraries. These libraries can be installed on a Linux system and linked to gfortran using the following commands:
> apt-get install libblas-dev liblapack-dev minpack-dev

> cp libseulex2.a /usr/lib/

> gfortran code.f90 -L/path/to/lib -llapack -lminpack -lseulex2
 
**Folder organization:**
- **FIG_2-S1/** : random reaction network involving two chemostats
- **FIG_3-4-S2/** : random reaction network involving three chemostats
- **FIG_5-6-S3/** : coarse grained metabolic network
  
Each folder contains Fortran and Gnuplot files to generate data and figures (fig*.pdf), as well as a subfolder named **DATA/** containing most data files.
