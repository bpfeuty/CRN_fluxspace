# metafluxspace_BJ2024
This repository contains the program, executable and data files necessary to reproduce the results of the scientific publication entitled "Free-energy transduction mechanisms shape the flux space of metabolic networks" by Benjamin Pfeuty, submitted to Biophysical Journal.

The file extensions are:
- ".f90" for fortran files
- ".py" for python files
- ".out" for compiled fortran files
- ".dat" for data files
- ".xxx" for fortran library file
  
The script file names are respectively:
- "manifold." for computing the feasible flux space of positive entropy production rate.
- "sampling." for computing the sampled flux solutions
- "figXY." for data file associated to FigX panel Y
  
The fortran compiler used is GNU fortran compiler. 

Below is the list of all files with their description
The repository contains the fortran and python codes (.f90 and .py), the executables () and data associated to the computation of thermodynamically-feasible flux spaces and the codes and data associated to sampling analysis.

To compute the  of diverse metabolic networks and chemostatting conditions, the repository contains:
- the fortran script: manifold_fig3-4.f90  and manifoled_fig5.f90
- the data files:  manifold_fig3.dat, manifold_fig4.dat, manifold_fig5.dat

Regarding the kinetic sampling of the feasible flux space of diverse metabolic networks, the repository contains:
- the fortran script: sampling_fig2.f90, sampling_fig3-4.f90  and sampling_fig5.f90
  The compilation requires three fortran libraries : lapack minpack and seulex.
- the data files: sampling_fig2.dat, sampling_fig3.dat, sampling_fig4.dat, sampling_fig5.dat

Besides these core programs that performs manifold analysis and sampling analysis, the repository also contains programs that performs specific analysis.

