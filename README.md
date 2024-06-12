# metafluxspace_BJ2024
This repository contains the fortran scripts and data files allowing to generate the figures of the paper "Free-energy transduction mechanisms shape the flux space of metabolic networks" authored by Benjamin Pfeuty and submitted to Biophysical Journal.

The repository contains the codes and data associated to the computation of thermodynamically-feasible flux spaces and the codes and data associated to sampling analysis.

To compute the feasible flux space of diverse metabolic networks and chemostatting conditions, the repository contains:
- the fortran script: manifold_fig3-4.f90  and manifoled_fig5.f90
- the data files:  manifold_fig3.dat, manifold_fig4.dat, manifold_fig5.dat

Regarding the kinetic sampling of the feasible flux space of diverse metabolic networks, the repository contains:
- the fortran script: sampling_fig2.f90, sampling_fig3-4.f90  and sampling_fig5.f90
  The compilation requires three fortran libraries : lapack minpack and seulex.
- the data files: sampling_fig2.dat, sampling_fig3.dat, sampling_fig4.dat, sampling_fig5.dat

Besides these core programs that performs manifold analysis and sampling analysis, the repository also contains programs that performs specific analysis.

