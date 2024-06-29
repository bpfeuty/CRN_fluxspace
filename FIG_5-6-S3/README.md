This folder contains files associated to Figures 5 and 6 and supplementary Figure S3 of our study.

**FEASIBLE SOLUTION SPACES IN FIGURES 5,6 AND S3**
- "manifold_fig56.f90" generates
  - boundary_fig5.dat (boundary manifold in high dimension space)
  - jbexglu_fig5.dat, jbexace_fig5.dat, jbexox_fig5.dat, jbexglu_fig5.dat (projection in 2D flux spaces)
- "nullspaces.py" computes:
  - left/right null space vectors
  - number of closed/emergent reaction cycles and number of broken/unbroken conservation laws

**SAMPLED SOLUTION SPACES IN FIGURES 5,6 AND S3**
- sampling_fig56.f90 generates (depending in the selected value i0)
  - fig56_samp.dat (default stoichiometry)
  - fig6b1.dat (high ATP yield)
  - fig6b3.dat (low ATP yield)
  - columns of those files reads: Ts,J(1-12),eta_gly,eta_ares,Ts/Ts*
  - these data are used in Fig5c, Fig6a, Fig6b and Fig6c (different columns are plotted)

**OPTIMAL BIOMASS YIELD IN FIGURES 6 AND S3**
- maxyield_fig6.f90 (that use a random optimization sceme) generates
  - fig6d.dat 
  - columns of this files reads: S_{atp,gly},S_{atp,ares},beta^{opt},eta_gly^{opt},eta,ares^{opt},ts^{opt}
  - data of figure S3 by varying a single feature.
