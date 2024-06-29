This folder contains files associated to Figures 3 and 4 and supplementary Figure S2 of our study.

**FEASIBLE SOLUTION SPACES IN FIGURES 3,4 AND S2**
- "manifold_fig34.f90" generates "manifold_fig3.dat" and "manifold_fig4.dat"
- columns of data files reads: J1ex,J2ex,j3ex,Ts,Ts/Ts*

**SAMPLED SOLUTION SPACES IN FIGURES 3,4 AND S2**
- sampling_fig34_nm1.f90 generates (depending in the selected value i0)
  - fig3b_nm1.dat
  - fig3e_nec2_nm1.dat (nbc=1) 
  - fig4b_nm1.dat
  - fig4e_ncc2_nm1.dat
- sampling_fig34_nm2.f90 generates (depending in the selected value i0)
  - fig3c_nm2.dat
  - fig3e_nec2_nm2.dat (nbc=1)
  - fig3e_nec1_nm2.dat (nbc=2)
  - fig3g_nm2_100K.dat
  - fig4c_nm2.dat
  - fig4e_ncc2_nm3.dat
  - fig4e_ncc0_nm2.dat
  - figS2.dat (need to adjust the hyperparameters and sampled parameter ranges)
- minimal_fig34.f90 generates (depending in the selected value ifig)
  - fig3g.dat (nr=5 must be selected)
  - values mu and J written in Fig3F.
  - fig4g.dat (nr=7 must be selected)
