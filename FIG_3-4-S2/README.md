This folder contains files associated to Figures 3 and 4 and supplementary Figure S2 of our study.

**FEASIBLE SOLUTION SPACES IN FIGURES 3,4 AND S2**
- use "manifold_fig34.f90" to generate "manifold_fig3.dat" and "manifold_fig4.dat"
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
  - fig4g.dat (nr=7 must be selected)

    
-- dis_n8_nm1.dat, dis_n8_kX.dat, dis_n16_kX.dat
- use "sampling_fig2_nm2.f90" to generate:
-- dis_n8_nm1.dat, dis_n8_kX.dat
  
**FIG2C:**
- use "sampling_fig2_nm1.f90" to generate:
-- dis_n8_nm1.dat, dis_n8_kX.dat, dis_n16_kX.dat
- use "sampling_fig2_nm2.f90" to generate:
-- dis_n8_nm1.dat, dis_n8_kX.dat
  
**FIGS1:**
- use "sampling_fig2_nm1.f90" to generate:
-- all distributions.
-- hyperparameters and ranges are iteratively changed
