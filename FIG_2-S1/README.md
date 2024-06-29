This folder contains files associated to Figure 2 and supplementary Figure S1 of our study.

**FIG2B:**
- "sampling_fig2_nm1.f90" generates:
  - dis_n8_nm1.dat, dis_n8_kX.dat, dis_n16_kX.dat
- "sampling_fig2_nm2.f90" generates:
  - dis_n8_nm2.dat
  
**FIG2C:**
- "sampling_fig2_nmX.f90" (X=1,2 for left/right panels) generates:
  - fig2C-mean_nmX.dat (columns are Ts,Ts^I,mu_1-mu_2)
  - fig2C-reac_nmX.dat (columns are Ts,Ts_r,|DG_r|,|J_r| where r are all reactions)
  - fig2C-reacb_nmX.dat (columns are Ts,Ts_r,|DG_r|,|J_r| where r are exchange reactions)
  - fig2C-mu_nmX.dat (columns are Ts,mu_i where i are all species)
  - fig2C-mub_nmX.dat (columns are Ts,mu_i where i are exchanged species)

**FIGS1:**
- "sampling_fig2_nm1.f90" generates:
-   distributions where hyperparameters and sampled parameter ranges are varied one at a time
