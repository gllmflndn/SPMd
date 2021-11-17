# SPMd toolbox

[![Build & Test](https://github.com/gllmflndn/SPMd/actions/workflows/main.yml/badge.svg)](https://github.com/gllmflndn/SPMd/actions/workflows/main.yml)

## Statistical Parametric Mapping Diagnosis

SPMd is a toolbox for SPM which you can use to establish the validity of inferences in fMRI modeling through diagnosis of linear model assumptions, and to characterize fMRI signal and artifacts through exploratory data analysis.

  * https://www.nisox.org/Software/spmd/
  * https://www.nisox.org/Software/spmd/spmdexample

## Publications

```
Wen-Lin Luo and Thomas E. Nichols
Diagnosis and exploration of massively univariate neuroimaging models.
NeuroImage 19:1014-1032, 2003.
```
```
Wen-Lin Luo and Thomas E. Nichols
Diagnosis and Exploration of Massively Univariate 4D Spatiotemporal Models.
Technical Report.
```

## Source code

### Interface

  * [`spm_diagnosis.m`](spm_diagnosis.m): (old) GUI
  * [`tbx_cfg_spmd.m`](tbx_cfg_spmd.m): (WIP) batch interface

### Computation

  * [`spmd_comp_MS.m`](spmd_comp_MS.m): see `spm_spm.m` and `spm_searchlight.m`
     * [`spmd_stat.m`](spmd_stat.m)
  * [`spmd_comp_SS.m`](spmd_comp_SS.m)
  * [`spmd_MkResid.m`](spmd_MkResid.m): see `spm_write_residuals.m`
    * [`spmd_resss.m`](spmd_resss.m): see `compat/spm_resss.m`

### Visualisation

  * [`spmd_MD.m`](spmd_MD.m)
    * [`spmd_MD_plot.m`](spmd_MD_plot.m)
      * [`spmd_Nplot.m`](spmd_Nplot.m)
  * [`spmd_MS.m`](spmd_MS.m)
  * [`spmd_SD.m`](spmd_SD.m)
  * [`spmd_SS.m`](spmd_SS.m)

### Utilities

  * [`spmd_check_registration.m`](spmd_check_registration.m): see `spm_check_registration.m`
    * [`spmd_SptlBar.m`](spmd_SptlBar.m)
  * [`spmd_mtsview.m`](spmd_mtsview.m)
  * [`spmd_orthviews.m`](spmd_orthviews.m): see `spm_orthviews.m`
  * [`spmd_pointer.m`](spmd_pointer.m): see `spm_XYZreg.m`
  * [`spmd_getTS.m`](spmd_getTS.m): (unused) see `spm_get_data.m`