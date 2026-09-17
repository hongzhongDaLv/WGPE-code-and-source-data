# September update execution records

The two ecology calculations were originally run with R 4.5.2. Concurrent-model spatial-block confidence intervals use 1,000 resamples and seed 20260904, set explicitly in `run_no_vpd.R`. The lag-3 classifications are deterministic given the canonical inputs. No new ecological OOF calculation is included.

The packaged calculation entrypoints remove workstation-specific R library and Rtools path overrides. `run_no_vpd.R` requires ggplot2 for its diagnostic plots; openxlsx is optional. The scientific calculations and `compute_f_no_vpd.R` use base R. The full approved figure chain additionally uses the packages and external inputs specified in its scripts and the existing environment ledgers.

For this packaging update, all 19 added R scripts were syntax-parsed with R 4.5.2. Python source-data checks were executed with the local bundled Python runtime and Pillow. This syntax/source-data verification is distinct from full raw-data recomputation or a fresh package-environment validation.
