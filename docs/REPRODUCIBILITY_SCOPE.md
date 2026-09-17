# Reproducibility scope

## Included and tested

- Reading and schema validation of compact final Source Data.
- Cross-checking headline global <z> and W-GPE trends against the public registry.
- Exact closure of the three W-GPE trend channels from included values.
- Verification of one Figure 2e OOF increment group and one ecology summary.
- Regeneration of a low-cost vertical-contribution test panel.

## Requires third-party data

Raw-data reconstruction of gridded ERA5/ERA5-Land diagnostics, JRA-55 and RHARM validation, GPP/LAI analyses, field significance and full-resolution maps requires the external products listed in `data/DATA_PROVENANCE.csv`.

## Not claimed

The complete multi-day raw-data analysis was not rerun during public packaging. The full-analysis wrapper validates inputs and prepares audited scripts, but does not claim an unattended end-to-end pass.
# September 2026 update

The current five-control ecology, paired VPD sensitivity and updated approved figure entrypoints are documented in [SEPTEMBER_2026_UPDATE.md](SEPTEMBER_2026_UPDATE.md). Source-data numerical tests are runnable from this repository. Refitting the ecology models requires external monthly canonical RDS inputs; exact full-map rendering additionally requires project adapters and legacy gridded inputs. The original scope below applies to the v1.0.0 modules.
