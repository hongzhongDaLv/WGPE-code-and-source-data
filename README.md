# WGPE code and Source Data

September 2026 development update (**v1.1.0-dev**), based on v1.0.0 (2026-08-11), for **Global moistening accompanies an upward shift in the atmospheric water-vapour centroid**. The v1.0.0 release remains unchanged; this branch does not assign a new release DOI.

The current ecology analysis uses five controls (precipitation, T2m, SSRD, root-zone soil moisture and CO2). New gridded results, paired VPD sensitivity, lag-3 classifications and approved figure scripts/PNGs are included. Read [the update and reproduction scope](docs/SEPTEMBER_2026_UPDATE.md) before combining current and historical results.

## Scientific purpose

The moisture-mass-weighted geopotential height, `<z>`, describes where atmospheric water-vapour mass is centred vertically. W-GPE is the gravitational-potential-energy-scaled first moment of column water vapour, `g × IWV × <z>`. The two diagnostics complement IWV by separating moisture amount and mean vertical position.

## Repository contents

- `code/`: original analysis/build scripts plus September five-control calculations and the approved figure dependency chain.
- `source_data/`: compact author-generated values used for figures, tables and headline checks.
- `data/`: third-party data provenance and acquisition/redistribution notices.
- `environment/`: recorded runtime, package and random-seed information.
- `docs/`: workflow, script-to-output map and reproduction boundaries.
- `tests/`: a low-cost Source Data smoke test.
- `figures/approved_20260916/`: unchanged approved Figure 1, Figure 2 and Supplementary Figure S1 PNGs.

## Five-minute start

1. Install Python 3.12 and Pillow.
2. Run `./run_source_data_workflow.ps1` on Windows or `sh run_source_data_workflow.sh` on POSIX systems.
3. Review `tests/output/SMOKE_TEST_REPORT.md`, `tests/output/SEPTEMBER_UPDATE_REPORT.md` and the regenerated test panel. The September checks recompute global/regional means and grid-versus-area percentages from the supplied grids.

## Full analysis

The included Source Data reproduce headline numerical checks without downloading ERA5. Full maps and raw-data analyses require the products in `data/DATA_PROVENANCE.csv` and the exact files listed in `config/config.example.yml`. Copy the example to `config/config.yml`, set local paths and run `./run_full_analysis.ps1 -Config config/config.yml` without `-Execute` for a dry-run. Use `-Execute` only after the preflight confirms that all inputs are present.

## Environment and licences

The recorded environment used R 4.5.2, Python 3.12.13, Node.js 24.14.0 and Rtools45/GCC 14.3.0 on Windows x64. Full recomputation is data- and storage-intensive. Code is MIT-licensed. Author-generated Source Data listed in `source_data/SOURCE_DATA_INDEX.csv` use CC BY 4.0. Third-party products retain their provider terms and are not redistributed.

## Citation

Use `CITATION.cff` and the identifier shown on the published release record.
