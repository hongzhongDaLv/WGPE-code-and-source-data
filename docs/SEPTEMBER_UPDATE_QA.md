# September update validation

- 379/379 additional numerical and semantic source-data checks passed.
- 7/7 original source-data smoke checks passed, including the explicitly historical six-control ecology check.
- 19/19 added R scripts parsed successfully with R 4.5.2.
- Published figures are unchanged copies of the approved Figure 1, aligned Figure 2 and Supplementary Figure S1 PNGs.
- Original current-panel files are preserved under `source_data/historical/v1.0.0/`; v1.0.0 tags/releases are untouched.

The additional tests reconstruct area-weighted means and significant-area fractions from valid gridded samples; verify fixed samples across adjustment sets; reconstruct temporal grid/area percentages; and verify current registry mapping and channel-interval ordering. Their machine-readable results are regenerated in `tests/output/SEPTEMBER_UPDATE_CHECKS.csv`.

These checks validate the packaged evidence. They do not assert that the full raw-data pipeline or all external plotting dependencies were rerun.
