# Full rerun guide

1. Obtain third-party products from official providers.
2. Copy `config/config.example.yml` to `config/config.yml` and set all paths.
3. Install the recorded R/Python packages and C++ toolchain described in `environment/README.md`.
4. Run `./run_full_analysis.ps1 -Config config/config.yml` for a dry run.
5. Resolve every `MISSING` entry. No formal output is produced when an input is absent.
6. Run with `-Execute` to validate and prepare the compatibility working tree.
7. Review `docs/WORKFLOW_STAGES.csv` and execute stages in order, monitoring storage and checkpoint outputs.

The release wrapper intentionally stops after safe preparation rather than silently launching a multi-day job. This boundary prevents partial outputs from being mistaken for a completed formal rerun.
