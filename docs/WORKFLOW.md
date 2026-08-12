# Workflow

The public code is organized by scientific stage. Original execution compatibility is preserved in `config/SCRIPT_COMPATIBILITY_MAP.csv`. The configuration helper reconstructs a temporary working tree and substitutes neutral path tokens; scientific statements, thresholds, estimands and seeds are not changed.

The lightweight workflow validates included Source Data. The full workflow first validates external data paths, reports missing inputs and supports a dry run. Public packaging did not execute the complete raw-data workflow.
