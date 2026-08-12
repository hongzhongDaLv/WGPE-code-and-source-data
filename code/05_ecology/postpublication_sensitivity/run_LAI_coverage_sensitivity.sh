#!/usr/bin/env bash
set -euo pipefail
export PATH="/x86_64-w64-mingw32.static.posix/bin:/usr/bin:/bin:${PATH}"
exec "/c/Program Files/R/R-4.5.2/bin/Rscript.exe" \
  "__WGPE_PROJECT_ROOT__/CANONICAL_WGPE_TOP100_v2/07_SCRIPTS/LAI_coverage_sensitivity_FULL.R"
