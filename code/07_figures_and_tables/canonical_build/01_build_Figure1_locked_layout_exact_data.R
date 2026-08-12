options(stringsAsFactors = FALSE, warn = 2)

# Canonical reproducibility wrapper. The final plotting script reads only
# canonical source tables; this wrapper redirects generated outputs to a
# dedicated reproduction area without touching the locked release figures.
ROOT <- "__WGPE_PROJECT_ROOT__"
CANONICAL <- file.path(ROOT, "CANONICAL_WGPE_TOP100_v2")
SCRIPT <- file.path(CANONICAL, "07_SCRIPTS", "build_Figure1_TOP100_FINAL_v5.R")
OUT <- file.path(CANONICAL, "08_REPRODUCTION_OUTPUT", "Figure1")

stopifnot(file.exists(SCRIPT))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TOP100_CLOSURE_ROOT = OUT)
if (!nzchar(Sys.getenv("TOP100_NO_PANEL_LABELS", unset = ""))) {
  Sys.setenv(TOP100_NO_PANEL_LABELS = "0")
}
source(SCRIPT, encoding = "UTF-8")
