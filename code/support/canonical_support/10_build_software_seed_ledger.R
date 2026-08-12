options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("__WGPE_R_LIBRARY__", .libPaths()))

ROOT <- "__WGPE_PROJECT_ROOT__"
CANON <- file.path(ROOT, "CANONICAL_WGPE_TOP100_v2")
PIPE <- file.path(ROOT, "scripts_TOP100_revision")
SCRIPT <- file.path(
  CANON, "07_SCRIPTS", "10_build_software_seed_ledger.R"
)
OUT <- file.path(
  CANON, "08_REPRODUCTION_OUTPUT", "reproducibility_ledger"
)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required.")

RUN_ALL <- file.path(PIPE, "run_all_TOP100_revision.R")
CONFIG <- file.path(PIPE, "00_config", "config_TOP100_revision.json")
SESSION <- file.path(
  ROOT, "output_TOP100_revision", "logs", "session_info.txt"
)

read_text <- function(path) {
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

runner_text <- read_text(RUN_ALL)
matches <- regmatches(
  runner_text,
  gregexpr(
    '"[^"]+\\.(?:R|py|mjs)"',
    runner_text,
    perl = TRUE
  )
)[[1]]
runner_rel <- unique(gsub('^"|"$', "", matches))
runner_paths <- file.path(PIPE, runner_rel)
runner_paths <- runner_paths[file.exists(runner_paths)]

canonical_paths <- list.files(
  file.path(CANON, "07_SCRIPTS"),
  pattern = "\\.(R|cpp|ps1|sh)$",
  recursive = TRUE,
  full.names = TRUE
)
formal_scripts <- unique(c(RUN_ALL, runner_paths, canonical_paths))
formal_scripts <- formal_scripts[file.exists(formal_scripts)]

language_from_path <- function(path) {
  ext <- tolower(tools::file_ext(path))
  switch(
    ext,
    r = "R",
    py = "Python",
    cpp = "C++",
    mjs = "Node.js",
    ps1 = "PowerShell",
    sh = "shell",
    ext
  )
}

script_ledger <- data.frame(
  script = normalizePath(formal_scripts, winslash = "/", mustWork = TRUE),
  language = vapply(formal_scripts, language_from_path, character(1)),
  size_bytes = file.info(formal_scripts)$size,
  sha256 = vapply(
    formal_scripts,
    digest::digest,
    character(1),
    algo = "sha256",
    file = TRUE
  ),
  stringsAsFactors = FALSE
)

seed_patterns <- paste(
  c(
    "set\\.seed\\s*\\(",
    "\\bSEED\\b\\s*(?:<-|=)",
    "\\bseed\\b\\s*(?:<-|=)",
    "random_seed",
    "random_state",
    "default_rng\\s*\\(",
    "np\\.random\\.seed\\s*\\(",
    "np\\.random\\.default_rng\\s*\\("
  ),
  collapse = "|"
)

seed_rows <- list()
for (path in formal_scripts) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  for (i in seq_along(lines)) {
    z <- trimws(lines[i])
    if (!nzchar(z) || grepl("^(#|//)", z)) next
    if (grepl(seed_patterns, z, perl = TRUE, ignore.case = TRUE)) {
      seed_rows[[length(seed_rows) + 1L]] <- data.frame(
        script = normalizePath(path, winslash = "/", mustWork = TRUE),
        language = language_from_path(path),
        line = i,
        seed_expression = z,
        stringsAsFactors = FALSE
      )
    }
  }
}
config_lines <- readLines(CONFIG, warn = FALSE, encoding = "UTF-8")
for (i in grep("random_seed", config_lines, fixed = TRUE)) {
  seed_rows[[length(seed_rows) + 1L]] <- data.frame(
    script = normalizePath(CONFIG, winslash = "/", mustWork = TRUE),
    language = "JSON",
    line = i,
    seed_expression = trimws(config_lines[i]),
    stringsAsFactors = FALSE
  )
}
seed_ledger <- if (length(seed_rows)) {
  unique(do.call(rbind, seed_rows))
} else {
  data.frame(
    script = character(),
    language = character(),
    line = integer(),
    seed_expression = character()
  )
}

script_ledger$seed_record_count <- vapply(
  script_ledger$script,
  function(path) sum(seed_ledger$script == path),
  integer(1)
)

r_paths <- formal_scripts[tolower(tools::file_ext(formal_scripts)) == "r"]
r_packages <- character()
for (path in r_paths) {
  txt <- read_text(path)
  ns <- regmatches(
    txt,
    gregexpr(
      "\\b[A-Za-z][A-Za-z0-9.]*::",
      txt,
      perl = TRUE
    )
  )[[1]]
  if (length(ns)) r_packages <- c(r_packages, sub("::$", "", ns))
  lib <- regmatches(
    txt,
    gregexpr(
      "(?:library|require)\\s*\\(\\s*[\"']?[A-Za-z][A-Za-z0-9.]*",
      txt,
      perl = TRUE
    )
  )[[1]]
  if (length(lib)) {
    r_packages <- c(
      r_packages,
      sub(
        ".*\\(\\s*[\"']?",
        "",
        lib
      )
    )
  }
  req <- regmatches(
    txt,
    gregexpr(
      "requireNamespace\\s*\\(\\s*[\"'][A-Za-z][A-Za-z0-9.]*[\"']",
      txt,
      perl = TRUE
    )
  )[[1]]
  if (length(req)) {
    r_packages <- c(
      r_packages,
      sub(
        ".*[\"']([A-Za-z][A-Za-z0-9.]*)[\"']$",
        "\\1",
        req
      )
    )
  }
}
r_packages <- sort(unique(r_packages))
r_packages <- r_packages[
  nzchar(r_packages) &
    grepl("^[A-Za-z][A-Za-z0-9.]*$", r_packages)
]
r_packages <- setdiff(
  r_packages,
  c(
    "base", "compiler", "datasets", "graphics", "grDevices", "grid",
    "methods", "parallel", "splines", "stats", "stats4", "tools", "utils"
  )
)
r_package_rows <- lapply(r_packages, function(pkg) {
  installed <- requireNamespace(pkg, quietly = TRUE)
  data.frame(
    language = "R",
    package = pkg,
    version = if (installed) as.character(packageVersion(pkg)) else NA_character_,
    status = if (installed) "available_in_archival_R_environment" else
      "not_available_in_current_archival_R_environment",
    environment = paste(.libPaths(), collapse = ";"),
    stringsAsFactors = FALSE
  )
})

dist_roots <- c(
  "__WGPE_PYTHON_SITE_PACKAGES__",
  file.path(ROOT, "scripts_core_validation", "_pydeps_netcdf_readable"),
  file.path(ROOT, "scripts_core_validation", "_pydeps_grib_readable")
)
dist_rows <- list()
for (dist_root in dist_roots[file.exists(dist_roots)]) {
  dirs <- list.dirs(dist_root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("\\.dist-info$", dirs, ignore.case = TRUE)]
  for (d in dirs) {
    nm <- basename(d)
    stem <- sub("\\.dist-info$", "", nm, ignore.case = TRUE)
    pos <- regexpr("-[0-9]", stem)
    if (pos < 1) next
    pkg <- substr(stem, 1, pos - 1L)
    ver <- substr(stem, pos + 1L, nchar(stem))
    dist_rows[[length(dist_rows) + 1L]] <- data.frame(
      language = "Python",
      package = pkg,
      version = ver,
      status = "version_recovered_from_dist_info",
      environment = normalizePath(
        dist_root,
        winslash = "/",
        mustWork = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }
}
python_package_rows <- if (length(dist_rows)) {
  unique(do.call(rbind, dist_rows))
} else {
  data.frame(
    language = character(),
    package = character(),
    version = character(),
    status = character(),
    environment = character()
  )
}

package_ledger <- rbind(
  do.call(rbind, r_package_rows),
  python_package_rows
)
package_ledger <- package_ledger[
  order(package_ledger$language, tolower(package_ledger$package),
        package_ledger$environment),
]

cmd_version <- function(command, args = "--version") {
  if (!file.exists(command)) return(NA_character_)
  z <- tryCatch(
    system2(command, args, stdout = TRUE, stderr = TRUE),
    error = function(e) NA_character_
  )
  if (!length(z)) NA_character_ else z[1]
}

python_exe <- "__WGPE_PYTHON__"
node_exe <- "__WGPE_NODE__"
gpp_exe <- "__WGPE_RTOOLS__/x86_64-w64-mingw32.static.posix/bin/g++.exe"
session_lines <- if (file.exists(SESSION)) {
  readLines(SESSION, warn = FALSE, encoding = "UTF-8")
} else {
  character()
}
recorded_python <- sub(
  "^python=",
  "",
  session_lines[grepl("^python=", session_lines)]
)

runtime_ledger <- data.frame(
  component = c(
    "R", "Python", "Node.js", "C++ compiler", "Operating system"
  ),
  version = c(
    R.version.string,
    if (length(recorded_python)) recorded_python[1] else
      cmd_version(python_exe),
    cmd_version(node_exe),
    cmd_version(gpp_exe),
    paste(Sys.info()[c("sysname", "release", "version")], collapse = " ")
  ),
  provenance = c(
    "current archival runtime and analysis QA files",
    "output_TOP100_revision/logs/session_info.txt",
    "bundled archival runtime",
    "Rtools45 compiler used by Rcpp",
    "current archival workstation"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  script_ledger,
  file.path(OUT, "ANALYSIS_SCRIPT_LEDGER.csv"),
  row.names = FALSE
)
write.csv(
  seed_ledger,
  file.path(OUT, "RANDOM_SEED_LEDGER.csv"),
  row.names = FALSE
)
write.csv(
  package_ledger,
  file.path(OUT, "SOFTWARE_PACKAGE_VERSION_LEDGER.csv"),
  row.names = FALSE
)
write.csv(
  runtime_ledger,
  file.path(OUT, "SOFTWARE_RUNTIME_LEDGER.csv"),
  row.names = FALSE
)

stochastic_pattern <- paste(
  c(
    "\\bsample\\s*\\(", "\\bsample\\.int\\s*\\(",
    "\\brunif\\s*\\(", "\\brnorm\\s*\\(",
    "np\\.random", "default_rng\\s*\\(",
    "\\.permutation\\s*\\(", "\\.shuffle\\s*\\("
  ),
  collapse = "|"
)
stochastic_scripts <- vapply(
  formal_scripts,
  function(path) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    lines <- trimws(lines)
    lines <- lines[nzchar(lines) & !grepl("^(#|//)", lines)]
    grepl(
      stochastic_pattern,
      paste(lines, collapse = "\n"),
      ignore.case = TRUE,
      perl = TRUE
    )
  },
  logical(1)
)
stochastic_paths <- normalizePath(
  formal_scripts[stochastic_scripts],
  winslash = "/",
  mustWork = TRUE
)
unresolved_seed_scripts <- setdiff(
  stochastic_paths,
  unique(seed_ledger$script)
)
seed_audit <- data.frame(
  script = stochastic_paths,
  seed_recorded = stochastic_paths %in% seed_ledger$script,
  stringsAsFactors = FALSE
)
write.csv(
  seed_audit,
  file.path(OUT, "STOCHASTIC_SCRIPT_SEED_AUDIT.csv"),
  row.names = FALSE
)

qa_lines <- c(
  "# Software and random-seed ledger QA",
  "",
  paste0("- Formal scripts inventoried: ", nrow(script_ledger), "."),
  paste0("- Seed expressions recorded: ", nrow(seed_ledger), "."),
  paste0("- Package/version records: ", nrow(package_ledger), "."),
  paste0("- Runtime records: ", nrow(runtime_ledger), "."),
  paste0("- Scripts containing stochastic operations: ",
         length(stochastic_paths), "."),
  paste0("- Stochastic scripts without an in-script seed expression: ",
         length(unresolved_seed_scripts), "."),
  "",
  paste0(
    "Python versions were recovered from the recorded session snapshot and ",
    "project/bundled `.dist-info` metadata. R package versions were queried ",
    "from the archival R 4.5.2 library."
  )
)
if (length(unresolved_seed_scripts)) {
  qa_lines <- c(
    qa_lines,
    "",
    "## Stochastic scripts requiring contextual review",
    "",
    paste0("- `", unresolved_seed_scripts, "`")
  )
}
qa_path <- file.path(OUT, "SOFTWARE_SEED_LEDGER_QA.md")
writeLines(qa_lines, qa_path, useBytes = TRUE)

manifest_path <- file.path(OUT, "REPRODUCIBILITY_LEDGER_SHA256.csv")
output_files <- setdiff(list.files(OUT, full.names = TRUE), manifest_path)
all_files <- c(SCRIPT, RUN_ALL, CONFIG, SESSION, output_files)
all_files <- all_files[file.exists(all_files)]
manifest <- data.frame(
  path = all_files,
  size_bytes = file.info(all_files)$size,
  sha256 = vapply(
    all_files,
    digest::digest,
    character(1),
    algo = "sha256",
    file = TRUE
  ),
  stringsAsFactors = FALSE
)
write.csv(manifest, manifest_path, row.names = FALSE)

message("REPRODUCIBILITY_LEDGER=COMPLETE")
message("OUTPUT=", OUT)
