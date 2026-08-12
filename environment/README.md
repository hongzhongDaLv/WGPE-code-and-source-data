# Software environment

The recorded analysis environment used R 4.5.2, Python 3.12.13, Node.js 24.14.0, Rtools45/GCC 14.3.0 and Windows x64. Exact package and session records are provided in this directory.

`R_package_versions.csv` is the R package record. A complete validated `renv.lock` is not available; restore the recorded versions and verify platform-specific system requirements. `requirements.txt` pins Python packages only where an unambiguous version was recorded. C++ helpers are compiled through Rcpp and require a C++17-capable compiler. A full rerun may require more than 100 GB for provider data and intermediate files.
