#!env bash
set -o pipefail
set -e # Exit on error

R --vanilla -q -e 'devtools::install_github("morinlab/ascat/ASCAT", dependencies = TRUE, upgrade = "never")'
R --vanilla -q -e 'devtools::install_github("morinlab/battenberg", dependencies = FALSE, upgrade = "never")'