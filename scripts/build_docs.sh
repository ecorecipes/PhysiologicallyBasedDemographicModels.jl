#!/usr/bin/env bash
# Build the Documenter.jl HTML site.
#
# Usage:
#   bash scripts/build_docs.sh                 # HTML only (default)
#   PBDM_BUILD_PDF=true  bash scripts/build_docs.sh   # also build a PDF
#   PBDM_RERENDER=true   bash scripts/build_docs.sh   # re-run quarto render
#                                                     # for each vignette
#
# Output: docs/build/  (open docs/build/index.html in a browser)
set -euo pipefail
cd "$(dirname "$0")/.."

if [ ! -f docs/Manifest.toml ]; then
  julia --project=docs -e '
    using Pkg
    Pkg.develop([Pkg.PackageSpec(path="../StructuredPopulationCore.jl"),
                 Pkg.PackageSpec(path=".")])
    Pkg.add("Documenter")
    Pkg.instantiate()'
fi

julia --project=docs docs/make.jl
echo
echo "Docs built: docs/build/index.html"
