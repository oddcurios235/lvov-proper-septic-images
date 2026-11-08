#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
# Fixed source date makes the LaTeX PDF reproducible across clean builds.
export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-1794096000}"
export FORCE_SOURCE_DATE=1
if ! command -v latexmk >/dev/null 2>&1; then
  echo 'A TeX installation with latexmk is required.' >&2
  exit 1
fi
if command -v bibtex >/dev/null 2>&1; then
  TEX_BIBTEX="$(command -v bibtex)"
elif command -v bibtex.original >/dev/null 2>&1; then
  TEX_BIBTEX="$(command -v bibtex.original)"
else
  echo 'This build requires bibtex or bibtex.original.' >&2
  exit 1
fi
# Use the selected BibTeX executable with latexmk.
latexmk -e "\$bibtex = '$TEX_BIBTEX %O %B';" \
  -pdf -interaction=nonstopmode -halt-on-error main.tex
