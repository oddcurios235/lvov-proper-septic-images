# Reproducing the computational checks

The `main.tex` and `references.bib` files produce `main.pdf`; the
verification scripts reconstruct the finite matrix and representation data
used in the manuscript, from the included integer inputs.

## Build the article

From the archive root, with a LaTeX distribution providing `pdflatex`, `bibtex`,
`amsart`, `hyperref`, `needspace`, `microtype`, and standard AMS packages, together with `latexmk`:

```sh
bash build.sh
```

The build script runs `latexmk` and BibTeX. Alternatively, compile with `pdflatex`, run `bibtex`, then run `pdflatex` twice.

## Run the verification code

Install Python 3 and NumPy, and from the archive root run:

```sh
python3 -m pip install -r requirements.txt
bash run_checks.sh
```

The shell script invokes all five programs using `python3` and writes
five JSON records to `results/`, asserting the expected values:

- `code/verify_septic_trace.py`: 1260 free-word generators independent modulo
  32003, exact rational cyclic-trace rank 211, plus fixed repeated-substitution
  comparison ranks.
- `code/verify_rep_filtration.py`: S7 character, irreducible multiplicities,
  Young-subgroup invariants, and 762-dimensional all-assignment kernel.
- `code/verify_deep_minors.py`: the fixed integer input matrices in
  `data/deep_core_inputs.json`, reduced modulo 101, endpoint determinants and
  lexicographically selected nonzero normal minors at N=3,4,5, and one mixed
  witness at N=7. Exact selected zero-based rows are saved to
  `results/deep_minors.json`.
- `code/verify_deep_jacobians.py`: direct derivatives of the full seven-input
  polynomial and the expected traceless ranks for n=4,5,6.
- `code/verify_deep_m3.py`: exact matrix-unit computation of
  `diag(2,-2,0)` and `tr(D^2)=8` at the frozen n=3 tuple.

A nonzero determinant residue modulo 101 certifies a nonzero **integer**
determinant, and hence nonvanishing in characteristic zero. The modular
residues are not numerical claims about the full unreduced integer minors.
The universal all-size conclusion also requires the non-computational proofs
in `main.tex`.
