# A Zariski-open family of proper septic polynomials with full traceless images

**Y. T. Miah** · Research manuscript · 8 November 2026

[Read the manuscript (PDF)](main.pdf) · [LaTeX source](main.tex) · [Reproducing the finite checks](REPRODUCIBILITY.md)

## About

This manuscript studies actual sets of values of multilinear matrix polynomials in the L’vov–Kaplansky problem. Over an algebraically closed field $K$ of characteristic zero, it identifies a **single nonempty Zariski-open subset** of a 762-dimensional coefficient space such that every polynomial $p$ in that subset satisfies

$$
\{p(A_0,\ldots,A_6): A_i\in M_n(K)\}=\mathfrak{sl}_n(K)
\qquad\text{for every }n\ge 3.
$$

These are equalities of sets of **individual evaluations**. The coefficient space consists of proper multilinear degree-seven polynomials in the third power of the commutator ideal, with vanishing formal cyclic trace and annihilation under *every* repeated-input assignment of multiplicities $(3,3,1)$.

The argument combines symmetric-group representation calculations, explicit finite determinant witnesses, a size-transfer argument and a separate treatment of $3\times3$ matrices. It builds on matrix-image work of Kanel-Belov, Malev and Rowen, image criteria of Fagundes and Koshlukov, results of Pan and Tang, and Miah’s distinct, [previously published standard-polynomial paper](https://hexagonmath.org/2611.00184v1). The explicit polynomial $q_0$ witnesses nonemptiness of the open set; the full-image property of that individual polynomial also follows from Pan’s degree-six actual-image theorem and the classical single-commutator theorem.

## Files

- `main.pdf`, `main.tex`, `main.bbl`, `references.bib` — manuscript, LaTeX source and bibliography.
- `code/verify_*.py` — five finite verifiers for the representation and cyclic-trace calculations, common kernel, normal-map determinants, full Jacobians and matrix-unit computation.
- `data/deep_core_inputs.json` — fixed integer seed matrices and vectors.
- `results/*.json` — five recorded finite-check certificates.
- `REPRODUCIBILITY.md`, `requirements.txt`, `run_checks.sh`, `build.sh` — detailed instructions and runnable commands.

## Build and checks

With a TeX distribution providing `latexmk`, pdfLaTeX, BibTeX, `amsart`, `hyperref`, `needspace`, `microtype` and standard AMS packages, run from the repository root:

```sh
bash build.sh
```

With Python 3 and NumPy installed, run the five finite verifiers:

```sh
python3 -m pip install -r requirements.txt
bash run_checks.sh
```

They reconstruct the recorded representation, kernel, determinant, derivative and matrix-unit data, including selected minor residues $8$, $91$ and $37$ modulo $101$ for the three seed sizes. The proofs in the manuscript establish the transfer to all matrix sizes. See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the individual checks and their precise scope.

## Citation

Y. T. Miah, *A Zariski-open family of proper septic polynomials with full traceless images*, manuscript dated 8 November 2026. Machine-readable citation details appear in [CITATION.cff](CITATION.cff).

## Licensing

Original scholarly material and data are licensed under [CC BY 4.0](LICENSES/CC-BY-4.0.txt); original verification code and related scripts are under [BSD 3-Clause](LICENSES/BSD-3-Clause.txt). See [LICENSE.md](LICENSE.md) for the exact file scope.
