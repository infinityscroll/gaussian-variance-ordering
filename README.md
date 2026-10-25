# Gaussian variance-ordering counterexamples

Explicit constructions and proofs for alpha-divergence Gaussian variational inference, by Jitendra Prajapati.

[Read the paper](output/pdf/variance-reversals-and-entropy-ordering.pdf) · [LaTeX source](main.tex) · [Exact finite-example verifier](verify_exact.py)

For a fixed Gaussian target, increasing the alpha-divergence parameter above one need not increase every variance of the optimal factorized Gaussian approximation. This concerns the extension conjecture in Margossian, Pillaud-Vivien and Saul's [JMLR paper](https://jmlr.org/papers/v26/24-0878.html); its proved ordering theorem below one remains unaffected.

## Results in the note

- For every prescribed `a > 1`, an explicit positive-definite target has a unique global optimum with mixed-sign variance derivatives at `a`. The construction fixes the target before varying alpha. Its dimension is `ceil(8a/(a-1)) + 2`; a separate ten-dimensional example works at `a = 2`.
- A rational three-dimensional target has optimal covariance `D(9) = I`. Its second and third optimal variances are below `19/20` for every `alpha >= 10000`. Thus the explicit pair `(9, 10000)` gives a decrease exceeding 5% for the same target.
- Dimension three is minimal for any reversal. The two-dimensional monotonicity antecedent is credited to Li and Turner (2016), Appendix B.
- An entropy-ordering statement survives for fixed linear Gaussian precision families. The paper includes a projection proof and a short derivation from classical log-determinant barrier inequalities. It is presented as an elementary consequence of established tools, not a new fundamental barrier principle.

## Reproduce the finite algebra

Run with Python 3; no external packages are required:

```sh
python3 verify_exact.py
```

The expected output is `PASS_EXACT_FINITE_ALGEBRA` with 89 checks. All scientific arithmetic uses `fractions.Fraction`. Checks cover positive-definiteness minors, the exact stationarity and derivative identities at alpha 9, and a rational root bracket at alpha 10000. The saved output is [verification/exact-results.json](verification/exact-results.json).

The verifier is **not a formal proof checker**. Gaussian integration, global optimality, the intermediate value theorem, the arbitrary-alpha construction, and the general entropy theorem are analytic arguments in the paper. No floating-point optimizer or trained model was used.

## Scope and provenance

This is a public research note, first released on October 25, 2026. It makes no first-discovery or publication-priority claim. The dated private source material is distinguished from this public release in [PROVENANCE.json](PROVENANCE.json).

AI tools were used extensively for the literature search, derivations, writing, and checking. Separate AI proof reviews and exact algebra checks are not independent human peer review. No practical inference improvement, real-data miscalibration, or research-job outcome is established by these examples. A decreasing variance can still exceed the true target variance.

Please report mathematical errors with the relevant equation and a reproducible calculation or argument.

## Build the paper

```sh
pdflatex -interaction=nonstopmode -halt-on-error -output-directory output/pdf main.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory output/pdf main.tex
```

This produces `output/pdf/main.pdf`; the release copy is named `variance-reversals-and-entropy-ordering.pdf`.
