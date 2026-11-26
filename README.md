# Sharp exponential rates for permanent-on-top and anomalous boson bunching at fixed rank

Version 0.1.0-candidate, 26 November 2026. DOI [10.5281/zenodo.24550760](https://doi.org/10.5281/zenodo.24550760)
(all versions: [10.5281/zenodo.24550759](https://doi.org/10.5281/zenodo.24550759)). Author: Anonymous.
Release page: https://evidencepress.org/releases/permanent-on-top-fixed-rank-rates/

**Status: unrefereed research candidate.** The results are written proofs with computational checks. They were checked
by internal AI-assisted review, proof audits and replays only; see [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). This is not
peer review, independent reproduction or formal verification.

## What this is about

The *permanent* of a square matrix is the determinant without minus signs. For a positive semidefinite matrix A (a
table of inner products of n vectors), a much larger companion matrix, the *Schur power* Π(A) of size n! by n!, always
has the permanent as one of its eigenvalues. The *permanent-on-top* conjecture (Soules, 1966) said it is the largest
one. If true, it would have implied a whole family of permanent inequalities. Shchesnovich showed in 2016 that it is
false, which left a quantitative question: by how much can it fail?

The size of the failure matters in three places: it bounds how badly that family of inequalities can fail; it appears
in the theory of approximating permanents, where Anari, Gurvits, Oveis Gharan and Saberi asked for the best constant;
and in linear optics it limits *anomalous bunching*, where photons that are partly distinguishable (for example by
polarisation) are more likely than identical photons to leave an interferometer together through a few chosen ports.

## What the paper shows

Fix the rank d of A (the dimension of the vectors) and let the order n grow.

- **Sharp base at every fixed rank.** The worst ratio ‖Π(A)‖/per A grows like an exact exponential base to the power n:
  e^{H_{d-1}}/d for complex matrices and 2·exp{ψ(d/2) − ψ(1)}/d for real ones. That is e/2 ≈ 1.36 at complex rank two
  and e²/6 ≈ 1.23 at real rank three; both rise to e^γ ≈ 1.78 as d grows (γ is Euler's constant).
- **The general constant.** The best constant c* with ‖Π(A)‖ ≤ (c*)^n per A for all matrices satisfies
  e^γ ≤ c* ≤ e^{1+γ} (before, only c* > 1 was known; the upper bound is Anari et al.'s). When the rank grows more slowly
  than √n the exact rate is e^γ. An explicit bound holds at every finite order and rank.
- **Permanental compounds.** Balanced compounds have the same sharp bases, so an inequality stated by Pate (2008) for
  0/1 vectors fails exponentially.
- **Optics.** In the standard lossless model, two output ports and ordinary polarisation already give exponential
  bunching enhancement with sharp base e/2. With the six standard polarisations (H, V, D, A, R, L) the model enhancement
  is 62731/28800 ≈ 2.18 at twelve photons and exceeds 24,000 at sixty photons.
- **Rigidity.** The full base is reached only when the vectors' directions become evenly spread.

**Credit.** The complex value of the key optimisation gap (products of quadratic forms over mixed versus pure states)
and its sharpness are due to Yuan and Parrilo, building on Pappas and Révész; the link from the Schur power to that gap
uses the postselection technique of Christandl, König and Renner. New here: the link itself, the real constant, the
rates, the lower bound c* ≥ e^γ, and the compound, optical and rigidity consequences with explicit designs.

**Limits.** The rates are asymptotic (fixed rank, n → ∞). Whether c* = e^γ at full rank is open. The optical numbers
are exact predictions of the standard model, not experiments, and the absolute probabilities are small (about 3×10⁻⁵
at twelve photons). Companion release on the first compound (logarithmic growth):
[doi:10.5281/zenodo.24509057](https://doi.org/10.5281/zenodo.24509057).

## Files

| File | Contents |
|---|---|
| [manuscript.pdf](manuscript.pdf), [manuscript.tex](manuscript.tex) | The paper (14 pages) and its LaTeX source |
| [code/verify_exp_rates.py](code/verify_exp_rates.py) | 28 checks: constants to 30 digits, closed forms against Ryser's formula, exact ratios and probabilities, identities |
| [results/verification.json](results/verification.json) | Output of the program |
| [refs.bib](refs.bib) | Bibliography (23 entries) |
| [AI_INDEX.md](AI_INDEX.md) | Machine-oriented map of claims, evidence and limits |
| [RESEARCH_GATES.json](RESEARCH_GATES.json) | Prior-work map, sources read, limitation triage |
| [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md) | Public summary of the internal review, without scores |
| [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), [ENVIRONMENT.txt](ENVIRONMENT.txt) | Release replay and environment |
| [CITATION.cff](CITATION.cff), [LICENSE.md](LICENSE.md), [MANIFEST.sha256](MANIFEST.sha256) | Citation, licences, SHA-256 sums |

## Reproduce

From the package root (about 70 seconds, one CPU):

```sh
uv run --python 3.12 --with numpy==2.3.5 --with mpmath==1.3.0 python3 -I -B -O code/verify_exp_rates.py --output results/verification.json
```

Expected: exit 0, 28 PASS lines and "ALL CHECKS PASSED". The program does not use `assert`, so it behaves the same
under `python -O`. On the release host the output was byte-identical to the shipped file. Across platforms the
comparison policy in [ENVIRONMENT.txt](ENVIRONMENT.txt) applies: the `environment` field and each check's free-text
`detail` are excluded; check names, pass flags, strings (all exact rationals) and integers must match exactly;
floating-point fields must agree within absolute 1e-12 or relative 1e-9. GitHub CI applies this policy on Linux and
runs a negative control.

## Licence

Text, paper, results and records: CC0-1.0. Code: MIT. See [LICENSE.md](LICENSE.md).
