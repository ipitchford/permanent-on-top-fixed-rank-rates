# AI index: Sharp exponential rates for permanent-on-top and anomalous boson bunching at fixed rank (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) ("What
this is about") and the paper [manuscript.pdf](manuscript.pdf) (source [manuscript.tex](manuscript.tex); Section 1
states the question, why it matters, what was known, the results, what is known and new, and a roadmap; Section 8 lists
open problems; Appendix B describes the computational checks).
DOI 10.5281/zenodo.24550760 (concept 10.5281/zenodo.24550759). Companion: Evidence Press release
doi:10.5281/zenodo.24509057 (bapat-sunder-eigenvalue-ratio-growth; first compound, logarithmic scale).
Status: unrefereed candidate; internal AI-assisted review, proof audits and computation only.

## Notation

A = (v_i^* v_j): Gram matrix of unit vectors v_1..v_n in K^d, K = R or C ("configuration"; d = Gram rank). per A =
permanent. Π(A): Schur power, n!×n!, Π(A)_{σ,τ} = Π_i a_{σ(i),τ(i)}; per A is an eigenvalue. S_{K,d}(n) = sup
‖Π(A)‖/per A over configurations of n vectors in K^d. H_j harmonic number, ψ digamma, γ Euler's constant.
g_C(d) = e^{H_{d-1}}/d; g_R(1) = 1, g_R(d) = (2/d) exp{ψ(d/2) − ψ(1)}. M_t(v) = max over states X of rank ≤ t of
Π_i v_i^* X v_i; q(v) = M_1(v) (pure states, complex u even for real data). D_{n,j} = C(n+j−1, j−1).
c* = least constant with ‖Π(A)‖ ≤ (c*)^n per A for all n and all PSD A.

## Claims and where they are proved

| ID | Claim | Location | Status / dependence |
|---|---|---|---|
| P | q ≤ per A/n! ≤ D_{n,d} q and M_d ≤ ‖Π(A)‖/n! ≤ D_{n,d²} M_d | Proposition 2.1 | written proof; pure side standard (Yuan-Parrilo Cor 1.2), mixed upper side via Christandl-König-Renner postselection |
| G | M_t ≤ g_K(s)^n q with s = min(d,t) (C) or min(d,2t) (R), and the base is sharp | Proposition 3.1 | written proof; complex case and sharpness due to Yuan-Parrilo (Thms 4.6, 5.1, 5.2; Pappas-Révész); real constant (real data, complex pure states) new; Yuan-Parrilo's real constant is 2 g_R |
| L | g_C and g_R (d ≥ 2) strictly increase to e^γ; g_R(2t) = g_C(t) | Lemma 3.2 | written proof + 30-digit checks for d < 60 |
| T1 | lim S_{K,d}(n)^{1/n} = g_K(d) for every fixed d; also with Gram rank exactly d. g_C(2) = e/2 ≈ 1.3591, g_R(3) = e²/6 ≈ 1.2315, g_R(2) = 1 | Theorem 1.1 | written proof |
| T2a | ‖Π(A)‖ ≤ C(n+d²−1, n) g_K(d)^n per A for every PSD A of rank ≤ d, every n | Theorem 1.2(a) | written proof; spot-checked on 36 random configurations |
| T2b | e^γ ≤ c* ≤ e^{1+γ} (also e^γ for real matrices) | Theorem 1.2(b) | lower bound new; upper bound Anari-Gurvits-Oveis Gharan-Saberi Thm 2 |
| T2c | rank d(n) → ∞, d(n) = o(√n), per A > 0: sup (‖Π(A)‖/per A)^{1/n} → e^γ over C and R | Theorem 1.2(c) | written proof; asymptotic rate only, not c* |
| C | Balanced permanental compounds C_β (t parts): sharp base g_C(min(d,t)) / g_R(min(d,2t)), also for the best 0/1 test; t = 2: Pate's 0/1 inequality (Thm 7) fails exponentially, base e/2 at every complex d ≥ 2 | Theorem 5.1 | written proof; layer-cake lemma from the companion (Lemma 2.1) and Lev Thm 2 |
| O | Multimode bunching ratio max_S per(A∘S)/per A, internal dimension t: same sharp bases; S = conj(A) attains g_K(d) for t ≥ d; two outputs, polarisation: base e/2 | Corollary 6.1, eq. (6.1) | written proof; standard lossless model of partial distinguishability |
| D | Tetrahedron: ratios 7/9, 377/405, 102469/74115 ≈ 1.3826 (first > 1 at twelve photons), 50.26 at forty photons; limit 27/16 per four photons (1.1398 per photon); identical-photon probability 305/6377292 ≈ 4.78e-5 at twelve photons. Octahedron (H, V, D, A, R, L): 14/15, 62731/28800 ≈ 2.178 (twelve photons), 64385639/10348800 ≈ 6.22 (eighteen); > 650 at 42, > 24,000 at 60; limit √(3/2) ≈ 1.2247 per photon; probability ≈ 3.0e-5 at twelve photons | Section 7 | exact rational arithmetic + Ryser cross-checks |
| R | Root ratio → g_K(d) iff directions equidistribute; bounded number of directions gives a strictly smaller base | Theorem A.1 | written proof |

## Evidence

| Check | Program | Expected output | Arithmetic |
|---|---|---|---|
| Constants (L, T1), monotonicity d < 60, limit at d = 10^6; Theorem 1.2(a) bound on 36 random frames (seed 20261125); Shchesnovich order-five example ratio 1.0168; tetrahedral and octahedral closed forms against Ryser; exact ratios and probabilities (D); incidence identity of Theorem 5.1 for p+q ≤ 8; maximal pure products 1/27 and 1/216 on a Bloch grid | [code/verify_exp_rates.py](code/verify_exp_rates.py) | [results/verification.json](results/verification.json): 28 checks, "ALL CHECKS PASSED", empty failure list | mpmath 40 digits (30-digit tolerance); Python fractions and integers; numpy floats with declared tolerances |

Negative control: changing the exact tetrahedral ratio 102469/74115 to 102470/74115 makes the program exit 1 under
python -O ([REPLAY_RECEIPT.md](REPLAY_RECEIPT.md); also run in CI). The finite calculations test constants, closed
forms and identities; they do not certify the asymptotic theorems.

Numbers quoted on the release page: e/2 ≈ 1.36; e²/6 ≈ 1.23; e^γ ≈ 1.78; e^{1+γ}; 62731/28800 ≈ 2.18 at twelve
photons; more than 24,000 at sixty photons; √(3/2) ≈ 1.22 per photon; 102469/74115 ≈ 1.38; about 3×10⁻⁵; 28 checks.

## Limits

Rates are asymptotic (fixed rank, n → ∞); no uniformity in d beyond Theorem 1.2(c). Whether c* = e^γ is open, as are
unbalanced compounds, the best base for few polarisation states with only the internal states restricted, and explicit
finite-order bounds. Optical values are exact predictions of the standard lossless model, not experiments; absolute
probabilities are small and source and detector efficiencies are not modelled. Literature search of limited scope.

## Replay and records

Replay commands: [README.md](README.md), [ENVIRONMENT.txt](ENVIRONMENT.txt), [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), CI
[.github/workflows/replay.yml](.github/workflows/replay.yml). Review: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). Prior art
and limitation triage: [RESEARCH_GATES.json](RESEARCH_GATES.json). Bibliography: [refs.bib](refs.bib). Licences:
[LICENSE.md](LICENSE.md). Citation: [CITATION.cff](CITATION.cff). Hashes: [MANIFEST.sha256](MANIFEST.sha256).
