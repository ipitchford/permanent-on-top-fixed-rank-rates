# Review summary (internal, AI-assisted)

This is a public summary of the internal checking of version 0.1.0-candidate, without scores. It is not peer review,
independent reproduction or formal verification. No human referee has assessed the paper.

## What was done (25 November 2026)

1. **A combined five-perspective Claude review** (editor, methodology, domain, perspective, devil's advocate) in one
   pass, as is usual for short papers in this series. The methodology reviewer checked every proof line by line and ran
   independent scratch computations (constants, both comparisons on random small configurations, the incidence identity,
   the multinomial mode bound, the binary-test lemma, the tetrahedral closed forms). No mathematical error was found.
   The panel asked for revisions of attribution, positioning and presentation, and suggested extensions.
2. **Clean replays** of the check program from a fresh copy under `python -O`, with output identical to the shipped
   result, and a negative control.
3. **One revision round with extensions**: a context-first introduction with a "known and new" paragraph; precise credit
   to Yuan and Parrilo, Pappas and Révész, and Christandl, König and Renner; the companion release cited; the Pate
   paragraph rewritten from the full text; the gap proposition stated for every rank parameter; an explicit finite bound
   at every order; the lower bound c* ≥ e^γ and the exact rate e^γ for slowly growing rank; the six-polarisation design
   with exact values; a uniform rigidity gap; a new check program.
4. **Verification**: a re-review of everything that changed and a proof audit of all new and rewritten proofs. No
   critical issue. One major item was an understatement (a theorem stated more weakly than its own proof gives); it was
   strengthened to rank o(√n). Minor items were applied before release. All cited sources were obtained and read.

No Astra/Codex second-model review took place in this batch.

## Main objections and how they were handled

| Objection | Disposition |
|---|---|
| The complex constant and its sharpness are Yuan and Parrilo's; what is new was unclear | "Known and new" paragraph in Section 1; abstract says the complex value was known |
| The companion release on the first compound was not cited | Cited and positioned (logarithmic versus exponential scale); its layer-cake lemma and rank-two appendix credited |
| Pate's inequality and conjecture were described before reading the full text | Paragraph rewritten from the full text (Theorem 7, Conjecture 2, Corollary 3) |
| The open question of Anari et al. on the best constant was not addressed | Theorem 1.2(b): e^γ ≤ c* ≤ e^{1+γ}; Theorem 1.2(c): exact rate e^γ for rank o(√n); full-rank question left open |
| The gap proposition was used beyond its stated range | Restated for every rank parameter t, with the effective rank in both fields |
| The experimental meaning of the optical result was not discussed | "Physical scale" paragraph: model predictions only, small absolute probabilities, efficiencies not modelled |
| Readability for non-specialists | New introduction (question, why it matters, what was known, results, roadmap); section openers; open problems |

## Limits that remain

The rates are asymptotic. Whether c* = e^γ is open, as are unbalanced compounds and finite-order bounds. The optical
values are predictions of the standard lossless model. The literature search had limited scope. The computational
checks test constants, closed forms and identities; the theorems rest on the written proofs.
