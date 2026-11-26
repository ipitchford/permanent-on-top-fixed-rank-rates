# Replay receipt (v0.1.0-candidate, 26 November 2026)

Release replay by the publisher from a fresh copy of the approved candidate (code, requirements.txt), macOS arm64,
working directory = package root, one CPU, about 70 seconds. Interpreter as printed inside the run (see
[ENVIRONMENT.txt](ENVIRONMENT.txt)): CPython 3.12.11 with numpy 2.3.5 and mpmath 1.3.0.

| Program | Flags | Exit | Output |
|---|---|---|---|
| code/verify_exp_rates.py | -I -B -O | 0 | results/verification.json byte-identical to the shipped file; 28 PASS lines; "ALL CHECKS PASSED" |
| Negative control | the same program with the exact tetrahedral ratio 102469/74115 changed to 102470/74115, -I -B -O | 1 | the tetrahedral-ratio check fails |

stderr was empty for the passing run. The output matches the lead's replay of 25 November 2026. After the replay only
the title-page DOI line of manuscript.tex, the rebuilt PDF and the release records changed; code, results and
requirements.txt are byte-identical to the approved candidate. GitHub CI on Linux reruns the program under the
comparison policy in ENVIRONMENT.txt, with the negative control.
