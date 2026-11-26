"""Checks for "Sharp exponential rates for permanent-on-top and anomalous boson bunching at fixed rank".

Run from the package root:
    python3 code/verify_exp_rates.py --output results/verification.json

The main results of the paper are proved analytically. This script checks the
constants, the exact finite examples of Section 7 and several identities used in
the proofs. Exact checks use integer or rational arithmetic; floating-point checks
state their tolerance. The script exits with status 1 if any check fails; it does
not use assert, so it behaves the same under python -O.
"""
import argparse
import itertools
import json
import math
import platform
import sys
from fractions import Fraction as F
from math import comb, factorial

import mpmath as mp
import numpy as np

mp.mp.dps = 40
RESULTS = {}
FAILURES = []


def check(name, ok, detail):
    RESULTS.setdefault("checks", []).append({"name": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + ": " + str(detail))
    if not ok:
        FAILURES.append(name)


# ---------------------------------------------------------------- constants
def g_C(d):
    return mp.e ** mp.harmonic(d - 1) / d


def g_R(d):
    if d == 1:
        return mp.mpf(1)
    return 2 * mp.e ** (mp.digamma(mp.mpf(d) / 2) - mp.digamma(1)) / d


def constants():
    eg = mp.e ** mp.euler
    table = []
    for d in range(1, 13):
        table.append({"d": d, "g_C": mp.nstr(g_C(d), 15), "g_R": mp.nstr(g_R(d), 15)})
    RESULTS["constants"] = {"table": table, "exp_gamma": mp.nstr(eg, 15),
                            "e_over_2": mp.nstr(mp.e / 2, 15), "e2_over_6": mp.nstr(mp.e ** 2 / 6, 15)}
    tol = mp.mpf(10) ** -30
    check("g_C(2)=e/2", abs(g_C(2) - mp.e / 2) < tol, mp.nstr(g_C(2), 15))
    check("g_R(3)=e^2/6", abs(g_R(3) - mp.e ** 2 / 6) < tol, mp.nstr(g_R(3), 15))
    check("g_R(2)=1", abs(g_R(2) - 1) < tol, mp.nstr(g_R(2), 15))
    check("g_R(2t)=g_C(t), t<=8", all(abs(g_R(2 * t) - g_C(t)) < tol for t in range(1, 9)), "t=1..8")
    inc_C = all(g_C(d) < g_C(d + 1) for d in range(1, 60))
    inc_R = all(g_R(d) < g_R(d + 1) for d in range(2, 60))
    check("g_C, g_R strictly increasing (d<60) and below e^gamma",
          inc_C and inc_R and g_C(60) < eg and g_R(60) < eg,
          "g_C(60)=%s g_R(60)=%s e^gamma=%s" % (mp.nstr(g_C(60), 10), mp.nstr(g_R(60), 10), mp.nstr(eg, 10)))
    check("g_C(d), g_R(d) -> e^gamma (d=10^6)",
          abs(g_C(10 ** 6) - eg) < 1e-5 and abs(g_R(10 ** 6) - eg) < 1e-5,
          "%s %s" % (mp.nstr(g_C(10 ** 6), 12), mp.nstr(g_R(10 ** 6), 12)))


# ---------------------------------------------------------------- permanents
def ryser_int(M):
    """Exact permanent of an integer matrix (Ryser with Gray code)."""
    n = len(M)
    rowsums = [0] * n
    total = 0
    sign_n = -1 if n % 2 else 1
    prev = 0
    for k in range(1, 1 << n):
        gray = k ^ (k >> 1)
        diff = gray ^ prev
        j = diff.bit_length() - 1
        add = 1 if gray & diff else -1
        for i in range(n):
            rowsums[i] += add * M[i][j]
        prev = gray
        p = 1
        for r in rowsums:
            p *= r
            if p == 0:
                break
        bits = bin(gray).count("1")
        total += (-1) ** bits * p
    return sign_n * total


def ryser_complex(A):
    n = A.shape[0]
    tot = 0j
    for S in range(1, 1 << n):
        cols = [j for j in range(n) if S >> j & 1]
        tot += (-1) ** len(cols) * np.prod(A[:, cols].sum(axis=1))
    return (-1) ** n * tot


# tetrahedral design (Section 7.1)
def tetra_vectors():
    w = np.exp(2j * np.pi / 3)
    return [np.array([1, 0], complex)] + [np.array([1, math.sqrt(2) * w ** j]) / math.sqrt(3) for j in range(3)]


def tetra_perA(m):
    return F(sum(comb(m, j) ** 2 * 8 ** j * factorial(4 * m - 3 * j) * factorial(3 * j) for j in range(m + 1)),
             3 ** (3 * m))


def tetra_perM(m):
    s = 0
    for h in itertools.product(range(m + 1), repeat=4):
        s += factorial(4 * m - sum(h)) * math.prod(comb(m, x) ** 2 * factorial(x) * 2 ** x for x in h)
    return F(s, 3 ** (4 * m))


def tetrahedron():
    vs = tetra_vectors()
    S = sum(np.outer(v, v.conj()) for v in vs)
    check("tetra: sum v v* = 2I", np.allclose(S, 2 * np.eye(2), atol=1e-13), "tol 1e-13")
    # pure product over the Bloch sphere: maximum 1/27 (grid + local refinement, float)
    best = 0.0
    for th in np.linspace(0, np.pi, 361):
        for ph in np.linspace(0, 2 * np.pi, 361):
            x = np.array([np.cos(th / 2), np.exp(1j * ph) * np.sin(th / 2)])
            best = max(best, np.prod([abs(x.conj() @ v) ** 2 for v in vs]))
    check("tetra: max pure product = 1/27 (grid, tol 1e-9)", abs(best - 1 / 27) < 1e-9 and best <= 1 / 27 + 1e-12,
          best)
    for m in (1, 2, 3, 4):
        V = np.column_stack([v for v in vs for _ in range(m)])
        A = V.conj().T @ V
        pa = ryser_complex(A).real
        pm = ryser_complex(A * A.conj()).real
        ok = abs(pa - float(tetra_perA(m))) < 1e-8 * float(tetra_perA(m)) and \
            abs(pm - float(tetra_perM(m))) < 1e-8 * float(tetra_perM(m))
        check("tetra: closed forms match Ryser, m=%d (rel tol 1e-8)" % m, ok,
              {"per_A": pa, "formula": float(tetra_perA(m)), "per_masked": pm, "formula_masked": float(tetra_perM(m))})
    ratios = {m: tetra_perM(m) / tetra_perA(m) for m in range(1, 41)}
    check("tetra: ratios m=1,2,3 = 7/9, 377/405, 102469/74115",
          ratios[1] == F(7, 9) and ratios[2] == F(377, 405) and ratios[3] == F(102469, 74115),
          [str(ratios[1]), str(ratios[2]), str(ratios[3])])
    check("tetra: ratio at m=10 exceeds 50.26", ratios[10] > F(5026, 100), float(ratios[10]))
    p3 = tetra_perA(3) / F(6 ** 12)
    p10 = tetra_perA(10) / F(20 ** 40)
    check("tetra: indistinguishable probability at m=3 is 305/6377292", p3 == F(305, 6377292), float(p3))
    check("tetra: indistinguishable probability at m=10 about 1.61e-18",
          F(160, 10 ** 20) < p10 < F(162, 10 ** 20), float(p10))
    W = np.column_stack([v / np.sqrt(2 * 3) for v in vs for _ in range(3)])
    check("tetra: rows of the 2 x 12 design matrix are orthonormal", np.allclose(W @ W.conj().T, np.eye(2), atol=1e-13),
          "tol 1e-13")
    roots = [float(ratios[m]) ** (1 / m) for m in (1, 2, 3, 5, 10, 20, 40)]
    RESULTS["tetrahedron"] = {
        "ratio": {str(m): str(ratios[m]) for m in (1, 2, 3)},
        "ratio_float": {str(m): float(ratios[m]) for m in (1, 2, 3, 4, 5, 10, 20, 40)},
        "ratio_root_per_repetition": dict(zip(["1", "2", "3", "5", "10", "20", "40"], roots)),
        "limit_per_repetition": 27 / 16, "limit_per_photon": (27 / 16) ** 0.25,
        "indistinguishable_probability": {"3": str(p3), "3_float": float(p3), "10_float": float(p10)},
    }
    check("tetra: first m with ratio > 1 is 3", ratios[2] < 1 < ratios[3], "")


# octahedral design (six standard polarisations; Section 7.2)
def poly_mul_linear(P, coeffs):
    """Multiply a polynomial {exponent tuple: complex int} by a linear form with Gaussian-integer coefficients."""
    Q = {}
    k = len(coeffs)
    for e, c in P.items():
        for i, a in enumerate(coeffs):
            if a == 0:
                continue
            e2 = e[:i] + (e[i] + 1,) + e[i + 1:]
            Q[e2] = Q.get(e2, 0) + c * a
    return {e: c for e, c in Q.items() if c != 0}


def fock_norm2(P):
    tot = 0
    for e, c in P.items():
        re, im = int(c.real), int(c.imag)
        tot += (re * re + im * im) * math.prod(factorial(x) for x in e)
    return tot


def octa_perA(m):
    return F(sum(comb(m, j) ** 2 * factorial(m + 4 * j) * factorial(5 * m - 4 * j) for j in range(m + 1)), 16 ** m)


def octa_masked_fock(mmax):
    # 2 (v ⊗ conj v) for H, V, D, A, R, L; coordinates (11, 12, 21, 22); Gaussian integers
    forms = [(2, 0, 0, 2 * 0), (0, 0, 0, 2), (1, 1, 1, 1), (1, -1, -1, 1), (1, -1j, 1j, 1), (1, 1j, -1j, 1)]
    forms = [tuple(complex(a) for a in f) for f in forms]
    out = {}
    P = {(0, 0, 0, 0): complex(1)}
    for m in range(1, mmax + 1):
        for f in forms:
            P = {e: complex(round(c.real), round(c.imag)) for e, c in poly_mul_linear(P, f).items()}
        out[m] = F(fock_norm2(P), 4 ** (6 * m))
    return out


def octahedron():
    s = 1 / math.sqrt(2)
    vs = [np.array([1, 0], complex), np.array([0, 1], complex), np.array([s, s], complex), np.array([s, -s], complex),
          np.array([s, 1j * s]), np.array([s, -1j * s])]
    S = sum(np.outer(v, v.conj()) for v in vs)
    check("octa: sum v v* = 3I", np.allclose(S, 3 * np.eye(2), atol=1e-13), "tol 1e-13")
    best = 0.0
    for th in np.linspace(0, np.pi, 361):
        for ph in np.linspace(0, 2 * np.pi, 361):
            x = np.array([np.cos(th / 2), np.exp(1j * ph) * np.sin(th / 2)])
            best = max(best, np.prod([abs(x.conj() @ v) ** 2 for v in vs]))
    check("octa: max pure product = 1/216 (grid, tol 1e-6)", abs(best - 1 / 216) < 1e-6 and best <= 1 / 216 + 1e-12,
          best)
    masked = octa_masked_fock(10)
    for m in (1, 2, 3):
        V = np.column_stack([v for v in vs for _ in range(m)])
        A = V.conj().T @ V
        pa = ryser_complex(A).real
        # exact masked permanent: 2 (A o conj A) has entries 2 (same vector), 0 (antipodal), 1 (otherwise)
        M2 = [[int(round(2 * abs(A[i, j]) ** 2)) for j in range(6 * m)] for i in range(6 * m)]
        pm = F(ryser_int(M2), 2 ** (6 * m))
        ok = abs(pa - float(octa_perA(m))) < 1e-8 * float(octa_perA(m)) and pm == masked[m]
        check("octa: per A_m closed form (Ryser, rel tol 1e-8) and exact masked permanent, m=%d" % m, ok,
              {"per_A": pa, "formula": float(octa_perA(m)), "masked_ryser": str(pm), "masked_fock": str(masked[m])})
    ratios = {m: masked[m] / octa_perA(m) for m in masked}
    RESULTS["octahedron"] = {
        "ratio": {str(m): str(ratios[m]) for m in (1, 2, 3)},
        "ratio_float": {str(m): float(ratios[m]) for m in ratios},
        "ratio_root_per_photon": {str(m): float(ratios[m]) ** (1 / (6 * m)) for m in ratios},
        "limit_per_repetition": 27 / 8, "limit_per_photon": math.sqrt(1.5),
        "indistinguishable_probability_float": {str(m): float(octa_perA(m) / F((3 * m) ** (6 * m))) for m in ratios},
    }
    first = min(m for m in ratios if ratios[m] > 1)
    RESULTS["octahedron"]["first_m_ratio_above_one"] = first
    check("octa: ratio exceeds one from m=%d (6m photons)" % first, all(ratios[m] > 1 for m in ratios if m >= first),
          {str(m): round(float(ratios[m]), 6) for m in ratios})
    W = np.column_stack([v / np.sqrt(3) for v in vs])
    check("octa: rows of the 2 x 6 design matrix (scaled by 1/sqrt(3m)) are orthonormal",
          np.allclose(W @ W.conj().T, np.eye(2), atol=1e-13), "tol 1e-13")


# ---------------------------------------------------------------- Schur power, small orders
def schur_power(A):
    n = A.shape[0]
    perms = list(itertools.permutations(range(n)))
    P = np.array([[np.prod([A[s[i], t[i]] for i in range(n)]) for t in perms] for s in perms])
    return P


def small_orders():
    rng = np.random.default_rng(20261125)
    worst = 0.0
    rows = []
    for field in ("C", "R"):
        for d in (2, 3):
            for n in (4, 5, 6):
                for _ in range(3):
                    if field == "C":
                        V = rng.normal(size=(d, n)) + 1j * rng.normal(size=(d, n))
                    else:
                        V = rng.normal(size=(d, n)) + 0j
                    A = V.conj().T @ V
                    P = schur_power(A)
                    top = np.linalg.eigvalsh(P)[-1]
                    per = P.sum(axis=1)[0].real
                    g = float(g_C(d) if field == "C" else g_R(d))
                    bound = comb(n + d * d - 1, n) * g ** n
                    r = top / per
                    worst = max(worst, r / bound)
                    rows.append({"field": field, "d": d, "n": n, "ratio": float(r), "bound": bound})
    check("Theorem 1.2(a) bound ||Pi(A)|| <= C(n+d^2-1,n) g(d)^n per A on 36 random frames (float)", worst <= 1,
          "max ratio/bound = %.4g" % worst)
    # the order-5 rank-2 complex counterexample used in the archive check
    V = np.array([[4 + 2j, 2 - 3j, -4 - 4j, -3 + 4j, 1], [2 - 4j, -3j, 2 - 4j, -3j, -5 - 7j]])
    A = V.conj().T @ V
    P = schur_power(A)
    top = np.linalg.eigvalsh(P)[-1]
    per = P.sum(axis=1)[0].real
    check("order-5 rank-2 complex example: ||Pi||/per > 1 (float, rel tol 1e-12)", top / per > 1 + 1e-6,
          {"per": per, "top": float(top), "ratio": float(top / per)})
    RESULTS["small_orders"] = {"random": rows, "order5_example_ratio": float(top / per)}


# ---------------------------------------------------------------- incidence identity (Theorem 5.1)
def incidence():
    ok = True
    for N in range(1, 9):
        for p in range(0, N):
            q = N - p
            if q <= p:
                continue
            Sp = list(itertools.combinations(range(N), p))
            Sp1 = list(itertools.combinations(range(N), p + 1))
            U = np.array([[1 if set(a) <= set(b) else 0 for a in Sp] for b in Sp1], dtype=np.int64)
            lhs = U.T @ U
            if p > 0:
                Sm = list(itertools.combinations(range(N), p - 1))
                Um = np.array([[1 if set(a) <= set(b) else 0 for a in Sm] for b in Sp], dtype=np.int64)
                lhs = lhs - Um @ Um.T
            ok = ok and np.array_equal(lhs, (q - p) * np.eye(len(Sp), dtype=np.int64))
    check("incidence identity U_p^T U_p - U_{p-1} U_{p-1}^T = (q-p) I, p+q <= 8 (exact integers)", ok, "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default="results/verification.json")
    args = ap.parse_args()
    constants()
    tetrahedron()
    octahedron()
    small_orders()
    incidence()
    RESULTS["environment"] = {"python": sys.version.split()[0], "numpy": np.__version__, "mpmath": mp.__version__,
                              "platform": platform.platform()}
    RESULTS["failures"] = FAILURES
    with open(args.output, "w") as fh:
        json.dump(RESULTS, fh, indent=1, sort_keys=True)
    print("ALL CHECKS PASSED" if not FAILURES else "FAILED: " + ", ".join(FAILURES))
    sys.exit(1 if FAILURES else 0)


if __name__ == "__main__":
    main()
