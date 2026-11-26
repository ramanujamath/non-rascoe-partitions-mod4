"""Section 2: from partitions to the cone count.

  * P_inf / Q_inf = 1 + 2L (mod 4);
  * H_l (defined as Q_inf sum q^{r(r+l)}/(-q;q)_r) equals the Hecke--Rogers
    double series (7) exactly (integer coefficients, n < 2500), and
    H_l = T_l (mod 2);
  * the congruence (4): B_l = H_l + 2 T_l L + 2 M_l (mod 4);
  * (5): [q^n] C_l = C(40n + (2l+1)^2) for both definitions;
  * Proposition 2.1: B_0 = H_0 + 2 C_0 and B_1 = H_1 + 2 C_1 + 2 U (mod 4);
  * the coefficient of H_l vanishes when v_p(40n + (2l+1)^2) is odd
    (p = 29, 61, 101, 109), and [q^n] U = 0 at odd valuation;
  * (10/p) = -1, chi(-1) = 1, chi(2) = -1, chi(5) = 1 for p = 21, 29 mod 40.
"""
import math
import sys

import numpy as np

import nrlib as N

R = N.Report("check_reduction")
L = 2500
M = 64




R.section("Series identities")
f = N.euler_f(1, L, M)
P = N.fft_mul(N.euler_f(2, L, M), N.series_inv(f, M, L), M, L)  # (-q;q)_inf
Lser = np.zeros(L, dtype=np.int64)
for j in range(1, L):
    Lser[j::j] += 1
Lser %= M
ratio = N.fft_mul(P, N.series_inv(f, M, L), M, L)
one2L = (2 * Lser) % M
one2L[0] = (one2L[0] + 1) % M
R.check(all((ratio - one2L) % 4 == 0), "P_inf/Q_inf = 1 + 2L (mod 4)")


def sum_over_r(l, numer_power_fn, denom, sign):
    """sum_r sign^r * c_r * q^{r(r+l)} / denom_r, denom 'plus' = (-q;q)_r, 'minus' = (q;q)_r."""
    S = np.zeros(L, dtype=np.int64)
    term = np.zeros(L, dtype=np.int64)
    term[0] = 1
    r = 0
    while r * (r + l) < L:
        S = (S + numer_power_fn(r) * term) % M
        r += 1
        shift = 2 * r - 1 + l
        new = np.zeros(L, dtype=np.int64)
        if shift < L:
            new[shift:] = (sign * term[:L - shift]) % M
        if denom == "plus":
            term = N.div_one_plus(new, r, M)
        else:  # divide by (1 - q^r): cumulative sums along stride r
            rows = -(-L // r)
            y = np.zeros(rows * r, dtype=np.int64)
            y[:L] = new
            term = (np.cumsum(y.reshape(rows, r), axis=0) % M).reshape(-1)[:L].copy()
    return S


for l in (0, 1):
    H = N.fft_mul(f, sum_over_r(l, lambda r: 1, "plus", 1), M, L)
    HR = N.H_hecke_rogers(l, L)
    R.check(all((int(H[n]) - HR[n]) % M == 0 for n in range(L)), f"H_{l}: Q_inf sum q^(r(r+l))/(-q;q)_r agrees with the Hecke-Rogers series (7) (mod 64, n < {L})")
    T = np.zeros(L, dtype=np.int64)
    for m in range(-200, 201):
        e = m * (5 * m + 2 * l + 1) // 2
        if 0 <= e < L:
            T[e] = (T[e] + (-1) ** (m % 2)) % M
    R.check(all((H - T) % 2 == 0), f"H_{l} = T_{l} (mod 2)")
    Mser = N.fft_mul(f, sum_over_r(l, lambda r: r, "minus", 1), M, L)
    Bl = N.B_series(l, L, M)
    rhs = (H + 2 * N.fft_mul(T, Lser, M, L) + 2 * Mser) % M
    R.check(all((Bl - rhs) % 4 == 0), f"(4): B_{l} = H_{l} + 2 T_{l} L + 2 M_{l} (mod 4)")
    # C_l from its definition
    C = np.zeros(L, dtype=np.int64)
    if l == 0:
        a = 1
        while a * (5 * a - 1) // 2 < L:
            for base, upto in ((a * (5 * a - 1) // 2, a - 1), (a * (5 * a + 1) // 2, a)):
                for u in range(1, upto + 1):
                    for k in range(1, L):
                        e = base + u * k
                        if e >= L:
                            break
                        C[e] += 1
            a += 1
        Useries = np.zeros(L, dtype=np.int64)
    else:
        a = 0
        while a * (5 * a + 3) // 2 < L:
            for base in (a * (5 * a + 3) // 2, (a + 1) * (5 * a + 2) // 2):
                for u in range(1, a + 1):
                    for k in range(1, L):
                        e = base + u * k
                        if e >= L:
                            break
                        C[e] += 1
            a += 1
        Useries = np.array(N.U_series(L), dtype=np.int64)
    Csieve = N.cone_count_sieve(40 * L + 9)
    R.check(all(C[n] == Csieve[40 * n + (2 * l + 1) ** 2] for n in range(L)),
            f"(5): [q^n] C_{l} = C(40n + {(2 * l + 1) ** 2}) for n < {L}")
    R.check(all((Bl - H - 2 * C - 2 * Useries) % 4 == 0), f"Proposition 2.1: B_{l} = H_{l} + 2 C_{l}" + (" + 2U" if l else "") + " (mod 4)")
    ok = True
    for p in (29, 61, 101, 109):
        for n in range(L):
            if N.vp(40 * n + (2 * l + 1) ** 2, p) % 2:
                ok &= HR[n] == 0 and Useries[n] == 0
    R.check(ok, f"[q^n] H_{l} = 0 and [q^n] U = 0 when v_p(40n + {(2 * l + 1) ** 2}) is odd (p = 29, 61, 101, 109)")

R.section("Legendre symbols")
ps = [p for p in range(3, 3000) if N.is_prime(p) and p % 40 in (21, 29)]
R.check(all(N.legendre(10, p) == -1 and N.legendre(-1, p) == 1 and N.legendre(2, p) == -1 and N.legendre(5, p) == 1 for p in ps),
        f"(10/p) = -1, chi(-1) = 1, chi(2) = -1, chi(5) = 1 for all {len(ps)} primes p = 21, 29 mod 40 below 3000")
R.check(all(N.legendre(-2, p) == -1 and N.legendre(-10, p) == -1 for p in ps), "(-2/p) = (-10/p) = -1 for the same primes")

sys.exit(0 if R.done() else 1)
