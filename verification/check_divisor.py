"""Sections 6-8: square discriminants.

  * Proposition 6.1 (gamma(m) = C(m^2) and the divisor transform A_m);
  * Proposition 6.2, exactly as an identity of integer power series;
  * P_{2j} = T_j + T_{j+5} - T_{2j}; Table 2 and its integer formula for
    2 w(m, s); the Liouville-type evaluations S_6 = 0, S_8 = J - H_5,
    S_10 = J - H_10; the integer identity 2W(N) = K_5(N) - Q(N); (11);
  * the characters of Section 7 (lambda(2) = i, lambda(-1) = lambda(5) = -1,
    lambda^2 = Legendre symbol mod 29) and the congruence (12);
  * Theorem 7.1 for every odd N divisible by 29 below 2.1*10^6 (FFT);
  * Proposition 8.1 ingredients: R(N) as a representation count,
    R(N) = R(N/29^2), q_0(d) + q_0(29d) = 2 mod 4, Q(N) - Q(N/29^2) = 2 tau(m);
  * Theorem 8.2 for every odd r coprime to 5 with r <= 345;
  * the group index 276480 and the Sturm bound 46080.
"""
import math
import sys

import numpy as np

import nrlib as N

R = N.Report("check_divisor")
NMAX = 2001
P = 29


def mobius(n):
    res = 1
    d = 2
    while d * d <= n:
        if n % d == 0:
            n //= d
            if n % d == 0:
                return 0
            res = -res
        d += 1
    return -res if n > 1 else res


def tau(n):
    return len(N.divisors(n))


def d5(n):
    return sum(1 for d in N.divisors(n) if d % 5)


def nu5(n):
    return sum(1 for d in N.divisors(n) if d % 5 in (2, 3))


R.section("Proposition 6.1")
Csieve = N.cone_count_sieve(NMAX ** 2)
gamma = {m: int(Csieve[m * m]) for m in range(1, NMAX) if m % 2 and m % 5}


def lambert(d, L):
    s = np.zeros(L, dtype=np.int64)
    s[d::d] = 1
    return s


Aser = np.zeros(NMAX, dtype=np.int64)
for a in range(5, NMAX, 5):
    for s in range(1, NMAX):
        if 2 * a * s >= NMAX:
            break
        Lrest = NMAX - 2 * a * s
        prod = np.convolve(lambert(a, Lrest), lambert(2 * s, Lrest))[:Lrest] + \
            np.convolve(lambert(2 * a, Lrest), lambert(s, Lrest))[:Lrest]
        Aser[2 * a * s:] += prod
ok1 = all(Aser[m] == sum(gamma[d] for d in N.divisors(m)) for m in gamma)
ok2 = all(gamma[m] == sum(mobius(d) * Aser[m // d] for d in N.divisors(m)) for m in gamma)
R.check(ok1 and ok2, f"A_m = sum_(d|m) gamma(d) and gamma(m) = sum mu(d) A_(m/d) for all odd m < {NMAX} coprime to 5")

R.section("Proposition 6.2 (exact)")
Lm = NMAX


def Ud(d):
    s = np.zeros(Lm, dtype=np.int64)
    s[d::2 * d] = 1
    return s


U = sum(Ud(d) for d in range(1, Lm, 2))
U5 = sum(Ud(d) for d in range(5, Lm, 10))
V = sum(lambert(2 * d, Lm) for d in range(1, Lm))
V5 = sum(lambert(10 * d, Lm) for d in range(1, Lm))
Pj = {}
for j in range(1, 11):
    s = np.zeros(Lm, dtype=np.int64)
    for r in range(1, Lm, 2):
        s[j * r::10 * r] += 1
    Pj[j] = s
Pj_alt = {j: sum(Ud(d) for d in range(1, Lm) if d % 10 == j % 10) for j in range(1, 11)}
R.check(all((Pj[j] == Pj_alt[j]).all() for j in range(1, 11)), "P_j = sum_(d = j mod 10) U_d (alternative definition)")
Ta = {a: sum(lambert(2 * d, Lm) for d in range(1, Lm) if d % 10 == a % 10) for a in range(1, 11)}
R.check(all((Pj[2 * j] == Ta[j] + Ta[j + 5] - Ta[2 * j if 2 * j <= 10 else 2 * j - 10]).all() for j in range(1, 6)),
        "P_(2j) = T_j + T_(j+5) - T_(2j) for 1 <= j <= 5")
conv = lambda x, y: np.convolve(x, y)[:Lm]
F5 = sum(conv(Pj[a], Ta[10 - a]) for a in (1, 3, 5, 7, 9))
G5 = sum(conv(Pj[j], Pj[j + 5]) for j in range(1, 6))
Aodd = Aser.copy()
Aodd[0::2] = 0
rhs = conv(U5, V) + conv(U, V5) - F5 - G5
R.check((Aodd == rhs).all(), f"A_odd = U5 V + U V5 - F5 - G5 exactly (n < {Lm})")
Rser = conv(U5, V) + conv(U, V5)

R.section("Proposition 6.3")
W2 = {1: {3, 6, 8, 9}, 3: {4, 7, 8, 9}, 7: {1, 2, 3, 6}, 9: {1, 2, 4, 7}}


def w(m, s):
    return int(s % 10 in W2.get(m % 10, set()))


delta = lambda z: int(z % 5 == 0)
Delta = lambda a, b: delta(a - b) - delta(a + b)
ok = all(2 * w(m, s) == int((m * s) % 5 != 0) + Delta(m, 2 * s) + (-1) ** s * Delta(m, s)
         for m in range(1, 41, 2) for s in range(1, 41))
R.check(ok, "Table 2: 2w(m,s) = 1[5 does not divide ms] + Delta(m,2s) + (-1)^s Delta(m,s)")
# F5 + G5 mod 2 equals the weighted quadruple count
# direct: coefficient at N of sum over m,r odd, s,t >= 1 of w(m,s) q^{mr + 2st}
Wd = np.zeros(Lm, dtype=np.int64)
for ma, sset in W2.items():
    A_ = np.zeros(Lm, dtype=np.int64)
    for m in range(ma, Lm, 10):
        A_[m::2 * m] += 1                      # q^{m r}, r odd
    for sb in sset:
        B_ = np.zeros(Lm, dtype=np.int64)
        for s_ in range(sb if sb else 10, Lm, 10):
            B_[2 * s_::2 * s_] += 1            # q^{2 s t}, t >= 1
        Wd += conv(A_, B_)
R.check(((F5 + G5 - Wd) % 2 == 0).all(), "F5 + G5 = generating function of w(m,s) q^(mr+2st) (mod 2)")

d5arr = np.zeros(Lm, dtype=np.int64)
for d in range(1, Lm):
    if d % 5:
        d5arr[d::d] += 1


def K5(n):
    return sum(int(d5arr[n - 2 * t] * d5arr[t]) for t in range(1, (n - 1) // 2 + 1))


def q0(d):
    return (d - 1) // 2 + 2 * (d // 5) - 3 * (d // 10)


def Qf(n):
    return sum(q0(d) for d in N.divisors(n))


def sigma(n):
    return sum(N.divisors(n))


def Jf(n):
    return (sigma(n) - tau(n)) // 2


def Hj(n, j):
    return sum(d // j for d in N.divisors(n))


ok_even = ok_11 = ok_W = ok_Q = True
for n in range(1, Lm, 2):
    if n % 5 == 0:
        continue
    k5, qq = K5(n), Qf(n)
    ok_even &= (k5 - qq) % 2 == 0
    ok_11 &= (int(Aser[n]) - int(Rser[n]) - (k5 - qq) // 2) % 2 == 0
    W = 0
    for m in range(1, n, 2):
        for r in range(1, n, 2):
            rest = n - m * r
            if rest <= 0:
                break
            if rest % 2:
                continue
            st = rest // 2
            for s in N.divisors(st):
                W += w(m, s)
    ok_W &= 2 * W == k5 - qq
    ok_Q &= qq == Jf(n) + 2 * Hj(n, 5) - 3 * Hj(n, 10)
    if n > 400:
        break
R.check(ok_even, "K5(N) - Q(N) is even (odd N coprime to 5, N < 2000)")
R.check(ok_11, "(11): A_N = R(N) + (K5(N) - Q(N))/2 (mod 2)")
R.check(ok_W, "integer identity 2W(N) = K5(N) - Q(N) (odd N < 400 coprime to 5)")
R.check(ok_Q, "Q(N) = J(N) + 2H_5(N) - 3H_10(N), J = (sigma - tau)/2")


def Ssum(n, kind):
    tot = 0
    for a in range(1, n + 1):
        for x in range(1, n // a + 1):
            rest = n - a * x
            if rest <= 0:
                break
            for b in N.divisors(rest):
                y = rest // b
                if kind == 6 and not (a % 2 and x % 2 and y % 2):
                    continue
                if kind == 8 and not (a % 2 and x % 2):
                    continue
                if kind == 10 and not (a % 2 and b % 2 and x % 2):
                    continue
                tot += Delta(a, b)
    return tot


ok = all(Ssum(n, 6) == 0 and Ssum(n, 8) == Jf(n) - Hj(n, 5) and Ssum(n, 10) == Jf(n) - Hj(n, 10)
         for n in range(1, 160, 2) if n % 5)
R.check(ok, "S6 = 0, S8 = J - H5, S10 = J - H10 (Liouville-type identities, odd N < 160 coprime to 5)")

R.section("Characters of Section 7 and the congruence (12)")
log2 = {pow(2, j, P): j for j in range(28)}
R.check(len(log2) == 28, "2 has order 28 modulo 29 (so lambda(2) = i determines lambda)")
lam = lambda a: 0 if a % P == 0 else [1, 1j, -1, -1j][log2[a % P] % 4]
R.check(lam(-1) == -1 and lam(5) == -1 and all(abs(lam(a) ** 2 - N.legendre(a, P)) < 1e-12 for a in range(1, 29)),
        "lambda(-1) = -1, lambda(5) = -1, lambda^2 = (./29)")
ok = True
for n in range(1, 30000):
    if N.legendre(n, P) != -1:
        continue
    dv = N.divisors(n)
    d5n = sum(1 for d in dv if d % 5)
    Ecoef = sum(lam(d) for d in dv if d % 5).real
    ok &= d5n % 2 == 0 and (d5n // 2 - round(Ecoef)) % 2 == 0
R.check(ok, "(12): if (n/29) = -1 then d5(n) is even and d5(n)/2 = [q^n]E (mod 2), n < 30000")

R.section("Theorem 7.1 (FFT, odd N divisible by 29, N < 2.1e6)")
NM = 2_100_000
d5big = np.zeros(NM + 1, dtype=np.int64)
nu5big = np.zeros(NM + 1, dtype=np.int64)
for d in range(1, NM + 1):
    if d % 5:
        d5big[d::d] += 1
    if d % 5 in (2, 3):
        nu5big[d::d] += 1
Dq2 = np.zeros(NM + 1, dtype=np.float64)
Dq2[0::2] = d5big[: NM // 2 + 1][: len(Dq2[0::2])]
Dq2[0] = 0
n = 1
while n < 2 * (NM + 1):
    n *= 2
K = np.rint(np.fft.irfft(np.fft.rfft(d5big.astype(np.float64), n) * np.fft.rfft(Dq2, n), n)[:NM + 1]).astype(np.int64)
R.check(all(K[m] == K5(m) for m in range(1, 400, 2)), "FFT convolution agrees with the direct K5 for odd N < 400")
ok = True
cnt = 0
for Nn in range(P, NM + 1, 2 * P):
    m = Nn
    while m % P == 0:
        m //= P
    prev = int(K[Nn // (P * P)]) if Nn % (P * P) == 0 else 0
    ok &= (int(K[Nn]) - prev - 2 * (int(nu5big[m]) + int(d5big[m]))) % 4 == 0
    cnt += 1
R.check(ok, f"K5(N) - K5(N/29^2) = 2(nu5(m) + d5(m)) mod 4 for all {cnt} odd N divisible by 29 below {NM}")

R.section("Proposition 8.1 ingredients")
ok = True
for n in range(1, Lm, 2):
    if n % 5 == 0:
        continue
    rep = sum(1 for x in range(1, math.isqrt(n) + 1, 2) for y in range(1, math.isqrt(n) + 1)
              if 5 * x * x + 2 * y * y == n or x * x + 10 * y * y == n)
    ok &= (rep - int(Rser[n])) % 2 == 0
R.check(ok, "R(N) = #{N = 5x^2 + 2y^2 or x^2 + 10y^2, x odd, y >= 1} (mod 2), odd N < 2000 coprime to 5")
R.check(all((q0(d) - q0(d + 40)) % 4 == 0 for d in range(1, 2000) if d % 2 and d % 5), "q0(d) mod 4 is periodic modulo 40")
R.check(all((q0(d) + q0(29 * d)) % 4 == 2 for d in range(1, 40) if d % 2 and d % 5), "q0(d) + q0(29d) = 2 (mod 4) on the 16 unit classes mod 40")
ok = True
for Nn in range(P, 200000, 2 * P):
    if Nn % 5 == 0:
        continue
    m = Nn
    while m % P == 0:
        m //= P
    prev = Qf(Nn // (P * P)) if Nn % (P * P) == 0 else 0
    ok &= (Qf(Nn) - prev - 2 * tau(m)) % 4 == 0
R.check(ok, "Q(N) - Q(N/29^2) = 2 tau(m) (mod 4) for odd N < 2e5 divisible by 29, coprime to 5")

R.section("Theorem 8.2")
spf = N.spf_sieve(2_000_000)
ok = True
cnt = 0
for r in range(1, 346, 2):
    if r % 5 == 0:
        continue
    g1 = N.cone_count((P * r) ** 2, spf)
    g0 = N.cone_count(r * r, spf)
    expect = 1 if r % 5 in (2, 3) else 0
    ok &= (g1 + g0) % 2 == expect
    cnt += 1
R.check(ok, f"gamma(29r) + gamma(r) = 1[r = 2,3 mod 5] (mod 2) for all {cnt} odd r <= 345 coprime to 5")

R.section("Index and Sturm bound")
base = 160 * P
for q in (2, 5, P):
    base = base * (q + 1) // q
R.check(base == 8640 and base * 32 == 276480 and 2 * 276480 // 12 == 46080 and 29 * 46080 == 1336320,
        "[SL2(Z):Gamma_0(4640)] = 8640, index 276480, Sturm bound 46080, expansion length 1336320")

sys.exit(0 if R.done() else 1)
