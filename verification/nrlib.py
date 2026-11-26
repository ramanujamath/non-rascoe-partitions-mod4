"""Shared helpers (Python standard library + numpy).

* Report                      pass/fail bookkeeping
* truncated power series mod M (numpy; FFT products, exact for small M)
* b(n), b1(n) by brute force and by the generating function
* the Hecke--Rogers series H_l, the cone count C(D), the series U
* arithmetic helpers: Legendre symbol, divisors, smallest-prime-factor sieve
"""
import math
import sys

import numpy as np


class Report:
    def __init__(self, name):
        self.name = name
        self.npass = 0
        self.nfail = 0

    def check(self, cond, msg):
        if cond:
            self.npass += 1
            print(f"  [ok]   {msg}", flush=True)
        else:
            self.nfail += 1
            print(f"  [FAIL] {msg}", flush=True)
        return bool(cond)

    def info(self, msg):
        print(f"  [info] {msg}", flush=True)

    def section(self, msg):
        print(f"\n== {msg}", flush=True)

    def done(self):
        print(f"\n-- {self.name}: {self.npass} passed, {self.nfail} failed")
        return self.nfail == 0


# ------------------------------------------------------------ arithmetic

def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def vp(n, p):
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def is_square(n):
    return n >= 0 and math.isqrt(n) ** 2 == n


def spf_sieve(n):
    spf = np.zeros(n + 1, dtype=np.int64)
    for i in range(2, n + 1):
        if spf[i] == 0:
            spf[i::i][spf[i::i] == 0] = i
    return spf


def divisors_from_spf(n, spf):
    ds = [1]
    while n > 1:
        p = int(spf[n])
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        ds = [d * p ** k for d in ds for k in range(e + 1)]
    return ds


def divisors(n):
    small = [d for d in range(1, math.isqrt(n) + 1) if n % d == 0]
    return sorted(set(small + [n // d for d in small]))


# ------------------------------------------------------------ series mod M

def fft_mul(a, b, M, L):
    """product of two integer arrays modulo M (M small), truncated to length L.
    Exact: entries < M, so every convolution value is < M^2 * L < 2^53."""
    a = np.asarray(a[:L], dtype=np.int64) % M
    b = np.asarray(b[:L], dtype=np.int64) % M
    assert (M - 1) ** 2 * L < 2 ** 50
    n = 1
    while n < len(a) + len(b):
        n *= 2
    fa = np.fft.rfft(a.astype(np.float64), n)
    fb = np.fft.rfft(b.astype(np.float64), n)
    c = np.rint(np.fft.irfft(fa * fb, n)[:L]).astype(np.int64)
    out = np.zeros(L, dtype=np.int64)
    out[:len(c)] = c % M
    return out


def series_inv(a, M, L):
    """inverse of a series with constant term 1 modulo M (Newton iteration)."""
    a = np.asarray(a, dtype=np.int64) % M
    assert a[0] == 1
    inv = np.array([1], dtype=np.int64)
    n = 1
    while n < L:
        n = min(2 * n, L)
        t = fft_mul(a[:n], inv, M, n)           # a*inv
        t = (-t) % M
        t[0] = (t[0] + 2) % M                    # 2 - a*inv
        inv = fft_mul(inv, t, M, n)
    return inv[:L]


def euler_f(k, L, M):
    """f(q^k) = prod (1 - q^{kn}) modulo M, via the pentagonal number theorem."""
    out = np.zeros(L, dtype=np.int64)
    m = 0
    while True:
        added = False
        exps = (0,) if m == 0 else (m * (3 * m - 1) // 2, m * (3 * m + 1) // 2)
        for g in exps:
            if k * g < L:
                out[k * g] = (out[k * g] + (-1) ** m) % M
                added = True
        if not added and m > 0:
            break
        m += 1
    return out


def div_one_plus(x, r, M):
    """x / (1 + q^r) modulo M."""
    n = len(x)
    rows = -(-n // r)
    y = np.zeros(rows * r, dtype=np.int64)
    y[:n] = x
    y = y.reshape(rows, r)
    sgn = np.where(np.arange(rows) % 2 == 0, 1, -1)[:, None]
    y = (np.cumsum(y * sgn, axis=0) * sgn) % M
    return y.reshape(-1)[:n].copy()


def B_series(l, L, M=64):
    """B_l(q) = (-q;q)_inf sum_r (-1)^r q^{r(r+l)} / (-q;q)_r  modulo M."""
    S = np.zeros(L, dtype=np.int64)
    term = np.zeros(L, dtype=np.int64)
    term[0] = 1
    r = 0
    while True:
        S = (S + term) % M
        r += 1
        shift = 2 * r - 1 + l  # q^{r(r+l)} / q^{(r-1)(r-1+l)}
        if r * (r + l) >= L:
            break
        new = np.zeros(L, dtype=np.int64)
        new[shift:] = (-term[:L - shift]) % M
        term = div_one_plus(new, r, M)
    Pinf = fft_mul(euler_f(2, L, M), series_inv(euler_f(1, L, M), M, L), M, L)
    return fft_mul(Pinf, S, M, L)


def b_bruteforce(nmax, l):
    """b_l(n) for n <= nmax straight from the definition (distinct parts, the
    number of parts k is such that k + l is not a part)."""
    counts = [0] * (nmax + 1)

    # enumerate by total: generate all distinct partitions with sum <= nmax
    def gen(start_sum, maxpart, parts):
        k = len(parts)
        if (k + l) not in parts:
            counts[start_sum] += 1
        for p in range(min(maxpart, nmax - start_sum), 0, -1):
            parts.add(p)
            gen(start_sum + p, p - 1, parts)
            parts.discard(p)

    gen(0, nmax, set())
    return counts


def H_hecke_rogers(l, L):
    """coefficients of (7): sum sgn(x) (-1)^j q^{(x^2 - 40 j^2 - (2l+1)^2)/40}."""
    out = [0] * L
    c = (2 * l + 1) ** 2
    X = math.isqrt(40 * L + c) + 20
    # x^2 - 40 j^2 with 10|j| < |x| is at least x^2 - 0.4 x^2 = 0.6 x^2
    xmax = math.isqrt(int((40 * L + c) / 0.6)) + 20
    for x in range(-xmax, xmax + 1):
        if x % 10 != (2 * l + 1) % 10 or x == 0:
            continue
        s = 1 if x > 0 else -1
        jmax = (abs(x) - 1) // 10
        for j in range(-jmax, jmax + 1):
            num = x * x - 40 * j * j - c
            if num < 0 or num % 40:
                continue
            n = num // 40
            if n < L:
                out[n] += s * (-1) ** (j % 2)
    return out


def cone_count_sieve(Dmax):
    """array C[D] = #{(x,u,v) in Z_{>0}^3 : 10u < x, x^2 + 40uv = D}, D <= Dmax."""
    C = np.zeros(Dmax + 1, dtype=np.int32)
    x = 1
    while x * x + 40 < Dmax + 1:
        for u in range(1, (x - 1) // 10 + 1):
            start = x * x + 40 * u
            if start > Dmax:
                break
            C[start::40 * u] += 1
        x += 1
    return C


def cone_count(D, spf=None):
    """C(D) for a single D (divisor enumeration)."""
    tot = 0
    x = 1
    while x * x < D:
        M, rem = divmod(D - x * x, 40)
        if rem == 0 and M > 0:
            ds = divisors_from_spf(M, spf) if spf is not None and M < len(spf) else divisors(M)
            tot += sum(1 for u in ds if 10 * u < x)
        x += 1
    return tot


def U_series(L):
    """U(q) = sum_{a>=0} q^{eps_a}, eps_a = (a+1)(5a+2)/2."""
    out = [0] * L
    a = 0
    while (a + 1) * (5 * a + 2) // 2 < L:
        out[(a + 1) * (5 * a + 2) // 2] += 1
        a += 1
    return out
