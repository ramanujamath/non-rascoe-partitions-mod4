"""Sections 3-5: the cone count C(D) and the cancellation argument.

  * Theorem 1.3 for every nonsquare D <= 4*10^6, D = 1, 9 (mod 40), p | D,
    p in {29, 61, 101, 109, 149, 181}; also how often the parity relation
    fails at square D (it is not claimed there);
  * Jackson's map J (Table 1) on the actual finite sets S_p(D): it preserves
    x^2 + 2yz, stays positive, is an involution, maps branch b to
    (4,1,2,3,0)[b], and no triple lies on a branch boundary;
  * Lemma 3.1: the bit g defined from the Legendre symbol changes as stated
    under S, J (by (1,0,0,1,1)[branch]), T_2 and T_5 -- computed from the
    actual triples, not from the formulas of the proof;
  * Lemma 4.1 evaluated directly on every triple of S_p(D) (and on 3t for
    D = 9 mod 40), using the decoded certificate but not the
    Fourier--Motzkin enumeration;
  * the summation step: the number of marked triples is even.
"""
import json
import math
import os
import subprocess
import sys

import nrlib as N

R = N.Report("check_cone")
HERE = os.path.dirname(os.path.abspath(__file__))
PRIMES = [29, 61, 101, 109, 149, 181]

R.section("Theorem 1.3 on all D <= 4*10^6")
Dmax = 4 * 10 ** 6
C = N.cone_count_sieve(Dmax)
for p in PRIMES:
    ok = True
    cnt = 0
    sqfail = sqtot = 0
    for D in range(p, Dmax + 1, p):
        if D % 40 not in (1, 9):
            continue
        rel = (int(C[D]) - (int(C[D // (p * p)]) if D % (p * p) == 0 else 0)) % 2 == 0
        if N.is_square(D):
            sqtot += 1
            sqfail += not rel
        else:
            ok &= rel
            cnt += 1
    R.check(ok, f"p = {p}: C(D) = C(D/p^2) mod 2 for all {cnt} nonsquare D <= 4e6, D = 1,9 mod 40, p | D")
    R.info(f"p = {p}: at square D the relation fails in {sqfail} of {sqtot} cases")

R.section("Jackson's involution and the sign bit on actual sets S_p(D)")
JM = [((1, -2, 0), (2, -2, 1), (0, 1, 0)), ((-1, 2, 0), (0, 1, 0), (2, -2, 1)),
      ((3, -2, 2), (2, -1, 2), (-2, 2, -1)), ((-3, 2, -2), (-2, 2, -1), (2, -1, 2)),
      ((1, 0, 2), (0, 0, 1), (-2, 1, -2))]


def Lvals(t):
    x, y, z = t
    return (x - 2 * y, 2 * x - 2 * y + z, 3 * x - 2 * y + 2 * z, 4 * x - 2 * y + 4 * z)


def beta(t):
    """branch / refined label of a positive rational triple (Section 4)."""
    Ls = Lvals(t)
    for i, v in enumerate(Ls):
        if v == 0:
            return 5 + i
    return sum(1 for v in Ls if v < 0)


def J(t):
    bb = beta(t)
    return tuple(sum(a * v for a, v in zip(row, t)) for row in JM[bb]), bb


def S_set(D, p):
    out = []
    x = 1
    while x * x < D:
        M = D - x * x
        if M % 2 == 0:
            M //= 2
            for y in N.divisors(M):
                z = M // y
                if not (x % p == 0 and y % p == 0 and z % p == 0):
                    out.append((x, y, z))
        x += 1
    return out


def gbit(t, p):
    x, y, z = t
    eps = N.legendre(y, p) if y % p else N.legendre(-2 * z, p)
    assert eps in (1, -1)
    return 0 if eps == 1 else 1


Ds = []
for p in (29, 61, 101):
    k = 0
    D = p
    while len([d for d in Ds if d[1] == p]) < 12:
        if D % 40 in (1, 9) and not N.is_square(D):
            Ds.append((D, p))
        D += p
okJ = okBits = okBoundary = True
ntrip = 0
for D, p in Ds:
    St = S_set(D, p)
    Sset = set(St)
    ntrip += len(St)
    for t in St:
        bb = beta(t)
        okBoundary &= bb < 5
        jt, _ = J(t)
        okJ &= jt in Sset and jt[0] ** 2 + 2 * jt[1] * jt[2] == D and min(jt) > 0
        jjt, bj = J(jt)
        okJ &= jjt == t and bj == (4, 1, 2, 3, 0)[bb]
        g = gbit(t, p)
        okBits &= gbit(jt, p) == (g + (1, 0, 0, 1, 1)[bb]) % 2
        okBits &= gbit((t[0], t[2], t[1]), p) == (g + 1) % 2
        if t[2] % 2 == 0:
            okBits &= gbit((t[0], 2 * t[1], t[2] // 2), p) == (g + 1) % 2
        if t[2] % 5 == 0:
            okBits &= gbit((t[0], 5 * t[1], t[2] // 5), p) == g
R.check(okBoundary, f"no triple of S_p(D) lies on a branch boundary ({len(Ds)} nonsquare D, {ntrip} triples)")
R.check(okJ, "J preserves x^2 + 2yz and positivity, is an involution on S_p(D), and maps branch b to (4,1,2,3,0)[b]")
R.check(okBits, "Lemma 3.1: g(St) = g+1, g(Jt) = g + (1,0,0,1,1)[branch], g(T_2 t) = g+1, g(T_5 t) = g on actual triples")

R.section("Lemma 4.1 on actual triples (decoded certificate)")
cert_path = os.path.join(HERE, "appendix", "original_cone_stratified.json")
if not os.path.exists(cert_path):
    subprocess.run([sys.executable, "decode_certificate.py"], cwd=os.path.join(HERE, "appendix"), check=True)
records = json.load(open(cert_path))
ones = {(f, tuple(b), tuple(t), g) for f, b, t, g in records}
R.check(len(records) == 5485 and len(ones) == 5485, "certificate: 5485 distinct records")


def eta(t):
    x, y, z = t
    from fractions import Fraction as Fr
    return (beta(t), beta((Fr(x), Fr(y, 2), Fr(2 * z))))


def Fval(tag, t, g):
    return int((tag, eta(t), tuple(a % 20 for a in t), g) in ones)


ok = True
okEven = True
nchk = 0
for D, p in Ds:
    for t0 in S_set(D, p):
        t = t0 if D % 40 == 1 else tuple(3 * a for a in t0)
        g = gbit(t, p)
        x, y, z = t
        jt, _ = J(t)
        lhs = int(2 * y < x and y % 5 == 0 and z % 4 == 0)
        rhs = Fval("J", t, g) ^ Fval("J", jt, gbit(jt, p)) ^ Fval("S", t, g) ^ Fval("S", (x, z, y), gbit((x, z, y), p))
        if z % 2 == 0:
            rhs ^= Fval("T", t, g)
        if y % 2 == 0:
            a2 = (x, y // 2, 2 * z)
            rhs ^= Fval("T", a2, gbit(a2, p))
        if z % 5 == 0:
            rhs ^= Fval("V", t, g)
        if y % 5 == 0:
            a5 = (x, y // 5, 5 * z)
            rhs ^= Fval("V", a5, gbit(a5, p))
        ok &= lhs == rhs
        nchk += 1
    marked = sum(1 for (x, y, z) in S_set(D, p) if 2 * y < x and y % 5 == 0 and z % 4 == 0)
    okEven &= marked % 2 == 0
R.check(ok, f"identity (9) holds at every actual triple ({nchk} triples; D = 9 mod 40 handled via 3t)")
R.check(okEven, "the number of marked triples (2y < x, 5|y, 4|z) in S_p(D) is even for every tested D")

sys.exit(0 if R.done() else 1)
