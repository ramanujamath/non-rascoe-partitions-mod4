"""Exact, independent verifier for the refined original-cone certificate.
Geometry uses Fourier--Motzkin elimination with strict and weak inequalities;
it does not use the polygon subdivision in the search program.
"""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import json, time
if not __debug__:
    raise SystemExit('Run without Python -O: this verifier uses assertions.')
I = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
S = ((1, 0, 0), (0, 0, 1), (0, 1, 0))
H = ((1, 0, 0), (0, Q(1, 2), 0), (0, 0, 2))
V = ((1, 0, 0), (0, Q(1, 5), 0), (0, 0, 5))
JM = [
    ((1, -2, 0), (2, -2, 1), (0, 1, 0)),
    ((-1, 2, 0), (0, 1, 0), (2, -2, 1)),
    ((3, -2, 2), (2, -1, 2), (-2, 2, -1)),
    ((-3, 2, -2), (-2, 2, -1), (2, -1, 2)),
    ((1, 0, 2), (0, 0, 1), (-2, 1, -2)),
]
L = [(1, -2, 0), (2, -2, 1), (3, -2, 2), (4, -2, 4)]

def mm(A, B):
    return tuple(
        tuple(sum(a * b for a, b in zip(row, col)) for col in zip(*B))
        for row in A
    )

def linear(l, A):
    return tuple((sum((a * b for a, b in zip(l, col))) for col in zip(*A)))

def feasible(cond):
    lower = []
    upper = []
    last = []
    for a, b, c, strict in cond:
        if a > 0:
            lower.append((-Q(b, a), -Q(c, a), strict))
        elif a < 0:
            upper.append((-Q(b, a), -Q(c, a), strict))
        else:
            last.append((Q(b), Q(c), strict))
    for b, c, s in lower:
        for d, e, t in upper:
            last.append((d - b, e - c, s or t))
    lo = hi = None
    ls = hs = False
    for b, c, s in last:
        if b > 0:
            q = -c / b
            if lo is None or q > lo:
                lo, ls = (q, s)
            elif q == lo:
                ls = ls or s
        elif b < 0:
            q = -c / b



            if hi is None or q < hi:
                hi, hs = (q, s)
            elif q == hi:
                hs = hs or s
        elif c < 0 or (c == 0 and s):
            return False
    return lo is None or hi is None or lo < hi or (lo == hi and (not (ls or hs)))

def constraints(A, b):

    def row(k, sign, strict):
        c, a, z = linear(L[k], A)
        return (sign * a, sign * z, sign * c, strict)
    if b >= 5:
        return [row(b - 5, 1, False), row(b - 5, -1, False)]
    out = []
    if b > 0:
        out.append(row(b - 1, -1, True))
    if b < 4:
        out.append(row(b, 1, True))
    return out

def geometry():
    out = []
    for b in range(5):
        maps = [H, JM[b], mm(H, JM[b]), S, mm(H, S), H, mm(H, H), V, mm(H, V)]
        initial = [(1, 0, 0, True), (0, 1, 0, True)] + constraints(I, b)

        def rec(i, cond, sig):
            if i == len(maps):
                out.append(tuple((tuple(sig[j:j + 2]) for j in range(0, 10, 2))))
                return
            for c in range(9):
                cs = cond + constraints(maps[i], c)
                if feasible(cs):
                    rec(i + 1, cs, sig + [c])
        rec(0, initial, [b])
    return sorted(out)

def J(t, b):
    return tuple((sum((a * v for a, v in zip(row, t))) for row in JM[b]))

def verify_transformations():
    norm = ((1, 0, 0), (0, 0, 1), (0, 1, 0))
    image = [4, 1, 2, 3, 0]
    for b, A in enumerate(JM):
        assert mm(mm(tuple(zip(*A)), norm), A) == norm
        assert mm(JM[image[b]], A) == I
        base = [(1, 0, 0, True), (0, 1, 0, True)] + constraints(I, b)
        for c, a, z in A:
            assert not feasible(base + [(-a, -z, -c, False)])
        for a, z, c, strict in constraints(A, image[b]):
            assert strict
            assert not feasible(base + [(-a, -z, -c, False)])
    for A in [S, H, V]:
        assert mm(mm(tuple(zip(*A)), norm), A) == norm
    print('Verified all five branch maps: norm, positivity, '



         'branch image, and inverse.', flush=True)

def eligible(t, b):
    x, y, z = t
    if b[0] >= 5 or b[1] in [6, 8]:
        return False
    if b[1] == 5 and (y - x) % 20:
        return False
    if b[1] == 7 and (y - 3 * x - 4 * z) % 20:
        return False
    return True

def verify(path=None):
    start = time.time()
    verify_transformations()
    G = geometry()
    assert len(G) == 263, len(G)
    print('Independently enumerated', len(G),
          'geometric signatures, including boundaries.', flush=True)
    file = (Path(path) if path else
            Path(__file__).with_name('original_cone_stratified.json'))
    chosen = json.loads(file.read_text())
    assert len({(f, tuple(b), tuple(t), g) for f, b, t, g in chosen}) == len(chosen)
    assert all(
        f in ['J', 'S', 'T', 'V'] and len(b) == 2
        and b[0] in range(5) and b[1] in list(range(5)) + [5, 7]
        and len(t) == 3 and all(a in range(20) for a in t)
        and g in [0, 1]
        for f, b, t, g in chosen
    )
    ones = {(f, tuple(b), tuple(t), g) for f, b, t, g in chosen}

    def val(f, b, t, g):
        return int((f, tuple(b), tuple((a % 20 for a in t)), g) in ones)
    count = 0
    boundaries = 0
    for x, y, z in product(range(20), range(200), range(20)):
        if (x * x + 2 * y * z) % 40 != 1:
            continue
        t = (x, y, z)
        sw = (x, z, y)
        pre2 = (x, y // 2, 2 * z)
        pre5 = (x, y // 5, 5 * z)
        for b, bj, bs, b2, b5 in G:
            jt = J(t, b[0])
            if not all((eligible(a, c) for a, c in [(t, b), (jt, bj), (sw, bs)])):
                continue
            if y % 2 == 0 and (not eligible(pre2, b2)):
                continue
            if y % 5 == 0 and (not eligible(pre5, b5)):
                continue
            for g in range(2):
                rhs = (val('J', b, t, g)
                      ^ val('J', bj, jt, g ^ [1, 0, 0, 1, 1][b[0]]))
                rhs ^= val('S', b, t, g) ^ val('S', bs, sw, g ^ 1)
                if z % 2 == 0:
                    rhs ^= val('T', b, t, g)



                if y % 2 == 0:
                    rhs ^= val('T', b2, pre2, g ^ 1)
                if z % 5 == 0:
                    rhs ^= val('V', b, t, g)
                if y % 5 == 0:
                    rhs ^= val('V', b5, pre5, g)
                target = int(b[0] == 0 and y % 5 == 0 and (z % 4 == 0))
                assert rhs == target, (t, (b, bj, bs, b2, b5), g, rhs, target)
                count += 1
                if any((v >= 5 for bb in [b, bj, bs, b2, b5] for v in bb)):
                    boundaries += 1
    print('Verified', count, 'pointwise identities;', boundaries,
          'boundary-signature checks;', len(ones), 'nonzero entries;',
          round(time.time() - start, 2), 'seconds.', flush=True)
    return (G, count, boundaries)
if __name__ == '__main__':
    verify()
