"""Theorems 1.1 and 1.2 and Corollary 5.1, checked directly on b(n), b1(n).

  * the generating function of Section 2 (Theorems 1.2/4.6 of Dixit-Kumar-
    Srivastava) agrees with a brute-force count from the definition, n < 90;
  * b(n), b1(n) mod 64 are then computed from the generating function for
    n < 62000;
  * Theorem 1.1 for every prime p = 21, 29 (mod 40) below 400 and every n in
    range with v_p(40n+1) (resp. v_p(40n+9)) odd; the two progressions;
  * Theorem 1.2 (both congruences, all n with 841n + 189 < 62000);
  * Corollary 5.1 for p = 29, 61, 101, 109 at nonsquare discriminants; and,
    for comparison, whether the scaling also holds at square discriminants.
"""
import sys

import nrlib as N

R = N.Report("check_partitions")
L = 62000
R.section("Generating function vs definition")
B = {l: N.B_series(l, L, 64) for l in (0, 1)}
for l in (0, 1):
    br = N.b_bruteforce(89, l)
    R.check(all((br[n] - int(B[l][n])) % 64 == 0 for n in range(90)),
            f"b_{l}(n) from the definition equals the series coefficient (mod 64) for n < 90")
R.check(N.b_bruteforce(10, 0)[:6] == [1, 0, 1, 1, 2, 2] and B[0][0] == 1 and B[1][0] == 1, "b(0) = b1(0) = 1")
b = [int(x) % 4 for x in B[0]]
b1 = [int(x) % 4 for x in B[1]]

R.section("Theorem 1.1")
primes = [p for p in range(3, 400) if N.is_prime(p) and p % 40 in (21, 29)]
for p in primes:
    cnt0 = cnt1 = 0
    ok = True
    for n in range(L):
        if N.vp(40 * n + 1, p) % 2 == 1:
            ok &= b[n] == 0
            cnt0 += 1
        if N.vp(40 * n + 9, p) % 2 == 1:
            ok &= b1[n] == 0
            cnt1 += 1
    R.check(ok and cnt0 > 0, f"p = {p}: b(n) = 0 mod 4 whenever v_p(40n+1) odd ({cnt0} n), b1(n) = 0 mod 4 whenever v_p(40n+9) odd ({cnt1} n)")
ok = all(b[29 * k + 21] == 0 for k in range(1, (L - 21) // 29) if k % 29)
R.check(ok, f"b(29k + 21) = 0 mod 4 for 1 <= k < {(L - 21) // 29}, 29 does not divide k (Dixit-Kumar-Srivastava conjecture)")
ok = all(b1[29 * k + 15] == 0 for k in range(0, (L - 15) // 29) if k % 29 != 6)
R.check(ok, "b1(29k + 15) = 0 mod 4 for k != 6 (mod 29) in range")
R.info(f"b(21) mod 4 = {b[21]}  (k = 0, valuation 2: not covered)")

R.section("Theorem 1.2")
n0 = (L - 1 - 189) // 841
R.check(all((b[841 * n + 21] + b[n]) % 4 == 0 for n in range(n0 + 1)), f"b(841n + 21) + b(n) = 0 mod 4 for 0 <= n <= {n0}")
R.check(all((b1[841 * n + 189] + b1[n]) % 4 == 0 for n in range(n0 + 1)), f"b1(841n + 189) + b1(n) = 0 mod 4 for 0 <= n <= {n0}")

R.section("Corollary 5.1")
for p in (29, 61, 101, 109):
    sp = 1 if p % 10 == 1 else -1
    ok_ns = True
    sq_fail = 0
    sq_tot = 0
    for l, arr in ((0, b), (1, b1)):
        c = (2 * l + 1) ** 2
        shift = c * (p * p - 1) // 40
        for n in range(L):
            m = p * p * n + shift
            if m >= L:
                break
            holds = (arr[m] - sp * arr[n]) % 4 == 0
            if N.is_square(40 * n + c):
                sq_tot += 1
                sq_fail += not holds
            else:
                ok_ns &= holds
    R.check(ok_ns, f"p = {p}: b_l(p^2 n + (2l+1)^2 (p^2-1)/40) = {sp:+d} b_l(n) mod 4 at nonsquare 40n+(2l+1)^2")
    R.info(f"p = {p}: at square discriminants the same scaling fails in {sq_fail} of {sq_tot} cases in range")

sys.exit(0 if R.done() else 1)
