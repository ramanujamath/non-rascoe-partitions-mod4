"""Exact Sturm certificate for the 29-adic convolution recurrence.

Section 7 of non_rascoe_congruences.tex proves that U_29(F) is an integral
holomorphic modular form of weight 2 on a group of index 276480.
The Sturm bound is 46080. This script computes F modulo 2 through
29*46080 and checks every required coefficient of U_29(F).

All arithmetic is exact. Python integers encode polynomials over F_2:
XOR is addition and carry-free shift/XOR multiplication is multiplication.
No floating-point convolution, numerical modular forms, or external
packages are used. The accompanying mathematical proof is essential.
"""
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path
import json
if not __debug__:
    raise SystemExit('Run without Python -O: this verifier uses assertions.')
P = 29
WEIGHT = 2
LEVEL = 160 * P
CHARACTER_INDEX = 32
BASE_INDEX = 160 * (P + 1) * 9 // 5
GROUP_INDEX = BASE_INDEX * CHARACTER_INDEX
BOUND = WEIGHT * GROUP_INDEX // 12
N = P * BOUND
MASK = (1 << N + 1) - 1
assert (BASE_INDEX, GROUP_INDEX, BOUND, N) == (8640, 276480, 46080, 1336320)
log = {pow(2, j, P): j for j in range(28)}
assert len(log) == 28 and pow(5, 7, P) == P - 1



chi5_log = {1: 0, 2: 1, 3: 3, 4: 2}
character_image = set()
for d in range(1, 40 * P):
    if gcd(d, 40 * P) != 1:
        continue
    chi4_log = 0 if d % 4 == 1 else 2
    alpha = (chi4_log + log[d % P]) % 4
    beta = (chi4_log + chi5_log[d % 5]) % 4
    chi8_log = 0 if d % 8 in (1, 7) else 1
    character_image.add((alpha, beta, chi8_log))
assert len(character_image) == CHARACTER_INDEX

def bits(indices):
    packed = bytearray((N + 8) // 8)
    for n in indices:
        packed[n >> 3] ^= 1 << (n & 7)
    return int.from_bytes(packed, 'little')

def odd_theta(d):
    return bits((d * x * x for x in range(1, isqrt(N // d) + 1, 2)))

def positive_theta(d):
    return bits((d * x * x for x in range(1, isqrt(N // d) + 1)))

def multiply(a, b):
    """Truncated polynomial product over F_2, with no integer carries."""
    if a.bit_count() > b.bit_count():
        a, b = (b, a)
    out = 0
    while a:
        low = a & -a
        out ^= b << low.bit_length() - 1
        a ^= low
    return out & MASK
E = bytearray(N + 1)
quadratic_residues = {x * x % P for x in range(1, P)}
for d in range(1, N + 1):
    if d % P in quadratic_residues and d % 5:
        for n in range(d, N + 1, d):
             E[n] ^= 1
E2 = bits((2 * n for n in range(1, N // 2 + 1) if E[n] and n % P))
Eodd = bits((n for n in range(1, N + 1, 2) if E[n] and n % P))
Op = odd_theta(P) ^ odd_theta(5 * P)
Opp = odd_theta(P * P) ^ odd_theta(5 * P * P)
D = positive_theta(2) ^ positive_theta(10)
Dp = positive_theta(2 * P) ^ positive_theta(10 * P)
Dpp = positive_theta(2 * P * P) ^ positive_theta(10 * P * P)
M = (multiply(odd_theta(1) ^ odd_theta(5) ^ Opp, E2)
     ^ multiply(D ^ Dpp, Eodd)
     ^ multiply(Opp, Dp)
     ^ multiply(Op, Dpp))
nu = bytearray(BOUND + 1)
for d in range(1, BOUND + 1):
    if d % 5 in (2, 3):
        for n in range(d, BOUND + 1, d):
             nu[n] ^= 1
J = bits((P * n for n in range(1, N // P + 1, 2) if nu[n]))



J ^= bits((P * P * n for n in range(1, N // (P * P) + 1, 2) if nu[n]))
F = M ^ J ^ Op ^ Opp
packed_F = F.to_bytes((N + 8) // 8, 'little')

def coefficient(n):
    return packed_F[n >> 3] >> (n & 7) & 1
up_coefficients = [coefficient(P * n) for n in range(BOUND + 1)]
assert not any(up_coefficients), [
    n for n, a in enumerate(up_coefficients) if a
][:10]
assert all((coefficient(n) == 0 for n in range(0, N + 1, 2)))
certificate = {
    'prime': P,
    'weight': WEIGHT,
    'level_after_U': LEVEL,
    'character_group_order': CHARACTER_INDEX,
    'group_index': GROUP_INDEX,
    'sturm_bound': BOUND,
    'input_expansion_bound': N,
    'modulus': 2,
    'coefficient_order': 'n=0,1,...,46080 for U_29(F)',
    'coefficients_mod2': ''.join(map(str, up_coefficients)),
    'packed_F_encoding': (
        'little-endian bits; bit n is [q^n]F mod2; truncated after N'),
    'packed_F_sha256': sha256(packed_F).hexdigest(),
}
Path(__file__).with_name('square_scaling_29_certificate.json').write_text(
    json.dumps(certificate, indent=2) + '\n', encoding='utf-8'
)
print(f'Weight {WEIGHT}; group index {GROUP_INDEX}; Sturm bound {BOUND}.')
print(f'All {BOUND + 1} coefficients of U_29(F) through the bound vanish modulo 2.')
print(f'Input expansion through {N} uses exact bitwise arithmetic.')
print('SHA256 of complete packed F expansion:', certificate['packed_F_sha256'])
print('Together with Section 7, this proves the 29-adic convolution recurrence.')
print('Section 8 deduces the companion congruences at ALL square discriminants.')
print('Scope of this script: square discriminants only.')
print('For the separate nonsquare proof, see non_rascoe_congruences.pdf, '
      'Sections 3--5.')
