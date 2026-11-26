# Audit: *Congruences modulo four for non-Rascoe partitions*

Draft audited: `manuscript/original/Congruences_Modulo_Four_Non_Rascoe_Submission.pdf` (28 pp.).
Rewrite: `manuscript/main.tex` → `manuscript/main.pdf` (17 pp.).

## 1. Summary

**No mathematical errors were found.** Every theorem was checked numerically far beyond the paper's examples, and both computer-assisted steps were re-run:

- **The paper's own programs.** The three programs from Appendices A–C were extracted from the PDF. Their indentation, which PDF extraction damages, was repaired and nothing else was changed. They reproduce the paper's figures exactly:
  - the certificate decodes with the printed SHA-256 `6c1b022e…6d60` (5485 records);
  - the cone verifier reports 263 signatures (92 interior, 171 with a boundary label), 1 156 160 pointwise identities and 272 960 boundary checks;
  - the Sturm verifier reports index 276480 and bound 46080, and all 46081 coefficients of `U₂₉𝔉` vanish mod 2.
- **The cancellation identity (Lemma 4.1), checked independently.** The identity was evaluated directly on actual triples of `S_p(D)`, without the Fourier–Motzkin enumeration: 27 420 triples, for 36 nonsquare discriminants and `p = 29, 61, 101`. For each triple the sign bit was computed from the Legendre symbol itself, so this also checks Lemma 3.1 independently.
- **Theorem 7.1.** It was confirmed directly by FFT for every odd `N < 2.1·10⁶` divisible by 29 (36 207 values). This supports the modular-form step, whose level and character bookkeeping (Table 3) is asserted rather than derived in detail.

The prose was rewritten and made self-contained:
- notation is defined up front;
- the cited classical inputs are stated as Facts F1–F7;
- the jargon is replaced, and the heavily overloaded symbols are renamed (see §4).

## 2. Verification

Run `python3 verification/run_all.py` (standard library and numpy; about 25 s). Logs go to `verification/logs/`.

| Script | Result | What is checked |
|---|---|---|
| `appendix/decode_certificate.py` | SHA-256 matches | Appendix A, as printed in the paper. |
| `appendix/verify_original_certificate.py` | 1 156 160 identities | Appendix B: branch maps, 263 signatures and all pointwise identities of Lemma 4.1. |
| `appendix/prove_square_scaling_29.py` | all 46081 coefficients vanish | Appendix C: the Sturm computation for Theorem 7.1. |
| `check_partitions.py` | 21/21 | The generating function (Fact F1) against a brute-force count from the definition (n < 90), then `b`, `b₁` mod 64 for n < 62 000. Thm 1.1 for every prime p ≡ 21, 29 (mod 40) below 400, at every n in range. Both progressions, including the Dixit–Kumar–Srivastava conjecture for k < 2137. Thm 1.2 for n ≤ 73. Cor 5.1 for p = 29, 61, 101, 109. |
| `check_reduction.py` | 15/15 | `P∞/Q∞ ≡ 1 + 2L`; the Hecke–Rogers identity (7) against the definition of `H_ℓ`; `H_ℓ ≡ T_ℓ`; congruence (4); eq. (5); Prop 2.1; vanishing of `h_ℓ` and `U` at odd valuation; the Legendre-symbol facts. |
| `check_cone.py` | 12/12 | Thm 1.3 for every nonsquare D ≤ 4·10⁶ and p ∈ {29, 61, 101, 109, 149, 181}. Jackson's map on actual sets. Lemma 3.1 from actual Legendre symbols. Lemma 4.1 on actual triples. The marked count is even. |
| `check_divisor.py` | 22/22 | Prop 6.1 (m < 2001). Prop 6.2 as an exact series identity. Table 2 formula, `P_{2j}` identity, the Liouville evaluations `S₆, S₈, S₁₀`, `2W = K₅ − Q`, (11). The characters `λ` and (12) for n < 30000. Thm 7.1 by FFT. Prop 8.1 ingredients. Thm 8.2 for odd r ≤ 345. Index and Sturm numbers. |

## 3. Mathematical review

1. **No errors found** in the reduction (§2), the involution and sign (§3), the certificate argument (§§4–5), the divisor transform (§6), or the scaling deduction (§8).
2. **Corollary 5.1 at square discriminants.** For p = 61 and 109 the scaling fails at squares, as the verification shows. So the nonsquare hypothesis is genuinely needed in general. For p = 29 it holds at squares, which is exactly what Thm 1.2 adds. The rewrite records this in a remark.
3. **Theorem 1.3 at square D.** It does fail at squares: for p = 29, 14 of the 27 square D ≤ 4·10⁶ violate it. That matches Thm 8.2, since `γ(29r) + γ(r)` is odd for r ≡ 2, 3 (mod 5). It is consistent with the stated hypothesis.
4. **Prop 2.1 proof.** The draft's displayed derivative formula was garbled in places. The rewrite writes it out:
   - where the factor `(1+q^j)` comes from;
   - how the `q^{2j}/(1−q^j)` term becomes the `u = j` term of `C₀`;
   - that `(−1)^r ≡ 1 + 2r (mod 4)` is the step behind congruence (4).
5. **Lemma 3.1 proof.** The rewrite gives the three new second coordinates with their correct signs: `−(x−2y)²/2y`, `−(x−y)²/y`, `(x−2y)²/2y`. The draft's display omitted a sign; the conclusion was unaffected.
6. **Theorem 8.2 proof.** The rewrite states that both Möbius-inversion formulas are *exact* identities (the second telescopes), not merely congruences.
7. **Modular-form step (§7.4).** The characters and levels in Table 3, and the claim that odd projection raises the level to 4L, are stated without a detailed derivation.
   - These were not re-derived symbolically here.
   - The direct numerical confirmation of Theorem 7.1 up to 2.1·10⁶ is strong evidence that the argument is sound.
   - A referee may still ask for the character computations of Table 3 to be written out.
8. **Python version.** The Sturm verifier uses `int.bit_count`, which needs Python 3.10 or later. The paper says so.

### Points for the author
- **References not verified:** Cuthbertson (arXiv:2608.30180) and the exact theorem numbers in Alaca–Alaca–McAfee–Williams, Sills, Andrews and Stein.
- **The Liouville-identity specializations** (eq. (7.2) in the rewrite) were confirmed numerically for odd N < 160 prime to 5. They were not re-derived from the cited theorems.

## 4. Presentation changes

- **New §2 "Notation and standard facts."**
  - It defines: `[q^n]`, the q-Pochhammer symbols, odd part, `τ, σ, μ, d₅, ν₅`, Dirichlet convolution, the Legendre facts for `p ≡ 21, 29 (mod 40)`, `Γ₀(L)`, `M_k` and characters.
  - It states the outside inputs as Facts F1–F7: Dixit–Kumar–Srivastava, Sills, Hecke–Rogers, the Liouville identities, the weight-one Eisenstein series, theta series and Sturm.
- **Each technical section opens with a plain-language explanation**, for example what the certificate is and why a finite check suffices.
- **Symbol renames.** The draft reused many letters:

| Draft | Rewrite | Clash |
|---|---|---|
| `T_ℓ(q)` (theta series) | `Θ_ℓ` | partial bijections `T_ℓ`, Lambert `T_a` |
| `M_ℓ(q)` | `N_ℓ` | `M_k(Γ)`, `M_p` |
| `L` (Lambert) | `Λ` | `L_i(t)`, `L(s,ρ)`, `L_d` |
| `S_p(D)`, swap `S` | `𝒯_p(D)`, `𝒮` | `S_d`, `S₆…` |
| Jackson `J` | `𝒥` | `J(N)`, series `J` |
| partial bijections `T_ℓ` | `A_ℓ^{-1}` (one family `A_ℓ`) | — |
| label `η(t)` | `lab(t)` | — |
| bits `δ_β` | `c_β` | `δ(z)` |
| `ε_a` (exponent) | `e_a` | sign `ϵ(t)` |
| Section-6 `L_d, U_d, U, U₅, V, P_j, T_a, F₅, G₅` | `𝓛_d, 𝒪_d, 𝒪, 𝒪₅, 𝒱, 𝒫_j, 𝒯_a, ℱ₅, 𝒢₅` | `U(q)`, `U_p`, `P∞`, Θ |
| `H_j(N)` | `Φ_j(N)` | `H_ℓ` |
| quartic char. mod 5 `χ` | `ψ` | Legendre `χ` in §3 |
| characters `α, β, ξ` | `κ₁, κ₂, κ₃` | exponents `α_a`, branch `β` |
| `P_d, O_d, S_d, D_d` (theta) | `θ⁺_d, θ^odd_d, Ω_d, ϑ_d` | `P∞`, Lambert `P_j`, `D` |
| `D(q), D₀, Z, R, M_p, J, Y, F` | `𝒟, 𝒟₀, 𝒵, ℛ, ℳ, 𝒩, 𝒴, 𝔉` | various |

- **Layout.** The program listings were removed from the appendices. The appendices now describe the certificate encoding and the three programs, and give GitHub links to them in `verification/appendix/` (the SHA-256 of the decoded certificate is printed). There are no overfull lines.

## 5. Draft → rewrite map

| Draft | Rewrite |
|---|---|
| §1 (Thms 1.1–1.3) | §1 (Thms 1.1–1.3) |
| — | §2 Notation and standard facts (new) |
| §2 (Prop 2.1) | §3 (Prop 3.1) |
| §3 (Table 1, Lemma 3.1) | §4 (Table 1, Lemma 4.1) |
| §4 (Lemma 4.1) | §5 (Lemma 5.1) |
| §5 (proofs, Cor 5.1) | §6 (proofs, Cor 6.1, new Remark 6.2) |
| §6 (Props 6.1–6.3, Table 2) | §7 (Props 7.1–7.3, Table 2) |
| §7 (Thm 7.1, Table 3) | §8 (Thm 8.1, Table 3) |
| §8 (Prop 8.1, Thm 8.2, proof of Thm 1.2) | §9 (Prop 9.1, Thm 9.2) |
| §9 | §10 |
| Appendices A–C (code listings) | Appendices A–C (descriptions + GitHub links) |
