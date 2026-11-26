# Non-Rascoe partitions modulo 4

Manuscript, audit and verification code for

> **Congruences modulo four for non-Rascoe partitions**
> Sumit Kumar Jha

The functions are defined as follows:
- `b(n)` counts partitions of `n` into distinct parts in which the number of parts is not a part;
- `b₁(n)` is defined the same way, but with "the number of parts plus one" not a part.

The paper proves the following.

- **Prime-valuation criterion.** Let `p ≡ 21, 29 (mod 40)` be prime. If `v_p(40n+1)` is odd, then `b(n) ≡ 0 (mod 4)`; likewise for `b₁` with `40n+9`. This includes the Dixit–Kumar–Srivastava conjecture `b(29k+21) ≡ 0 (mod 4)` for `29 ∤ k`.
- **Scaling relations.** `b(841n+21) + b(n) ≡ 0` and `b₁(841n+189) + b₁(n) ≡ 0 (mod 4)` for all `n ≥ 0`. The first is Cuthbertson's conjecture.

The proof has two computer-assisted steps:
- a finite cancellation certificate for nonsquare discriminants (Jackson's involution and a Legendre-symbol sign);
- a Sturm-bound check of a weight-two modular form for square discriminants.

## Layout

```
manuscript/
  main.tex, sections/*.tex, references.tex   rewritten manuscript
  code/                                      the three certificate programs (included verbatim)
  main.pdf                                   compiled rewrite
  original/                                  the audited draft
verification/
  appendix/decode_certificate.py             paper's Appendix A (certificate + SHA-256 check)
  appendix/verify_original_certificate.py    paper's Appendix B (cone-identity verifier)
  appendix/prove_square_scaling_29.py        paper's Appendix C (Sturm-bound verifier)
  nrlib.py                                   helpers for the independent checks
  check_partitions.py                        Thms 1.1, 1.2, Cor 5.1 directly on b(n), b1(n)
  check_reduction.py                         Section 2 (reduction to the cone count)
  check_cone.py                              Thm 1.3, Jackson's map, Lemma 3.1, Lemma 4.1 on actual triples
  check_divisor.py                           Sections 6-8 (divisor transform, Lambert series, Thm 7.1, Thm 8.2)
  run_all.py                                 runs everything, writes logs/ and logs/SUMMARY.txt
AUDIT.md                                     audit report
```

The appendix programs were extracted from the PDF. Their indentation, which PDF extraction damages, was repaired; the code is otherwise unchanged. They reproduce the paper's printed counts and checksum exactly.

## Running the checks

```
cd verification
python3 run_all.py
```

This needs Python 3.10 or later and numpy, and takes about 25 seconds. All checks pass.

## Building the manuscript

```
cd manuscript
latexmk -pdf main.tex
```
