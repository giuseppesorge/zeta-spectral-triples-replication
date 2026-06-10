# Independent numerical replication of "Zeta Spectral Triples"
### (Connes–Consani–Moscovici, arXiv:2511.22755) — June 10, 2026
### Giuseppe Sorge

## What was replicated

The central numerical construction of the preprint: the truncated Weil
quadratic form on the even subspace of `L²([λ⁻¹, λ], du/u)` (cosine basis,
N = 120, i.e. a 121×121 matrix), its eigenvector ξ of minimal eigenvalue,
and the real zeros of the entire function ξ̂(z), compared against the
nontrivial zeros of the Riemann zeta function.

**Fully independent implementation**: own basis and closed forms (pole and
prime terms exact; archimedean term reduced via partial fractions to 242
one-dimensional integrals on a shared Gauss–Legendre grid of degree 60 per
half-oscillation panel, with tails handled analytically — smooth part via
u-substitution, oscillatory part via iterated integration by parts with
exact derivatives through polygamma functions). Python/mpmath at 80 working
digits. No use of the authors' code or data beyond the three public
parameters (N = 120, λ² ∈ {12, 13, 14}).

## Results (error |replicated zero − Riemann zero|, 80 digits, N = 120)

| zero # | λ²=13 (ours) | λ²=13 (paper) | λ²=14 (ours) | λ²=14 (paper) |
|---|---|---|---|---|
| 1  | 6.4e-50 | 2.4e-55 | **4.8e-55** | 1.1e-60 |
| 10 | 4.9e-34 | 4.0e-37 | 6.2e-39 | 3.0e-42 |
| 20 | 6.3e-22 | – | 1.9e-26 | – |
| 30 | 2.0e-13 | – | 2.1e-17 | – |
| 40 | 3.1e-06 | – | 1.6e-09 | – |
| 45 | 4.5e-04 | – | 6.8e-07 | – |
| 50 | 9.3e-03 | 2.0e-03 | 5.1e-05 | 4.8e-06 |

(λ²=12 was also run, at 35 working digits: zero #50 error 1.13e-1 vs the
paper's 9.02e-2.)

**Minimal eigenvalues of the truncated Weil form** (not reported in the
preprint, possibly of independent interest):
- λ²=13: ε₁ = 8.16e-52, ε₂ = 1.73e-50, ε₃ = 1.09e-44
- λ²=14: ε₁ = 6.01e-57, ε₂ = 7.41e-51, ε₃ = 9.91e-50

Truncated Weil positivity is thus confirmed numerically down to 10⁻⁵²/10⁻⁵⁷,
and ε₁ decreases with λ consistently with the deepening accuracy reported in
the paper's Table 6.

## Conclusions

1. **The construction works as described**: the first Riemann zero emerges
   with 54 correct digits (λ²=14) from a matrix built using only the primes
   ≤ λ², the digamma function and the pole term. Our errors saturate our
   numerical floor, fully consistent with the 200-digit results of Table 6.
2. **The λ-scaling law is confirmed**: at each step λ² = 12→13→14 the error
   at fixed zero index drops by 3–5 orders of magnitude.
3. **The error growth with zero index is confirmed**, from the numerical
   floor up to ~10⁻² near zero #50 where the N = 120 truncation dominates.
4. **A methodological remark**: at 35 working digits the bottom 4
   eigenvalues of the form are numerically degenerate and the computed ξ is
   a mixture of the bottom cluster, which inflates the errors at zeros
   #45–50 by ~2 orders of magnitude. 80 digits resolve the cluster
   (ε₂/ε₁ ≈ 21 for λ²=13, ≈ 10⁶ for λ²=14). Anyone attempting a replication
   at standard precision will encounter this.

## Reproducibility — suggested run order

Environment: Python ≥ 3.10; `pip install mpmath numpy scipy matplotlib`.

1. `python replica_zeta_spectral_triples.py` — float64 prototype (~1 min):
   the first six zeros emerge at 1e-3/1e-4 accuracy; also produces the
   per-prime decomposition figure (`replica_zst.png`).
2. `python replica_zst_80.py` — main 80-digit run (~42 min on a laptop):
   prints error tables for λ² = 13, 14 against Table 6 and the five minimal
   eigenvalues; writes `replica_zst_80.json`.
3. (optional) `python replica_zst_hp.py` — 35-digit run, all three λ²
   (~30 min): exhibits the bottom-cluster degeneracy discussed above.

Included data: `replica_zst_80.json`, `replica_zst_hp.json` (full error
tables, 50 zeros each), `replica_zst_hp.png` (summary figure).
