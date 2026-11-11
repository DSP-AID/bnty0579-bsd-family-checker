# CLAIMS ledger — what is proved, computed, cited, still open

Scope: every rank-0 elliptic curve over Q of conductor N <= 100 (284 curves,
58 conductors; completeness witnessed against all 306 curves N<=100 in 67
conductors via `count.txt` class counts), as an isogeny-invariant BSD identity
per isogeny class.

## PROVED (exact arithmetic in this repo, no unproved input)

1. **Conductor re-derivation.** For every curve and every bad prime p, `a_p =
   p + 1 - #E(F_p)` is recomputed by direct projective point counting; the exact
   sequence `0 -> E_1(F_p) -> E(F_p) -> E_ns(F_p) -> 0` (the formal group over the
   finite field contributes only `O`, so `#E_1(F_p) = 1`) then *forces* `a_p` in
   `{+1, -1, 0}` = split multiplicative / nonsplit multiplicative / additive —
   the type decision is a theorem, valid at p = 2, 3 as well (Hasse alone would
   not pin it there; the group structure does). The conductor exponent is `f_p =
   1` multiplicative, `f_p = 2` additive at p >= 5 (Tate's algorithm output,
   cross-checked against the count), and `f_p = v_p(N)` cited at wild additive
   primes; `prod_p p^f_p == N` exactly, and every root of Delta mod p divides N
   (minimality witness for the input model).
2. **Local Euler factors.** `a_p = p + 1 - #E(F_p)` recomputed by projective point
   counting for every prime p <= 97 (the stored 25-term aplist) matches the db for
   all 284 curves; the a_p are then determined for all p by the Euler product
   `sum a_n n^-s = prod_p (1 - a_p p^-s + p^(1-2s))^-1`.
3. **Rational torsion structure** for each curve: exhibited rational points of
   exact order (group law over Q in Fraction arithmetic), and the checker computes
   the exact order of the subgroup they generate (closure under addition, so
   non-cyclic groups like Z/2 x Z/4 are counted correctly). The upper bound
   `T | gcd_p #E(F_p)` is taken over odd primes of good reduction only: for p odd
   the reduction map is injective on E(Q)_tors (AEC VII.3.1, formal group torsion
   -free over Q_p); p=2 must be excluded because E_1(Q_2) can carry formal
   2-torsion (e.g. 15a1 has a point of order 2 reducing to O mod 2, and indeed
    #E(F_2)=4 does not divide T=8). Two certification tiers: PROVED when the
    exhibited subgroup attains the bound (this covers every curve whose db entry
    displays generator points); CITED-with-condition elsewhere (db entries that
    store no generators, e.g. 11a2/19a2 - then the bound must be an exact
    multiple of the cited order). Every report row marks its tier.
4. **Algebraic side of the identity** (`Omega_E`, Kodaira symbols, `c_p`, `prod c_p`)
   is exact/rational arithmetic given the algorithm outputs (period via mpmath
   quadrature to 30 digits — see COMPUTED). Multiplicative Tamagawa numbers are
   re-derived at *every* prime including wild p = 2, 3: `c_p = n` (split) or
   `gcd(2, n)` (nonsplit), from `#Phi(F_p)` for `I_n` fibers (Silverman, *Advanced
   Arithmetic*, VII.6.5a-b — the component group is `Z/nZ` / `Z/nZ` with inversion
   Galois action; the formula is valid in the wild too). Additive `c_p` is CITED
   (below), and the checker verifies the cited product `prod c_p` (db CP column)
   is an *exact multiple* of the re-derived multiplicative part — equality when
   there are no additive primes at all.

## COMPUTED (rigorous numerics with certified error bounds + internal self-checks)

1. **L(E,1)**: critical L-series `sum a_n n^-s e^{-2n/sqrt(N)}` with an explicit
   tail bound (`2.4 sqrt(N) A N^(1/4) e^-R`), evaluated at 30 digits; the completed
   L-value satisfies the Fricke functional equation (w_N = +1, forced by rank 0) to
   1e-12 relative under cutoff halving. This is a *computation of the analytic
   side*, not a proof that the Dirichlet series converges at 1 — convergence/
   analyticity of L(E,s) is CITED (modularity).
2. **Real period** `Omega_E`: complete-real-period convention (integral over all of
   E(R); = 2*K(m)/sqrt(rho) for Delta<0, = 4*K(k)/(sqrt(e1-e3))... see README
   conventions) computed by mpmath elliptic integral at 30 digits; validated to
   ~1e-15 against the OM column on all 284 curves.
3. **Torsion search**: rational points enumerated with |p| <= 300, q <= 24 in
   x = p/q; combined with the exact gcd upper bound this certifies T for all 284
   bound is attained (where it is not, the db displays no generators and the
   tier is CITED-with-condition, never OPEN; search radius itself is finite).
4. **BSD ratio** `L(E,1)*T^2/(Omega_E * prod_p c_p)` = predicted Sha order; rounded
   to nearest integer and compared with the SHA column: 284/284 match.

## CITED (external results used as reference / input)

- **Modularity & analytic continuation of L(E,s)** (Wiles et al.) — licenses the
  series evaluation at s=1.
- **BSD conjecture** is not used to *derive* anything here: we compare an
  independently computed analytic quantity with an independently stored algebraic
  quantity. Where Sha = 1 and Kolyvagin's theorem applies (root number +1, Heegner
  or Euler-system input), `Sha(E/Q)[p^infty] = 0` for all p is a *theorem*; the
  full finiteness of Sha is CITED from those cases only where recorded in the db.
- **Cremona eclib `ecdata`** (`https://github.com/JohnCremona/ecdata`, master) —
  reference T, CP, OM, L1, SHA, a_p columns; live-fetched by `scripts/fetch_cremona.sh`.
- **LMFDB conventions page** (`/api/1/ec_mwbsd`, curve pages) for the real-period /
  L1/r! normalization; spot-checked via `--live` (best effort — captcha wall).
- **Tate's algorithm, Cohen's course algorithms** (Cohen 7.5.1/7.5.2) as implemented
  in `src/tate.py` — used and relied on only for p >= 5 (tame), where its output is
  cross-checked against the point-count type decision; its wild p = 2, 3 port was
  found (see next bullet) to misclassify some fibers and is *not* relied on.
- **Additive `c_p`** (all additive primes) and **Kodaira symbol / `f_p` at wild
  additive p = 2, 3** (reported as `"?"` in the rows): taken from the Cremona CP
  column via the db `tamagawa_product`, pinned by the divisibility check in
  PROVED-4. Cross-checked against an independent oracle: PARI/GP
  `ellglobalred()` run over all 284 curves — the oracle agrees with the db on
  every `(p, f_p, c_p)` and exposed our `tate.py` wild-port misclassifications
  (e.g. 48a6 at p = 2: truth `I*_3`, c = 4; our port said otherwise). Discrepancies
  are reported, never patched.

## OPEN (explicitly not claimed)

- BSD rank-part formula for any curve here is a *verified identity*, not a proof of
  finiteness of Sha, unless covered by Kolyvagin/Kato theorems cited above.
- Curves of conductor > 100, and rank > 0 curves (the L1 column stores L/r! there;
  generalizing to r>0 needs Heegner-point input for T) — out of scope.
- A corrected *pure-Python* wild (`p = 2, 3`) Tate/Conradi port for additive
  fibers: found broken via the PARI oracle, reported, not patched (out of wrap-up
  scope); everything wild-additive is cited instead, as above.
- The Sha *square class* refinement (ratio is a square times 1/4-power): we only
  assert the rounded integer equals the stored Sha.
