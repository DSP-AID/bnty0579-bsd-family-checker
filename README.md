# bnty_0579 — BSD counterexample hunt: the rank-0 family for conductor N <= 100

Generalization of the `32a2` BSD-identity checker (bounty `bnty_0573`) to **every
rank-0 elliptic curve over Q of conductor N <= 100**: 284 rank-0 curves across 58
conductors (completeness witnessed against all 306 curves of conductor <= 100 in
67 conductors, class counts verified against `count.txt`), each cross-checked
end-to-end against the live reference database.

## What this repo contains

| path | what |
|---|---|
| `src/bsd_family.py` | the generalized checker: conductor re-derivation (a_p point-count type decision + tame Tate/Ogg), Fricke L-series with rigorous tail bound + functional-equation self-check, real period (both Delta signs), rational-torsion certification, Tamagawa product (multiplicative c_p re-derived at every p incl. wild; additive c_p cited, divisibility-checked), full BSD identity vs database |
| `src/tate.py` | Tate's algorithm (Cohen 7.5.1/7.5.2) — Kodaira type, conductor exponent, c_p; used for p >= 5 only (cross-checked against point counts); its wild p = 2,3 port misclassifies some additive fibers (found via a PARI `ellglobalred` oracle) — reported, not patched |
| `scripts/fetch_cremona.sh` | live-fetched Cremona eclib tables (`ecdata@master`, the upstream source of LMFDB EC pages) via HTTP Range heads |
| `scripts/parse_cremona.py` | builds `data/families.json` from the raw tables (allbsd + aplist + count) |
| `scripts/make_report.py` | regenerates `REPORT.md` from a results file (incl. completeness audit) |
| `scripts/run_parallel.sh` | shards the family over N workers on an offload box and merges results |
| `tests/mutation_tests.py` | 12 mutations of stored cross-check values (incl. regressions for the Z/2xZ/4 torsion and additive-type bugs found in the family run) + 4 positive controls — proves every check bites |
| `verify.sh` | entry point: whole family, one label, or `--live` LMFDB spot-check |
| `REPORT.md` | per-curve family table + completeness audit + mismatch list |
| `CLAIMS.md` | claim ledger: PROVED / COMPUTED / CITED / OPEN |
| `results_full.json` | raw per-curve evidence (every computed invariant, 16 digits) |
| `data/` | raw eclib heads (provenance) + parsed family input with db cross-check fields |

## Run it

```bash
./verify.sh                 # all 284 rank-0 curves; pure stdlib + mpmath, ~100 s/curve
./verify.sh 32a2            # one curve: per-invariant table incl. the explicit BSD identity
./verify.sh --refresh       # re-fetch + rebuild the family input from ecdata@master first
./verify.sh 11a1 --live     # optional LMFDB API spot-check (best effort: captcha walls)
sh scripts/run_parallel.sh 16   # sharded run (~40 min wall on a 36-thread Xeon)
```

## Result (this commit)

- **284/284 rank-0 curves PASS** at 1e-12 relative tolerance on: conductor
  re-derivation, every a_p (p <= 97), L(E,1)/Omega vs functional equation,
  real period, torsion, prod c_p, and the BSD ratio `L(E,1)*T^2/(Omega*prod c_p)`
  equalling the stored Sha.
- **0 verified mismatches** against the reference database — and no fabrication:
  two real checker bugs were caught by running over the family (lcm-of-orders
  undercounting Z/2xZ/4 torsion; additive II/III/IV types misread as
  multiplicative) and are pinned by dedicated regression mutations.
- Completeness: raw eclib heads contain every curve for all N <= 100; parsed
  isogeny-class counts match `count.txt` NNF for 67/67 conductors; the rank-0
  subset (Mordell-Weil rank column of `allbsd`) is exactly these 284 curves.

## Data & licenses

Code: MIT (`LICENSE`). Computed tables derived from John Cremona's eclib `ecdata`
(`https://github.com/JohnCremona/ecdata`, CC-BY-SA-4.0) — the derived JSON and this
report are therefore CC-BY-SA-4.0 as well; see `DATA_LICENSE`.
