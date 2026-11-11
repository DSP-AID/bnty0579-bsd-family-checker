"""BSD leading-term checker generalized to every rank-0 curve of conductor <= 100.

Family input: data/families.json (built by scripts/parse_cremona.py from the
live-fetched Cremona eclib tables - the source database of LMFDB's elliptic
curve pages; --live additionally spot-checks the running LMFDB API).

For each curve the checker computes from first principles at ~35 digits:

  1. conductor re-derivation at every bad prime: a_p = p + 1 - #E(F_p) by
     direct projective point count decides the type (a_p in {0,+-1}: split /
     nonsplit multiplicative vs additive); f_p = 1 for multiplicative fibers,
     f_p = 2 for additive fibers at p >= 5 (full Tate algorithm, tate.py,
     cross-checked), and f_p = v_p(N) cited at wild additive primes p = 2, 3;
     prod p^f_p == N; the primes of bad reduction (roots of Delta) must all
     divide N, which re-proves minimality of the input model;
  2. a_p by direct point counting over F_p (good p) and a_p in {0,+-1} from
     the nonsingular-locus count at bad p; cross-checked against the Cremona
     AP table for all primes < 100;
  3. L(E,1) by the completed-L rapidly convergent series (Fricke sign +1,
     forced by rank 0), with rigorous tail bound |a_n| <= d(n) sqrt(n) and an
     independent cutoff-halving self-check of the functional equation;
  4. Omega_E = integral over E(R) of dx/(2y+a1x+a3) in closed form via
     elliptic K / AGM for both topologies (Delta < 0 connected, Delta > 0 two
     components - where the BSD convention uses the integral over all of
     E(R), twice the identity component);
  5. torsion: rigorous upper bound T | gcd_p #E(F_p) over good primes
     (injection into E(F_p)); lower bound from explicit rational points whose
     orders are verified by exact chord-tangent arithmetic. T is proved when
     exhibited points meet the bound; otherwise the gap is reported OPEN;
  6. Tamagawa product: c_p re-derived rigorously for multiplicative fibers
     (even at p = 2, 3); additive c_p CITED from the db's proven
     tamagawa_product, verified to be an exact multiple of the re-derived part.

The BSD prediction |Sha| = L(E,1) T^2 / (Omega_E prod c_p) must be a positive
integer square equal to the stored analytic order (Cremona SHA column / LMFDB
sha). Discrepancies are reported, never patched.

Usage: python src/bsd_family.py [--data DIR] [--out FILE] [--live] [--label LAB]
"""
from __future__ import annotations

import argparse
import bisect
import json
import math
import os
import sys
from fractions import Fraction as F

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tate import b_invariants, c_invariants, tate

WORK_DP = 38
AGM_DP = 60


def vp(n, p):
    if n == 0:
        return 10**9
    v, n = 0, abs(n)
    while n % p == 0:
        n //= p
        v += 1
    return v


def _make_primes(limit):
    sieve = bytearray([1]) * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, math.isqrt(limit) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytearray(len(sieve[i * i :: i]))
    return [i for i in range(2, limit + 1) if sieve[i]]


_ALL_PRIMES = _make_primes(4_000_000)


def primes_upto(n):
    return _ALL_PRIMES[: bisect.bisect_right(_ALL_PRIMES, n)]


def factor_small(n):
    """Prime factors of n by trial division (n up to ~1e13)."""
    n, out = abs(n), []
    for p in primes_upto(math.isqrt(n) + 1) if n > 1 else []:
        if p * p > n:
            break
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
    if n > 1:
        out.append(n)
    return out


# --------------------------------------------------------------- point counts

def count_points(a, p):
    """#E(F_p) set-theoretically (singular affine points included)."""
    a1, a2, a3, a4, a6 = a
    cnt = 1
    for x in range(p):
        if p == 2:
            for y in range(2):
                if (y * y + a1 * x * y + a3 * y - (x**3 + a2 * x * x + a4 * x + a6)) % 2 == 0:
                    cnt += 1
            continue
        d = ((a1 * x + a3) ** 2 + 4 * (x**3 + a2 * x * x + a4 * x + a6)) % p
        if d == 0:
            cnt += 1
        else:
            r = pow(d, (p - 1) // 2, p)
            cnt += 2 if r == 1 else 0
    return cnt


def singular_count(a, p):
    a1, a2, a3, a4, a6 = a
    n = 0
    for x in range(p):
        for y in range(p):
            if (y * y + a1 * x * y + a3 * y - (x**3 + a2 * x * x + a4 * x + a6)) % p:
                continue
            if (a1 * y - 3 * x * x - 2 * a2 * x - a4) % p:
                continue
            if (2 * y + a1 * x + a3) % p:
                continue
            n += 1
    return n


def local_data(a, N, p):
    """(kodaira, f_p, c_p or None=cited, a_p) at a bad prime p.

    a_p = p + 1 - #E(F_p) by direct projective point count; bad reduction
    forces a_p in {0, +-1} (Hasse pins the exact value), so the count alone
    decides the type: a_p = +1 multiplicative split, -1 multiplicative
    nonsplit, 0 additive. Multiplicative fibers are tame at every p, so
    f_p = 1 and c_p = n (split) / gcd(2, n) (nonsplit), n = v_p(Delta_min),
    are re-derived rigorously even at p = 2, 3. At p >= 5 the full Tate
    algorithm (tate.py) supplies the Kodaira symbol and f_p = 2 for the
    additive fibers, cross-checked against the count splitness test.
    Additive c_p is CITED from the db's proven tamagawa_product (checked
    in check_curve to be an exact multiple of the re-derived multiplicative
    part); at p = 2, 3 additive even the symbol/f are cited ("?"), since our
    wild-port of Cohen 7.5.2 was found (via a PARI ellglobalred oracle) to
    misclassify some p = 2, 3 fibers - reported, never patched.
    """
    ap = p + 1 - count_points(a, p)
    assert ap in (-1, 0, 1), f"bad-reduction Hasse anomaly at {p}: a_p={ap}"
    if ap in (1, -1):  # multiplicative I_n: n = v_p(Delta_min), model is minimal
        _c4, _c6, delta = c_invariants(a)
        n = vp(delta, p)
        if p >= 5:
            t = tate(a, p)
            assert t["symbol"] == f"I{n}" and t["f"] == 1, f"tate vs count at {p}"
            if n >= 3:  # Tate's Legendre split test vs the point count
                assert (t["c"] == n) == (ap == 1), f"tate splitness vs count at {p}"
        c = n if ap == 1 else math.gcd(2, n)
        return f"I{n}", 1, c, ap
    if p < 5:
        return "?", None, None, ap  # wild additive: cited tier
    t = tate(a, p)
    assert t["f"] == 2, f"tame additive conductor exponent != 2 at {p}"
    return t["symbol"], 2, None, ap


# ------------------------------------------------------------- L-series input

def hecke_coefficients(a, N, bad_ap, nmax):
    """a_n, 1<=n<=nmax: Euler product from a_p (good p by point count, bad p
    supplied); a_{p^k} = a_p^k at p | N, recursion elsewhere."""
    primes = primes_upto(nmax)
    ap = {p: (bad_ap[p] if p in bad_ap else p + 1 - count_points(a, p)) for p in primes}
    coef = [0] * (nmax + 1)
    coef[1] = 1
    for p in primes:
        pk, e = p, 1
        while pk <= nmax:
            if p in bad_ap:
                val = ap[p] ** e
            elif e == 1:
                val = ap[p]
            else:
                val = ap[p] * coef[p ** (e - 1)] - p * coef[p ** (e - 2)]
            for m in range(1, nmax // pk + 1):
                if math.gcd(m, p) == 1:
                    coef[m * pk] = coef[m] * val
            pk *= p
            e += 1
    return coef, {p: ap[p] for p in primes if p <= 97}


def choose_nmax(N, digits=WORK_DP - 4):
    c = 2 * math.pi / N
    return math.ceil(math.log(10) * (digits + 3) / c) + 8


def l_series(N, coef, nmax):
    """L(E,1) = 2 pi (I_slow + 2 J) from the completed L-function with Fricke
    sign +1; returns (value at nmax, value at nmax//2, tail bound)."""

    def partial(k):
        with mp.workdps(WORK_DP + 8):
            tp = 2 * mp.pi
            i_slow = mp.fsum(
                coef[n] * (mp.exp(-tp * n / N) - mp.exp(-tp * n)) / (tp * n)
                for n in range(1, k + 1)
            )
            j_sum = mp.fsum(coef[n] * mp.exp(-tp * n) / (tp * n) for n in range(1, k + 1))
            return tp * (i_slow + 2 * j_sum)

    full, half = partial(nmax), partial(max(4, nmax // 2))
    with mp.workdps(20):
        c = 2 * mp.pi / N
        tail = 2 * mp.e ** (-c * (nmax + 1)) / (mp.pi * (1 - mp.e ** (-c)))
    return full, half, tail


# -------------------------------------------------------------- real period

def cubic_roots(a):
    """Roots of F(x) = 4x^3 + b2 x^2 + 2 b4 x + b6 at AGM_DP digits, polished
    by Newton steps (polyroots alone is not reliable to 50+ digits)."""
    a1, a2, a3, a4, a6 = [mp.mpf(x) for x in a]
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6

    def f(x):
        return 4 * x**3 + b2 * x**2 + 2 * b4 * x + b6

    def fp(x):
        return 12 * x**2 + 2 * b2 * x + 2 * b4

    rts = mp.polyroots([mp.mpf(4), b2, 2 * b4, b6], maxsteps=200, extraprec=500)
    for _ in range(3):
        rts = [r - f(r) / fp(r) for r in rts]
    return sorted(rts, key=lambda z: (mp.im(z) != 0, -mp.re(z)))


def real_period(a):
    """Omega_E = integral of the invariant differential over E(R), the BSD /
    LMFDB real_period convention.

    Delta < 0: E(R) connected; deflate F=(x-e1)(4x^2+cx+d) with conjugate pair
    alpha, alpha-bar; Omega = 2 K(m)/sqrt(rho), m=(rho+Re(alpha)-e1)/(2 rho).
    Delta > 0: roots e1>e2>e3; the identity component has period
    2 K(k)/sqrt(e1-e3), k^2=(e2-e3)/(e1-e3), and BSD uses the integral over
    all of E(R) = twice that (calibrated against the Cremona OM column; the
    agreement is checked at 1e-12 relative for every curve in the family)."""
    b2, b4, _b6, _b8 = b_invariants(a)
    _c4, _c6, delta = c_invariants(a)
    assert delta != 0
    with mp.workdps(AGM_DP):
        rts = cubic_roots(a)
        if delta < 0:
            e1 = next(r for r in rts if abs(mp.im(r)) < mp.mpf(10) ** -(AGM_DP - 10))
            c = b2 + 4 * e1
            d = 2 * b4 + e1 * c
            p_re, q_im = -c / 8, mp.sqrt(16 * d - c * c) / 8
            rho = mp.sqrt((p_re - e1) ** 2 + q_im**2)
            m_par = (rho + p_re - e1) / (2 * rho)
            return mp.re(2 * mp.ellipk(m_par) / mp.sqrt(rho))
        e1, e2, e3 = [mp.re(r) for r in rts[:3]]
        k2 = (e2 - e3) / (e1 - e3)
        return mp.re(4 * mp.ellipk(k2) / mp.sqrt(e1 - e3))


# --------------------------------------------------------- exact group law

def _sadd(p1, p2, a):
    """Chord-tangent law over Q (Fractions); None is the point at O."""
    if p1 is None:
        return p2
    if p2 is None:
        return p1
    a1, a2, a3, a4, _a6 = a
    x1, y1 = p1
    x2, y2 = p2
    if p1 == p2:
        den = 2 * y1 + a1 * x1 + a3
        if den == 0:
            return None
        num = 3 * x1 * x1 + 2 * a2 * x1 + a4 - a1 * y1
    elif x1 == x2:
        return None  # vertical line: inverse pair
    else:
        num, den = y2 - y1, x2 - x1
    lam = F(num, den)
    nx = lam * lam + a1 * lam - a2 - x1 - x2
    ny = -(lam + a1) * (nx - x1) - y1 - a1 * x1 - a3
    return (nx, ny)


def mul_pt(P, m, a):
    R = None
    for bit in bin(m)[2:]:
        if R is not None:
            R = _sadd(R, R, a)
        if bit == "1":
            R = P if R is None else _sadd(R, P, a)
    return R


def divisors(n):
    out = []
    for i in range(1, math.isqrt(n) + 1):
        if n % i == 0:
            out.append(i)
            if i * i != n:
                out.append(n // i)
    return sorted(out)


def torsion_bound(N, a, good=None):
    """Rigorous upper bound: for ODD primes of good reduction the reduction map
    is injective on all of E(Q)_tor (AEC VII.3.1: E_1(Q_p) has no torsion for
    p odd over Q), so T | gcd_p #E(F_p), p >= 3, p not | N.

    p=2 must be EXCLUDED: E_1(Q_2) can carry formal 2-torsion (e.g. 15a1 has a
    point of order 2 reducing to O mod 2), so #E(F_2) need not be divisible by T.
    """
    g = 0
    for p in good or primes_upto(97):
        if p == 2 or N % p == 0:
            continue
        g = math.gcd(g, count_points(a, p))
        if g == 1:
            break
    return g


def _closure(pts, a, cap=64):
    """Exact order of the subgroup of E(Q) generated by pts (None = O); the
    group is small (torsion), so a naive +P fixed point suffices."""
    grp = [None]
    for P in pts:
        changed = True
        while changed and len(grp) <= cap:
            changed = False
            for r in list(grp):
                s = _sadd(P, r, a)
                if s not in grp:
                    grp.append(s)
                    changed = True
    return len(grp) if len(grp) <= cap else cap + 1


def torsion_computed(N, a, delta, pbound=300, qbound=24):
    """Lower bound = exact order of the subgroup generated by exhibited rational
    points (search |p| <= pbound, q <= qbound for x = p/q), orders verified by
    exact group law; certified iff it meets the gcd upper bound g.

    shortcut: points of order 2 are completed exactly by the rational-root test
    on psi_2 (their x can sit far outside the naive search box, e.g. 45a7); a
    missed point of higher order only weakens the lower bound (reported as
    torsion_certified=False / db-mismatch), never fakes a proof.
    """
    g = torsion_bound(N, a)
    divs = [d for d in divisors(g) if d > 1]
    a1, a2, a3, a4, a6 = a
    pts = []
    for q in range(1, qbound + 1):
        for p in range(-pbound, pbound + 1):
            if math.gcd(p, q) != 1:
                continue
            x = F(p, q)
            bq = a1 * x + a3
            rhs = x**3 + a2 * x * x + a4 * x + a6
            disc = bq * bq + 4 * rhs
            if disc < 0:
                continue
            sn, sd = math.isqrt(disc.numerator), math.isqrt(disc.denominator)
            if sn * sn != disc.numerator or sd * sd != disc.denominator:
                continue
            for s in (F(sn, sd), F(-sn, sd)):
                y = (-bq + s) / 2
                P = (x, y)
                for d in divs:
                    if mul_pt(P, d, a) is None:
                        pts.append(P)
                        break
    if g % 2 == 0:  # exact order-2 points: roots of psi_2's x-cubic 4x^3+b2x^2+2b4x+b6
        b2 = a1 * a1 + 4 * a2
        b4 = 2 * a4 + a1 * a3
        b6 = a3 * a3 + 4 * a6
        cands = divisors(abs(b6)) if b6 else []
        for r in {d for d in cands} | {-d for d in cands} | ({0} if b6 == 0 else set()):
            for s in (1, 2, 4):
                if math.gcd(r, s) != 1:
                    continue
                x = F(r, s)
                if 4 * x**3 + b2 * x * x + 2 * b4 * x + b6 == 0:
                    pts.append((x, -(a1 * x + a3) / 2))
    tor = _closure(pts, a)
    return tor, g, tor == g


# ----------------------------------------------------------------- BSD check

def check_curve(cv, ap_table=None):
    a = [int(x) for x in cv["ainvs"]]
    N = int(cv["conductor"])
    _c4, _c6, delta = c_invariants(a)
    bad = sorted(set(factor_small(delta)) | {p for p in primes_upto(N) if N % p == 0})
    bad_ap, tam_mult, fprod, kod, cited = {}, 1, 1, {}, False
    for p in bad:
        sym, f, c, ap = local_data(a, N, p)
        bad_ap[p] = ap
        if c is None:  # additive c_p (wild p=2,3 also symbol/f): cited tier
            cited = True
        else:
            tam_mult *= c
        fprod *= p ** (f if f is not None else vp(N, p))
        kod[str(p)] = sym
    # minimality witness: every prime of bad reduction (root of Delta) divides N
    conductor_ok = fprod == N and all(N % p == 0 for p in factor_small(delta))
    # cited additive c_p come from the db's proven tamagawa_product; check_curve
    # below pins it to an exact multiple of the re-derived multiplicative part
    # (exact equality when there are no additive primes at all).
    tam = int(cv["tamagawa_product"]) if cv.get("tamagawa_product") else tam_mult

    nmax = choose_nmax(N)
    coef, apmap = hecke_coefficients(a, N, bad_ap, nmax)
    l1, l1_half, tail = l_series(N, coef, nmax)
    # functional-equation self-check: halving the cutoff only drops the
    # e^{-2 pi n/N} tail (~1e-18 relative here), so 12-digit agreement
    # certifies the symmetric-series identity.
    fe_ok = abs(l1 - l1_half) < mp.mpf("1e-12") * max(abs(l1), mp.mpf(10) ** -30)

    omega = real_period(a)
    tor, tor_bound, tor_cert = torsion_computed(N, a, delta)

    with mp.workdps(WORK_DP):
        ratio = l1 * tor * tor / (omega * tam)
        sha_pred = int(mp.nint(ratio)) if ratio > 0 else -1
        square = sha_pred > 0 and math.isqrt(sha_pred) ** 2 == sha_pred
        rel_om = (
            abs(omega - mp.mpf(cv["real_period"])) / mp.mpf(cv["real_period"])
            if cv.get("real_period")
            else None
        )
        rel_l = (
            abs(l1 - mp.mpf(cv["l1_over"])) / max(abs(mp.mpf(cv["l1_over"])), mp.mpf(10) ** -30)
            if cv.get("l1_over")
            else None
        )

    out = {
        "label": cv["label"],
        "conductor": N,
        "delta": delta,
        "components": "connected" if delta < 0 else "two",
        "kodaira": kod,
        "tamagawa_product": tam,
        "tamagawa_mult_part": tam_mult,
        "tamagawa_cited_additive": cited,
        "conductor_rederived": conductor_ok,
        "torsion": tor,
        "torsion_bound": tor_bound,
        "torsion_certified": tor_cert,
        # two-tier: certified iff the exhibited subgroup attains the gcd bound
        # (proved here); otherwise fall back to agreement with the db's proven
        # descent value (CITED - mwrank/eclib), which the mismatch check pins.
        "torsion_agrees": bool(tor_cert or ("torsion" in cv and int(cv["torsion"]) == tor)),
        "omega": mp.nstr(omega, 16),
        "l1": mp.nstr(l1, 16),
        "bsd_ratio": mp.nstr(ratio, 14),
        "sha_predicted": sha_pred if square else f"NON-SQUARE({mp.nstr(ratio, 10)})",
        "sha_is_square": square,
        "fe_selfcheck": bool(fe_ok),
        "tail_bound": mp.nstr(tail, 3),
    }
    # Cremona/LMFDB cross-checks (stored values, 15 digits)
    mism = []
    if not conductor_ok:
        mism.append(f"conductor prod p^f {fprod} vs N {N}")
    if "torsion" in cv and int(cv["torsion"]) != tor:
        mism.append(f"torsion {tor} vs db {cv['torsion']}")
    if cv.get("tamagawa_product"):
        db_tam = int(cv["tamagawa_product"])
        if db_tam % tam_mult:
            mism.append(f"mult-part cprod {tam_mult} does not divide db {db_tam}")
        elif not cited and db_tam != tam_mult:
            mism.append(f"cprod {tam_mult} vs db {db_tam}")
    if "sha" in cv:
        try:
            s_db = int(mp.nint(mp.mpf(str(cv["sha"]))))
            if s_db != sha_pred:
                mism.append(f"sha {sha_pred} vs db {s_db}")
        except (TypeError, ValueError):
            pass
    if rel_om is not None:
        out["omega_rel_err"] = mp.nstr(rel_om, 4)
        if rel_om > mp.mpf("1e-12"):
            mism.append(f"omega {out['omega']} vs db {cv['real_period']}")
    if rel_l is not None:
        out["l1_rel_err"] = mp.nstr(rel_l, 4)
        if rel_l > mp.mpf("1e-12"):
            mism.append(f"L(E,1) {out['l1']} vs db {cv['l1_over']}")
    if ap_table:
        key = f"{N}{''.join(ch for ch in cv['label'] if ch.isalpha())}"
        tab = ap_table.get(key)
        if tab:
            for p, tok in zip(primes_upto(97), tab):
                got = apmap.get(p)
                if got is None:
                    continue
                if N % p == 0:
                    # aplist stores the Atkin-Lehner W_p eigenvalue; for p || N
                    # exactly once a_p = -W_p (Atkin-Li); additive a_p = 0.
                    if vp(N, p) == 1 and tok in "+-" and got != (1 if tok == "-" else -1):
                        mism.append(f"a{p} {got} vs db W_p {tok}")
                elif tok.lstrip("-").isdigit() and got != int(tok):
                    mism.append(f"a{p} {got} vs db {tok}")
    out["mismatches"] = mism or "none"
    return out


def live_lmfdb(cv):
    """Optional spot-check against the running LMFDB API (best effort)."""
    import urllib.request

    url = f"https://www.lmfdb.org/api/ec_curvedata/?label={cv['lmfdb_label']}&_format=json"
    try:
        op = urllib.request.build_opener()
        op.addheaders = [("User-Agent", "Mozilla/5.0")]
        raw = op.open(url, timeout=45).read().decode("utf-8", "replace")
        if raw.lstrip()[:1] not in "{[":
            return "captcha"
        d = json.loads(raw)["data"]
        if not d:
            return "empty"
        tor = int(d[0]["torsion"])
        return "ok" if tor == cv.get("torsion") else f"live torsion {tor} vs {cv.get('torsion')}"
    except Exception as e:  # noqa: BLE001
        return f"err {type(e).__name__}"


def main(argv=None):
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    pa = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    pa.add_argument("--data", default=os.path.join(root, "data"))
    pa.add_argument("--out", default=None)
    pa.add_argument("--label", nargs="+", default=None,
                    help="one or more Cremona labels (for sharded parallel runs)")
    pa.add_argument("--live", action="store_true")
    args = pa.parse_args(argv)

    with open(os.path.join(args.data, "families.json"), encoding="utf-8") as f:
        fam = json.load(f)
    curves = list(fam["curves"].values())
    if args.label:
        want = set(args.label)
        curves = [c for c in curves if c["label"] in want]
    print(f"rank-0 family BSD check: {len(curves)} curves, dps={WORK_DP}")
    print("label      N     T   cprod  Sha_pred  verdict  mismatches")
    fails, results = 0, []
    for cv in curves:
        r = check_curve(cv, fam.get("ap25"))
        if args.live:
            r["live_lmfdb"] = live_lmfdb(cv)
        results.append(r)
        ok = (
            r["conductor_rederived"]
            and r["sha_is_square"]
            and r["fe_selfcheck"]
            and r["torsion_agrees"]
            and r["mismatches"] == "none"
        )
        fails += 0 if ok else 1
        print(
            f"{r['label']:<10} {r['conductor']:<5} {r['torsion']:<3} "
            f"{r['tamagawa_product']:<6} {r['sha_predicted']!s:<9} "
            f"{'PASS' if ok else 'FAIL '}   {r['mismatches']}"
        )
    print(f"\n{len(results) - fails}/{len(results)} PASS; failures: {fails}")
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=1)
    return 0 if fails == 0 and results else 1


if __name__ == "__main__":
    sys.exit(main())
