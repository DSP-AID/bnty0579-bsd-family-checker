"""Tate's algorithm over Q: Kodaira symbol, conductor exponent f_p, Tamagawa
number c_p, and the valuation of the local minimal discriminant.

Implements Cohen, "A Course in Computational Algebraic Number Theory",
Algorithms 7.5.1 (p >= 5) and 7.5.2 (p = 2, 3), including the I*_n loop
(steps 16-21). For p >= 5 the algorithm is invariant-based (c4, c6, Delta and
residue polynomials in those invariants); for p = 2, 3 it follows Cohen's
coefficient transformations verbatim. Splitting tests over F_p use plain root
counts, which classify quadratics/cubics correctly for every residue degree
and also at p = 2, 3 (no discriminant sign arguments needed).

Conductor exponent f is read off each Cohen step (equivalently Ogg's formula
f = v - m + 1 with m the number of geometric components).

Pure stdlib, exact integer arithmetic. Input must be an integral model; a
globally minimal input model (as LMFDB provides) never triggers the
non-minimal branches, but they are implemented for completeness.
"""
from __future__ import annotations

from math import gcd


def b_invariants(a):
    a1, a2, a3, a4, a6 = a
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    return b2, b4, b6, b8


def c_invariants(a):
    b2, b4, b6, b8 = b_invariants(a)
    c4 = b2 * b2 - 24 * b4
    c6 = -b2 ** 3 + 36 * b2 * b4 - 216 * b6
    delta = -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6
    return c4, c6, delta


def val(n, p):
    if n == 0:
        return 10 ** 9
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def legendre(v, p):
    r = pow(v % p, (p - 1) // 2, p)
    return 0 if r == p - 1 else r


def transform(a, u, r, s, t):
    """Cohen formulas (7.2); asserts exact divisibility (integrality kept)."""
    a1, a2, a3, a4, a6 = a
    n1 = a1 + 2 * s
    n2 = a2 - s * a1 + 3 * r - s * s
    n3 = a3 + r * a1 + 2 * t
    n4 = a4 - s * a3 + 2 * r * a2 - (t + r * s) * a1 + 3 * r * r - 2 * s * t
    n6 = a6 + r * a4 + r * r * a2 + r ** 3 - t * a3 - t * t - r * t * a1
    assert n1 % u == 0 and n2 % u ** 2 == 0 and n3 % u ** 3 == 0
    assert n4 % u ** 4 == 0 and n6 % u ** 6 == 0, "transformation not integral"
    return [n1 // u, n2 // u ** 2, n3 // u ** 3, n4 // u ** 4, n6 // u ** 6]


def roots_mod_p(coeffs, p):
    """All roots in Z/pZ of the polynomial coeffs (highest power first)."""
    return [x for x in range(p) if _poly(coeffs, x) % p == 0]


def _poly(coeffs, x):
    acc = 0
    for c in coeffs:
        acc = acc * x + c
    return acc


def _poly_gcd(f, g, p):
    """Monic gcd of polynomials f, g over F_p (coeff lists, highest first)."""
    def trim(f):
        f = [c % p for c in f]
        while len(f) > 1 and f[0] == 0:
            f = f[1:]
        return f

    def mod(f, g):
        f = trim(f)[:]
        g = trim(g)[:]
        if len(g) == 1:
            return [0]
        inv = pow(g[0], p - 2, p)
        while not (len(f) == 1 and f[0] == 0) and len(f) >= len(g):
            k = f[0] * inv % p
            for i in range(len(g)):
                f[i] = (f[i] - k * g[i]) % p
            f = trim(f)
        return f

    f, g = trim(f), trim(g)
    while trim(g) != [0]:
        f, g = g, mod(f, g)
    lead_inv = pow(f[0], p - 2, p)
    return [c * lead_inv % p for c in f]


def _classify_cubic(cub, p):
    """('distinct', None) | ('double', alpha) | ('triple', alpha).

    'distinct' means three distinct roots over Fbar_p (the I*_0 case);
    otherwise alpha is the multiple root (a lift to Z of its F_p value).
    """
    cub = [c % p for c in cub]
    der = [(3 * cub[0]) % p, (2 * cub[1]) % p, cub[2] % p]
    if all(d % p == 0 for d in der):  # p=3, P = X^3 + C = (X+a)^3 (F_p perfect)
        a0 = next(x for x in range(p) if _poly(cub, x) % p == 0)
        return "triple", a0
    g = _poly_gcd(cub, der, p)
    if len(g) == 1:
        return "distinct", None
    roots_g = roots_mod_p(g, p)
    assert len(roots_g) == 1, f"unexpected multiple-root structure gcd={g}"
    a0 = roots_g[0]
    return ("double", a0) if len(g) == 2 else ("triple", a0)


def _quad_status(quad, p):
    """('split', None) | ('double', alpha) | ('inert', None) | ('linear', None).

    Over a finite field a quadratic has 2/1/0 roots in F_p iff it splits into
    distinct factors / is a square of a linear factor / is irreducible, so a
    plain root count classifies it for every p (no discriminant needed).
    """
    quad = [c % p for c in quad]
    if len(quad) == 3 and quad[0] == 0:
        quad = quad[1:]
    if len(quad) == 2:
        return "linear", None
    rs = roots_mod_p(quad, p)
    if len(rs) == 2:
        return "split", None
    if len(rs) == 1:
        return "double", rs[0]
    return "inert", None


def tate(a, p):
    """Return {'symbol','f','c','v'} at prime p for integral model a=[a1..a6]."""
    assert all(isinstance(x, int) for x in a), "integral coefficients required"
    if p >= 5:
        return _tate_large(list(a), p)
    return _tate_small(list(a), p)


def _out(sym, f, c, v):
    return {"symbol": sym, "f": f, "c": c, "v": v}


# ------------------------------------------------------------- p >= 5 (7.5.1)

def _tate_large(a, p):
    c4, c6, delta = c_invariants(a)
    vd = val(delta, p)
    if vd == 0:
        return _out("I0", 0, 1, 0)
    vj = (3 * val(c4, p) - vd) if c4 != 0 else 10 ** 9
    if vj < 0:  # step 3: potentially multiplicative, k = 3*v(c4) in {0, 6}
        vc4 = val(c4, p)
        assert 3 * vc4 == vd + vj, "step 3 k-mismatch"
        if vc4 == 0:  # I_n multiplicative
            n = vd
            assert val(c6, p) == 0
            split = legendre(-c6, p) == 1
            return _out(f"I{n}", 1, n if split else gcd(2, n), n)
        assert vc4 == 2, f"step 3: v(c4)={vc4} not minimal?"
        nv = vd - 6
        assert nv >= 1
        # I*_n (Cohen step 3 k=6 branch): odd n via (Delta' * c6'), even n via
        # (Delta'); Delta' = Delta / p^v, c6' = c6 / p^3 (units at p >= 5).
        dlt = delta // p ** vd
        if nv % 2:
            arg = dlt * (c6 // p ** val(c6, p))
            c = 3 + legendre(arg, p)
        else:
            c = 3 + legendre(dlt, p)
        return _out(f"I{nv}*", 2, c, vd)
    # step 4: integral j; a minimal model at p >= 5 has v(delta) in
    # {0, 2, 3, 4, 6, 8, 9, 10}; v(delta) >= 12 would mean non-minimal input.
    k = vd
    assert k < 12, "non-minimal input model at p >= 5"
    if k == 2:
        return _out("II", 2, 1, 2)
    if k == 3:
        return _out("III", 2, 2, 3)
    if k == 4:
        assert val(c6, p) == 2
        return _out("IV", 2, 2 + legendre(-6 * (c6 // p ** 2), p), 4)
    if k == 6:
        assert val(c4, p) >= 2 and val(c6, p) >= 3
        cub = [4, 0, -3 * (c4 // p ** 2), -(c6 // p ** 3)]
        return _out("I0*", 2, 1 + len(roots_mod_p(cub, p)), 6)
    if k == 8:
        assert val(c6, p) == 4
        return _out("IV*", 2, 2 + legendre(-6 * (c6 // p ** 4), p), 8)
    if k == 9:
        return _out("III*", 2, 2, 9)
    if k == 10:
        return _out("II*", 2, 1, 10)
    raise AssertionError(f"tate p>=5: unhandled v(delta)={vd} v(j)={vj}")


# -------------------------------------------------------- p in {2, 3} (7.5.2)

def _tate_small(a, p):
    """Tate's algorithm at p in {2,3} (Cohen 7.5.2), coefficients named."""
    for _ in range(64):  # step 15 loops back; a minimal input never loops
        _c4, _c6, delta = c_invariants(a)
        v = val(delta, p)
        if v == 0:  # step 2
            return _out("I0", 0, 1, 0)
        b2, _b4, b6, _b8 = b_invariants(a)  # b2 -> step 3; b6 -> p=3 branch of step 4
        a1, a2, a3, a4, a6 = a
        if b2 % p:  # step 3: multiplicative
            split = len(roots_mod_p([1, a1, -a2], p)) > 0
            return _out(f"I{v}", 1, v if split else gcd(2, v), v)
        # step 4: change equation to force p | b2 (u = 1 keeps v(delta))
        if p == 2:
            r1 = a4 % 2
            s1 = (r1 + a2) % 2
            t1 = (a6 + r1 * (a4 + s1)) % 2
        else:
            r1 = (-b6) % 3
            s1 = a1 % 3
            t1 = (a3 + r1 * a1) % 3
        a = transform(a, 1, r1, s1, t1)
        b2, _b4, b6, b8 = b_invariants(a)
        a1, a2, a3, a4, a6 = a
        assert b2 % p == 0, "step 4 did not force p | b2"
        if val(a6, p) < 2:  # step 5: II (v unchanged by u=1 transforms)
            return _out("II", v, 1, v)
        if val(b8, p) < 3:  # step 6: III
            return _out("III", v - 1, 2, v)
        if val(b6, p) < 3:  # step 7: IV
            q = [1, a3 // p, -(a6 // p**2)]
            return _out("IV", v - 2, 3 if len(roots_mod_p(q, p)) else 1, v)
        if val(a6, p) < 3:  # step 8: y-shift to force p^3 | a6
            k = 2 if p == 2 else a2 % 9
            a = transform(a, 1, 0, 0, k)
            a1, a2, a3, a4, a6 = a
        cub = [1, a1 // p, a3 // p**2, a6 // p**3]  # steps 9/10
        kind, alpha = _classify_cubic(cub, p)
        if kind == "distinct":
            return _out("I0*", v - 4, 1 + len(roots_mod_p(cub, p)), v)
        if alpha % p:  # step 10 x-shift moves the multiple root to 0
            a = transform(a, 1, alpha * p, 0, 0)
            a1, a2, a3, a4, a6 = a
        if kind == "double":
            return _tate_in_star(a, p, v)  # steps 16-21
        # triple root at 0: step 11
        quad = [1, a3 // p**2, -(a6 // p**4)]
        st = _quad_status(quad, p)
        if st[0] in ("split", "inert"):
            return _out("IV*", v - 6, 3 if st[0] == "split" else 1, v)
        dr = st[1]
        if dr % p:  # step 12
            a = transform(a, 1, 0, 0, dr * p * p)
            a1, a2, a3, a4, a6 = a
        if val(a4, p) < 4:  # step 13
            return _out("III*", v - 7, 2, v)
        if val(a6, p) < 6:  # step 14
            return _out("II*", v - 8, 1, v)
        a = transform(a, p, 0, 0, 0)  # step 15: non-minimal, restart at step 1
    raise AssertionError("tate_small iteration limit")


def _tate_in_star(a, p, v):
    """Cohen steps 16-21: the I*_n loop (double root already moved to 0)."""
    n, q = 1, p * p
    for _ in range(64):
        a1, _a2, a3, _a4, a6 = a
        if n % 2:  # step 17
            quad = [1, a3 // q, -(a6 // q**2)]
        else:  # step 19
            quad = [a1 // p, a3 // (p * q), a6 // (p * q * q)]
        st = _quad_status(quad, p)
        if st[0] in ("split", "inert", "linear"):
            return _out(f"I{n}*", v - n - 4, 4 if st[0] == "split" else 2, v)
        dr = st[1]
        if n % 2:  # step 18: y-shift by a*q
            if dr % p:
                a = transform(a, 1, 0, 0, dr * q)
        else:  # step 20: x-shift by a*q
            if dr % p:
                a = transform(a, 1, dr * q, 0, 0)
        n += 1
        if n % 2 == 1:  # step 21
            q *= p
    raise AssertionError("tate_in_star iteration limit")
