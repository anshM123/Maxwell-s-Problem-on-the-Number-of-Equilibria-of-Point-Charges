"""AUDIT_2 item 4d (continued): independent rigorous re-certification of the BAND and CORNER charts in mpmath.iv
(correctly rounded; no libm assumption), with the auditor's own implementation of the formulas (derived and checked in
AUDIT_2 sections 2-4), on the same initial grids as code/mixed_cert_all.py, with local bisection on failure.

band:   A12 = e in [0, e0], lambda in [0,1], rho in [-1,1]; quantity I/e^4 in scaled variables (psi = e psih,
        S = e Sh, beta' = e bph, c_i = e c_ih, Delta = e Dh) and (a1 - a2)/e by the mean value theorem.
corner: A1 = pi - r, A2 = theta r, A3 = (1 - theta) r, r in [0, r0]; vertex-1 quantity J1/U2^2, both arcs.
Rigorous helpers: sinc, asinc = arcsin(t)/t, asinhc = arcsinh(t)/t through monotonicity in |t| with endpoint values
computed in mp.iv (arcsin p = atan2(p, sqrt(1-p^2)), arcsinh p = log(p + sqrt(p^2+1))); sinh/cosh/acos as in b4e.
usage: python b4f_mpiv_band_corner.py [e0] [r0] [prec]"""
import math
import sys
import time

import mpmath as mp
import numpy as np

e0 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.2
r0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
prec = int(sys.argv[3]) if len(sys.argv) > 3 else 80
mp.mp.prec = prec
iv = mp.iv
iv.prec = prec
FM, FP = iv.mpf(1) / 4, iv.mpf(9) / 4
PI = iv.pi
pt = lambda v: iv.mpf([v, v])
hull = lambda *xs: iv.mpf([min(x.a for x in xs), max(x.b for x in xs)])


def absrange(x):
    a, b = x.a, x.b
    lo = 0 if a <= 0 <= b else min(abs(a), abs(b))
    return lo, max(abs(a), abs(b))


def sinhI(x):
    lo = (iv.exp(pt(x.a)) - iv.exp(-pt(x.a))) / 2
    hi = (iv.exp(pt(x.b)) - iv.exp(-pt(x.b))) / 2
    return iv.mpf([lo.a, hi.b])


def coshI(x):
    m_lo, m_hi = absrange(x)
    lo = (iv.exp(pt(m_lo)) + iv.exp(-pt(m_lo))) / 2
    hi = (iv.exp(pt(m_hi)) + iv.exp(-pt(m_hi))) / 2
    return iv.mpf([max(lo.a, 1), hi.b])


def acosI(x):
    a, b = max(x.a, -1), min(x.b, 1)
    f = lambda p: iv.atan2(iv.sqrt(iv.mpf([max((1 - pt(p) ** 2).a, 0), (1 - pt(p) ** 2).b])), pt(p))
    return iv.mpf([f(b).a, f(a).b])


def atanI(x):
    return iv.mpf([iv.atan2(pt(x.a), pt(1)).a, iv.atan2(pt(x.b), pt(1)).b])


def ratio_at(f, p):
    """f(p)/p as a rigorous interval at a point p > 0; value 1 at p = 0"""
    if p == 0:
        return pt(1)
    return f(pt(p)) / pt(p)


def sincI(x):
    lo_a, hi_a = absrange(x)
    assert hi_a < math.pi
    return iv.mpf([ratio_at(iv.sin, hi_a).a, min(ratio_at(iv.sin, lo_a).b, 1)])


def asincI(x):
    lo_a, hi_a = absrange(x)
    assert hi_a < 1
    asn = lambda P: iv.atan2(P, iv.sqrt(1 - P ** 2))
    return iv.mpf([max(ratio_at(asn, lo_a).a, 1), ratio_at(asn, hi_a).b])


def asinhcI(x):
    lo_a, hi_a = absrange(x)
    ash = lambda P: iv.log(P + iv.sqrt(P ** 2 + 1))
    return iv.mpf([ratio_at(ash, hi_a).a, min(ratio_at(ash, lo_a).b, 1)])


def band(e, lam, rho):
    q = iv.mpf(1) / 3                                        # sqrt((F-1)/(F+1)) = 1/3
    he = e / 2
    y = q * he * sincI(he)                                   # q sin(e/2)
    pmh = q * sincI(he) * asincI(y)                          # psi_max / e
    psih = rho * pmh
    xp, xm = (1 + rho) * pmh / 2, (1 - rho) * pmh / 2
    bm1h = FP * (xp * sincI(e * xp)) * (xm * sincI(e * xm))  # (beta - 1)/e^2
    bm1h = iv.mpf([max(bm1h.a, 0), bm1h.b])
    beta = 1 + e ** 2 * bm1h
    Sh = iv.sqrt(bm1h * (beta + 1))
    bph = -FP * psih * sincI(e * psih)
    P = 2 * FP * iv.cos(e * psih)
    eS = e * Sh
    Dh = Sh * asinhcI(eS)
    emD = 1 / (beta + eS)
    u1 = (1 - psih) / 2
    qa = emD * u1 * sincI(e * u1)
    l2 = 2 * lam - 1
    v1, v2 = (l2 - psih) / 2, (-l2 - psih) / 2
    qb = v1 * sincI(e * v1) + emD * v2 * sincI(e * v2)
    w1 = (1 + psih) / 2
    qc = w1 * sincI(e * w1)
    disc = qb * qb + 4 * qa * qc
    disc = iv.mpf([max(disc.a, 0), disc.b])
    den = qb + iv.sqrt(disc)
    if den.a <= 0:
        raise ZeroDivisionError
    zeta = 2 * qc / den
    if zeta.a <= 0:
        raise ValueError
    L1 = iv.log(zeta)
    L2 = L1 - e * Dh
    lam2 = 1 - lam
    A1, A2 = lam * e, lam2 * e
    cA1, cA2 = iv.cos(A1), iv.cos(A2)
    ch1, ch2 = coshI(L1), coshI(L2)
    h1, h2 = ch1 + cA1, ch2 + cA2
    c1h = lam * sincI(A1) / (2 * h1)
    c2h = lam2 * sincI(A2) / (2 * h2)
    a1 = sinhI(L1) / h1
    a2 = sinhI(L2) / h2
    chx = coshI(hull(L1, L2))
    dadL = (1 + cA1 * chx) / (chx + cA1) ** 2
    dcos = sincI(he) * l2 * sincI(l2 * e / 2) / 2          # (cos A2 - cos A1)/e^2
    part2 = e * sinhI(L2) * dcos / ((ch2 + cA1) * (ch2 + cA2))
    da = Dh * dadL + part2                                  # (a1 - a2)/e
    S2 = Sh ** 2
    return ((c1h + c2h) ** 2 * (P * S2 + beta * bph ** 2) + S2 * Sh * da + 2 * bph * S2 * (a1 * c2h + a2 * c1h)
            + e * Sh * bph ** 2 * (a1 * c2h ** 2 - a2 * c1h ** 2))


def corner(r, th, rho1, sig):
    A2, A3, A1 = th * r, (1 - th) * r, PI - r
    cA1 = -iv.cos(r)
    pm = acosI((2 + FM * cA1) / FP)
    psi = rho1 * pm
    t1 = (A1 + psi) / 2
    bm1 = FP * iv.sin((1 + rho1) * pm / 2) * iv.sin((1 - rho1) * pm / 2)
    bm1 = iv.mpf([max(bm1.a, 0), bm1.b])
    beta = 1 + bm1
    Sabs = iv.sqrt(bm1 * (bm1 + 2))
    S = Sabs if sig > 0 else -Sabs
    bp = -FP * iv.sin(psi)
    P = 2 * FP * iv.cos(psi)
    eD = beta + S
    A1mt1 = (A1 - psi) / 2
    qa = eD * iv.sin(A1mt1)
    qb = eD * iv.sin(A2 + A1mt1) - iv.sin(A2 + t1)
    qc = iv.sin(t1)
    disc = qb * qb + 4 * qa * qc
    disc = iv.mpf([max(disc.a, 0), disc.b])
    den = qb + iv.sqrt(disc)
    if den.a <= 0:
        raise ZeroDivisionError
    zeta = 2 * qc / den
    if zeta.a <= 0:
        raise ValueError
    L2 = iv.log(zeta)
    sA2, cA2 = iv.sin(A2), iv.cos(A2)
    t2 = atanI(zeta * sA2 / (1 + zeta * cA2))
    s = t1 + t2
    h2 = coshI(L2) + cA2
    c2 = sA2 / (2 * h2)
    a2 = sinhI(L2) / h2
    ss, sa = iv.sin(s), iv.sin(s + A3)
    U3 = iv.sin(A3) / (ss * sa)
    U3p = -U3 * (iv.cos(s) / ss + iv.cos(s + A3) / sa)
    return ((1 - U3 * c2) ** 2 * (P * S ** 2 + beta * bp ** 2)
            + S * (U3p * (bp * c2 - S) ** 2 - a2 * (S * U3 - bp) ** 2))


def certify(f, box, scale, depth=0, maxdepth=14):
    try:
        r = f(*[iv.mpf([a, b]) for a, b in box])
        if r.a > 0:
            return True, 1, r.a
    except (ZeroDivisionError, ValueError, AssertionError):
        pass
    if depth >= maxdepth:
        return False, 1, None
    k = int(np.argmax([(b - a) / s for (a, b), s in zip(box, scale)]))
    m = 0.5 * (box[k][0] + box[k][1])
    L, R = list(box), list(box)
    L[k], R[k] = (box[k][0], m), (m, box[k][1])
    o1, n1, m1 = certify(f, tuple(L), scale, depth + 1, maxdepth)
    o2, n2, m2 = certify(f, tuple(R), scale, depth + 1, maxdepth)
    ms = [v for v in (m1, m2) if v is not None]
    return o1 and o2, n1 + n2, (min(ms) if ms else None)


def run(name, f, domain, grid):
    t0 = time.time()
    edges = [np.linspace(a, b, n + 1) for (a, b), n in zip(domain, grid)]
    scale = [b - a for a, b in domain]
    nbox, fails, leaves, refined, minlo = 0, 0, 0, 0, None
    for i in range(grid[0]):
        for j in range(grid[1]):
            for k in range(grid[2]):
                box = ((edges[0][i], edges[0][i + 1]), (edges[1][j], edges[1][j + 1]), (edges[2][k], edges[2][k + 1]))
                ok, nl, mn = certify(f, box, scale)
                nbox += 1
                leaves += nl
                refined += nl > 1
                if not ok:
                    fails += 1
                    if fails <= 5:
                        print('   %s NOT certified: %s' % (name, box))
                elif mn is not None:
                    minlo = mn if minlo is None else min(minlo, mn)
    print('%-9s %d grid boxes: refined %d, leaf boxes %d, NOT certified %d, min lower bound %s (%.0fs)'
          % (name, nbox, refined, leaves, fails, mp.nstr(minlo, 4) if minlo is not None else None, time.time() - t0))


run('band', band, [(0.0, e0), (0.0, 1.0), (-1.0, 1.0)], (8, 16, 16))
run('corner+', lambda r, th, rh: corner(r, th, rh, +1), [(0.0, r0), (0.0, 1.0), (-1.0, 1.0)], (8, 8, 16))
run('corner-', lambda r, th, rh: corner(r, th, rh, -1), [(0.0, r0), (0.0, 1.0), (-1.0, 1.0)], (8, 8, 16))
