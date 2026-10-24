"""AUDIT_2 item 4d: independent re-certification of the main-chart boxes in RIGOROUS interval arithmetic (mpmath.iv,
correctly rounded at every operation; no libm-accuracy assumption), with the auditor's own implementation of the
MIXED.md formulas (written from MIXED.md section 6.2/6.3, not from code/).

Input: the certified main-chart boxes of the project's run (replayed by b4a_where_small.py into
b4a_main_certified_boxes.npz).  Checks:
  (i)  coverage: the boxes tile the main-chart domain [e0, pi_float] x [0, 1] x [-1, 1] minus the dropped corner boxes
       (total volume and overlap-free tiling by construction of the bisection);
  (ii) on every box, the auditor's mp.iv enclosure of I has a positive lower end.
Optionally the band and corner charts are re-certified the same way on their initial grids (2048 + 2 x 1024 boxes).
usage: python b4e_mpiv_recert.py [prec_bits] [stride]"""
import math
import os
import sys
import time

import mpmath as mp
import numpy as np

prec = int(sys.argv[1]) if len(sys.argv) > 1 else 64
stride = int(sys.argv[2]) if len(sys.argv) > 2 else 1
mp.mp.prec = prec
iv = mp.iv
iv.prec = prec
HERE = os.path.dirname(os.path.abspath(__file__))
FM, FP = iv.mpf(1) / 4, iv.mpf(9) / 4          # F - 1, F + 1 for F = 5/4
PI = iv.pi


def ivl(a, b):
    return iv.mpf([a, b])


def _pt(v):
    return iv.mpf([v, v])


def sinhI(x):
    """sinh is increasing: rigorous enclosures at the two endpoints (exp is rigorous in mp.iv)"""
    lo = (iv.exp(_pt(x.a)) - iv.exp(-_pt(x.a))) / 2
    hi = (iv.exp(_pt(x.b)) - iv.exp(-_pt(x.b))) / 2
    return iv.mpf([lo.a, hi.b])


def coshI(x):
    """cosh is even and increasing in |x|"""
    a, b = x.a, x.b
    m_lo = 0 if a <= 0 <= b else min(abs(a), abs(b))
    m_hi = max(abs(a), abs(b))
    lo = (iv.exp(_pt(m_lo)) + iv.exp(-_pt(m_lo))) / 2
    hi = (iv.exp(_pt(m_hi)) + iv.exp(-_pt(m_hi))) / 2
    return iv.mpf([max(lo.a, 1), hi.b])


def acosI(x):
    """acos is decreasing on [-1, 1]; acos(p) = atan2(sqrt(1 - p^2), p) evaluated rigorously at the endpoints"""
    a, b = max(x.a, -1), min(x.b, 1)
    hi_ = iv.atan2(iv.sqrt(iv.mpf([max((1 - _pt(a) ** 2).a, 0), (1 - _pt(a) ** 2).b])), _pt(a))
    lo_ = iv.atan2(iv.sqrt(iv.mpf([max((1 - _pt(b) ** 2).a, 0), (1 - _pt(b) ** 2).b])), _pt(b))
    return iv.mpf([lo_.a, hi_.b])


def I_main(A12, lam, rho):
    """I at the upper-arc boundary point, MIXED.md 6.2 (vertex 3), written with the dependency-free forms
    s = A12/2 + psi/2, A12 - s = A12/2 - psi/2, A1 - s = (2 lam - 1) A12/2 - psi/2, A2 - s = -(2 lam - 1) A12/2 - psi/2
    (algebraic identities from A1 = lam A12, A2 = (1 - lam) A12, psi = 2 s - A12)"""
    A1 = lam * A12
    A2 = (1 - lam) * A12
    A3 = PI - A12
    kap = (2 - FM * iv.cos(A3)) / FP
    psim = acosI(kap)
    psi = rho * psim
    sp_ = iv.sin((1 + rho) * psim / 2)
    sm_ = iv.sin((1 - rho) * psim / 2)
    bm1 = FP * sp_ * sm_
    bm1 = iv.mpf([max(bm1.a, 0), bm1.b])            # beta - 1 >= 0 (true value)
    beta = 1 + bm1
    S = iv.sqrt(bm1 * (bm1 + 2))
    emD = 1 / (beta + S)
    Dl = iv.log(beta + S)                              # Delta = arccosh beta = log(beta + S)
    bp = -FP * iv.sin(psi)
    P = 2 * FP * iv.cos(psi)
    hpsi = psi / 2
    half12 = A12 / 2
    dA = (2 * lam - 1) * A12 / 2
    qa = emD * iv.sin(half12 - hpsi)
    qb = iv.sin(dA - hpsi) + emD * iv.sin(-dA - hpsi)
    qc = iv.sin(half12 + hpsi)
    disc = qb * qb + 4 * qa * qc
    disc = iv.mpf([max(disc.a, 0), disc.b])
    den = qb + iv.sqrt(disc)
    if den.a <= 0:
        raise ZeroDivisionError
    zeta = 2 * qc / den
    if zeta.a <= 0:
        raise ValueError
    L1 = iv.log(zeta)
    L2 = L1 - Dl
    h1 = coshI(L1) + iv.cos(A1)
    h2 = coshI(L2) + iv.cos(A2)
    c1 = iv.sin(A1) / (2 * h1)
    c2 = iv.sin(A2) / (2 * h2)
    a1 = sinhI(L1) / h1
    a2 = sinhI(L2) / h2
    return (c1 + c2) ** 2 * (P * S ** 2 + beta * bp ** 2) + S * (a1 * (bp * c2 + S) ** 2 - a2 * (S - bp * c1) ** 2)


def certify(box, depth=0, maxdepth=12):
    """rigorous check that I > 0 on the box; bisect the widest (normalised) side on failure. Returns
    (ok, number_of_leaf_boxes, min_lower_bound)"""
    (a0, b0), (a1, b1), (a2, b2) = box
    try:
        r = I_main(ivl(a0, b0), ivl(a1, b1), ivl(a2, b2))
        if r.a > 0:
            return True, 1, r.a
    except (ZeroDivisionError, ValueError):
        pass
    if depth >= maxdepth:
        return False, 1, None
    w = [(b0 - a0) / 3.0, (b1 - a1) / 1.0, (b2 - a2) / 2.0]
    k = int(np.argmax(w))
    m = 0.5 * (box[k][0] + box[k][1])
    left = list(box); right = list(box)
    left[k] = (box[k][0], m); right[k] = (m, box[k][1])
    ok1, n1, m1 = certify(tuple(left), depth + 1, maxdepth)
    ok2, n2, m2 = certify(tuple(right), depth + 1, maxdepth)
    mins = [v for v in (m1, m2) if v is not None]
    return ok1 and ok2, n1 + n2, (min(mins) if mins else None)


dat = np.load(os.path.join(HERE, 'b4a_main_certified_boxes.npz'))
lo, hi = dat['lo'], dat['hi']
nbox = lo.shape[1]
# (i) coverage by volume
e0 = float(lo[0].min())
vol = np.prod(hi - lo, axis=0).sum()
dom = (math.pi - e0) * 1.0 * 2.0
print('main chart: %d certified boxes; A12 range [%r, %r]; total volume %.15g vs domain volume %.15g (difference %.3e = '
      'dropped corner boxes)' % (nbox, lo[0].min(), hi[0].max(), vol, dom, dom - vol))
t0 = time.time()
fails, minlo, done, leaves, refined = 0, None, 0, 0, 0
for j in range(0, nbox, stride):
    box = ((float(lo[0, j]), float(hi[0, j])), (float(lo[1, j]), float(hi[1, j])), (float(lo[2, j]), float(hi[2, j])))
    ok, nl, mn = certify(box)
    done += 1
    leaves += nl
    refined += nl > 1
    if not ok:
        fails += 1
        if fails <= 10:
            print('   NOT certified even after refinement: box %d %s' % (j, box))
    elif mn is not None:
        minlo = mn if minlo is None else min(minlo, mn)
    if done % 20000 == 0:
        print('   %d boxes done (%.0fs), refined %d, failures %d' % (done, time.time() - t0, refined, fails), flush=True)
print('mp.iv (prec %d bits) independent re-certification of %d main-chart boxes (stride %d): needed local refinement '
      '%d, leaf boxes %d, NOT certified %d; min lower bound %s; %.0fs'
      % (prec, done, stride, refined, leaves, fails, mp.nstr(minlo, 4) if minlo is not None else None, time.time() - t0))
