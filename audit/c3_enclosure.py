"""AUDIT_3, item 5 (and item 2): enclosure test of the three chart functions of code/mixed_certF.py with interval F,
against the auditor's own reference audit/c0_ref.py (E-derivatives, root finding on E; independent of code/ and of
audit/b0_common.py).

Charts and domains (as in the sweeps: e0 = 0.4, r0 = 0.05; F in [13/12, 13/5]):
  main    I_main(A12, lam, rho, F)          A12 in [0.4, pi], lam in [0, 1], rho in [-1, 1]
  band    I4_band(e, lam, rho, F)           e in [0, 0.4]
  corner  J1n_corner(r, th, rho1, F, sig)   r in [0, 0.05], th in [0, 1], sig = +1 and -1
Box families per chart: 'generic' (random boxes, widths 10^U(-7,-1) of every side, F included), 'point' (zero width),
'degenerate' (shape coordinates pushed to the singular edges: main A12 -> pi, lam -> 0/1; band e -> 0, lam -> 0/1;
corner r -> 0, th -> 0/1; rho near +-1 in a third of them). The reference is evaluated at the box centre and at 3 random
points of each box (points on degenerate faces, where the triangle does not exist, are moved inside by 1e-12 relative).
A failure is a reference value outside a finite returned enclosure. Non-finite enclosures are counted separately (they
can never certify, so they are harmless).
Part 2 compares, at the Coulomb level F = 5/4 (a point F-box), mixed_certF's chart functions with the audited Coulomb
functions mixed_cert.I_iv / mixed_charts.I4_band / mixed_charts.J1n_corner on the same boxes: the enclosures must
overlap, and their lower ends should agree to rounding level.
usage: python c3_enclosure.py [n_boxes_per_family]"""
import math
import os
import sys
import time

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'code'))
import c0_ref as R  # noqa: E402
import mixed_certF  # noqa: E402  (object under test)
import mixed_cert  # noqa: E402   (audited Coulomb code, for part 2)
import mixed_charts  # noqa: E402

rng = np.random.default_rng(20261024)
FD = (13 / 12, 13 / 5)


def sample_box(dom, family, kind):
    """dom: list of (a, b); returns lo, hi lists"""
    n = len(dom)
    if family == 'point':
        w = [0.0] * n
    else:
        w = [10 ** rng.uniform(-7, -1) * (b - a) for a, b in dom]
    lo = [rng.uniform(a, b - wi) for (a, b), wi in zip(dom, w)]
    if family == 'degenerate':
        # push shape coordinates to singular edges
        if kind == 'main':
            if rng.random() < 0.6:
                lo[0] = math.pi - w[0] - 10 ** rng.uniform(-12, -2)
            if rng.random() < 0.7:
                u = 10 ** rng.uniform(-9, -2)
                lo[1] = u if rng.random() < 0.5 else 1 - u - w[1]
        elif kind == 'band':
            if rng.random() < 0.7:
                lo[0] = 10 ** rng.uniform(-8, -3)
            if rng.random() < 0.7:
                u = 10 ** rng.uniform(-9, -2)
                lo[1] = u if rng.random() < 0.5 else 1 - u - w[1]
        else:
            if rng.random() < 0.7:
                lo[0] = 10 ** rng.uniform(-8, -3)
            if rng.random() < 0.7:
                u = 10 ** rng.uniform(-9, -2)
                lo[1] = u if rng.random() < 0.5 else 1 - u - w[1]
        if rng.random() < 0.33:
            u = 10 ** rng.uniform(-9, -2)
            lo[2] = 1 - u - w[2] if rng.random() < 0.5 else -1 + u
    hi = [l + wi for l, wi in zip(lo, w)]
    # clamp into the domain
    lo = [max(a, min(l, b)) for (a, b), l in zip(dom, lo)]
    hi = [max(a, min(h, b)) for (a, b), h in zip(dom, hi)]
    hi = [max(l, h) for l, h in zip(lo, hi)]
    return lo, hi


def inside_point(lo, hi, kind, centre):
    p = [(l + h) / 2 if centre else rng.uniform(l, h) if h > l else l for l, h in zip(lo, hi)]
    # keep strictly inside the set of genuine triangles / boundary parameters
    if kind == 'main':
        p[0] = min(p[0], math.pi * (1 - 1e-15))
    if kind in ('band', 'corner'):
        p[0] = max(p[0], 1e-12)
    p[1] = min(max(p[1], 1e-12), 1 - 1e-12)
    return p


def reference(kind, p, sig):
    A, lam, rho, F = (mp.mpf(x) for x in p)
    if kind == 'main':
        return R.I_ref_main(A, lam, rho, F)
    if kind == 'band':
        return R.I4_ref_band(A, lam, rho, F)
    return R.J1n_ref_corner(A, lam, rho, F, sig)


def evaluate(kind, lo, hi, sig):
    args = [(np.array([l]), np.array([h])) for l, h in zip(lo, hi)]
    with np.errstate(all='ignore'):
        if kind == 'main':
            out = mixed_certF.I_main(*args)
        elif kind == 'band':
            out = mixed_certF.I4_band(*args)
        else:
            out = mixed_certF.J1n_corner(*args, sig)
    return float(out[0][0]), float(out[1][0])


def part1(nb):
    charts = [('main', [(0.4, math.pi), (0.0, 1.0), (-1.0, 1.0), FD], None, 50),
              ('band', [(0.0, 0.4), (0.0, 1.0), (-1.0, 1.0), FD], None, 80),
              ('corner', [(0.0, 0.05), (0.0, 1.0), (-1.0, 1.0), FD], 1, 60),
              ('corner', [(0.0, 0.05), (0.0, 1.0), (-1.0, 1.0), FD], -1, 60)]
    total_fail = 0
    for kind, dom, sig, dps in charts:
        mp.mp.dps = dps
        name = kind + ('' if sig is None else ('+' if sig > 0 else '-'))
        for family in ('generic', 'point', 'degenerate'):
            t0 = time.time()
            nfail, ntest, nonfin, noref, slack = 0, 0, 0, 0, None
            neg_lo = 0
            for _ in range(nb):
                lo, hi = sample_box(dom, family, kind)
                olo, ohi = evaluate(kind, lo, hi, sig)
                if not (np.isfinite(olo) and np.isfinite(ohi)):
                    nonfin += 1
                    continue
                if olo <= 0:
                    neg_lo += 1
                for k in range(4):
                    p = inside_point(lo, hi, kind, centre=(k == 0))
                    v = reference(kind, p, sig)
                    if v is None:
                        noref += 1
                        continue
                    ntest += 1
                    if not (mp.mpf(olo) <= v <= mp.mpf(ohi)):
                        nfail += 1
                        if nfail <= 3:
                            print('   FAIL %s %s: box %s..%s -> [%r, %r]; ref at %s = %s'
                                  % (name, family, lo, hi, olo, ohi, p, mp.nstr(v, 20)), flush=True)
                    else:
                        sl = min(v - mp.mpf(olo), mp.mpf(ohi) - v) / (abs(v) + mp.mpf(10) ** -300)
                        slack = sl if slack is None else min(slack, sl)
            total_fail += nfail
            print('%-8s %-10s %4d boxes, %5d reference points: containment failures %d; non-finite enclosures %d; '
                  'enclosures with lower end <= 0: %d; reference undefined %d; min relative slack %s  (%.0fs)'
                  % (name, family, nb, ntest, nfail, nonfin, neg_lo, noref, mp.nstr(slack, 3) if slack is not None
                     else '-', time.time() - t0), flush=True)
    return total_fail


def part2(nb):
    """Coulomb consistency: mixed_certF at the point F = 5/4 versus the audited Coulomb chart functions"""
    print('part 2: mixed_certF at F = 5/4 (point F-box) versus mixed_cert / mixed_charts (audited Coulomb code)')
    worst = {}
    for kind, dom in (('main', [(0.4, math.pi), (0.0, 1.0), (-1.0, 1.0)]),
                      ('band', [(0.0, 0.4), (0.0, 1.0), (-1.0, 1.0)]),
                      ('corner+', [(0.0, 0.05), (0.0, 1.0), (-1.0, 1.0)]),
                      ('corner-', [(0.0, 0.05), (0.0, 1.0), (-1.0, 1.0)])):
        nover, nboth, rel = 0, 0, 0.0
        for k in range(nb):
            lo, hi = sample_box(dom, 'point' if k % 3 == 0 else 'generic', kind.rstrip('+-'))
            a = [(np.array([l]), np.array([h])) for l, h in zip(lo, hi)]
            Fb = (np.array([1.25]), np.array([1.25]))
            with np.errstate(all='ignore'):
                if kind == 'main':
                    x, y = mixed_certF.I_main(*a, Fb), mixed_cert.I_iv(*a)
                elif kind == 'band':
                    x, y = mixed_certF.I4_band(*a, Fb), mixed_charts.I4_band(*a)
                else:
                    sg = 1 if kind.endswith('+') else -1
                    x, y = mixed_certF.J1n_corner(*a, Fb, sg), mixed_charts.J1n_corner(*a, sg)
            x = (float(x[0][0]), float(x[1][0]))
            y = (float(y[0][0]), float(y[1][0]))
            if not all(np.isfinite(v) for v in x + y):
                continue
            nboth += 1
            if x[0] > y[1] or y[0] > x[1]:
                nover += 1
            sc = max(abs(y[0]), abs(y[1]), 1e-300)
            rel = max(rel, abs(x[0] - y[0]) / sc)
        worst[kind] = rel
        print('   %-8s %d finite pairs: non-overlapping enclosures %d; max |lower(certF) - lower(Coulomb)| / |value| %.2e'
              % (kind, nboth, nover, rel), flush=True)


if __name__ == '__main__':
    nb = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    print('c3_enclosure: mixed_certF chart functions, F in [13/12, 13/5], %d boxes per family' % nb, flush=True)
    tf = part1(nb)
    part2(400)
    print('TOTAL containment failures: %d' % tf)
