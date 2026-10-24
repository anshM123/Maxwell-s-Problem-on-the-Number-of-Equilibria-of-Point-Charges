"""AUDIT_2 item 4b: do the three chart functions of the certificate enclose the quantity they claim?

The project's chart functions (mixed_cert.I_iv, mixed_charts.I4_band, mixed_charts.J1n_corner) are the objects under
test (imported read-only).  The reference values are computed by the auditor, independently:
  - the boundary point is found by BISECTION on the defining relation (cosh Delta = beta(s) with Delta = Lambda1 -
    Lambda2 on the line t1 + t2 = s for the vertex-3 charts; Lambda3(s) - Lambda2(t2) = +-arccosh beta1(t1) for the
    vertex-1 chart), not by the quadratic formulas of the certificate;
  - I, I/e^4 and J1/U2^2 are then evaluated from t with the auditor's own code (b0_common), 40-60 digits.
Random boxes (widths 1e-7 .. 1e-1 of the domain, plus point boxes) are evaluated by the chart function; the reference
values at random points of the box (and its centre) must lie inside the returned interval.
usage: python b4c_chart_containment.py [n_boxes_per_chart]"""
import math
import os
import sys

import mpmath as mp
import numpy as np

CODE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'code'))
sys.path.insert(0, CODE)
import mixed_cert    # noqa: E402  (object under test)
import mixed_charts  # noqa: E402  (object under test)

import b0_common as C  # noqa: E402

rng = np.random.default_rng(4040)
nb = int(sys.argv[1]) if len(sys.argv) > 1 else 200
K = mp.mpf(9) / 2
F = K / 2 - 1


def s_line_point(A, s, Dl):
    A1, A2, A3 = A
    lo, hi = max(mp.mpf(0), s - A2), min(A1, s)
    for _ in range(mp.mp.prec + 20):
        m = (lo + hi) / 2
        if C.Lam(m, A1) - C.Lam(s - m, A2) < Dl:
            lo = m
        else:
            hi = m
    t1 = (lo + hi) / 2
    return t1, s - t1


def ref_main(A12, lam, rho):
    A12, lam, rho = mp.mpf(A12), mp.mpf(lam), mp.mpf(rho)
    A1, A2 = lam * A12, (1 - lam) * A12
    A = (A1, A2, mp.pi - A12)
    pm = mp.acos((2 - (F - 1) * mp.cos(A[2])) / (F + 1))
    psi = rho * pm
    s = (A12 + psi) / 2
    beta = ((F - 1) * mp.cos(A[2]) + (F + 1) * mp.cos(psi)) / 2
    Dl = mp.acosh(max(beta, mp.mpf(1)))
    t1, t2 = s_line_point(A, s, Dl)
    return C.I_vertex3(A, K, t1, t2)[0]


def ref_band(e, lam, rho):
    return ref_main(e, lam, rho) / mp.mpf(e) ** 4


def ref_corner(r, th, rho1, sig):
    r, th, rho1 = mp.mpf(r), mp.mpf(th), mp.mpf(rho1)
    A1, A2, A3 = mp.pi - r, th * r, (1 - th) * r
    A = (A1, A2, A3)
    pm1 = mp.acos((2 + (F - 1) * mp.cos(A1)) / (F + 1))
    psi1 = rho1 * pm1
    t1 = (A1 + psi1) / 2
    beta1 = ((F + 1) * mp.cos(psi1) - (F - 1) * mp.cos(A1)) / 2
    Dv = sig * mp.acosh(max(beta1, mp.mpf(1)))
    lo, hi = mp.mpf(0), A2              # Lambda3(t1 + t2) - Lambda2(t2) decreases from +inf to -inf on (0, A2)
    for _ in range(mp.mp.prec + 20):
        m = (lo + hi) / 2
        if C.Lam(t1 + m, A1 + A2) - C.Lam(m, A2) > Dv:
            lo = m
        else:
            hi = m
    t2 = (lo + hi) / 2
    return C.J1n_vertex1(A, K, t1, t2)[0]


def run(name, f, ref, domain, dps, extra=()):
    mp.mp.dps = dps
    fails, tested, nonfinite, worst_slack = 0, 0, 0, None
    for k in range(nb):
        width = [10 ** rng.uniform(-7, -1) * (b - a) for a, b in domain] if k % 4 else [0.0] * len(domain)
        lo = [rng.uniform(a, b - w) for (a, b), w in zip(domain, width)]
        hi = [l + w for l, w in zip(lo, width)]
        args = [(np.array([l]), np.array([h])) for l, h in zip(lo, hi)]
        with np.errstate(all='ignore'):
            out = f(*args, *extra)
        olo, ohi = float(out[0][0]), float(out[1][0])
        if not (np.isfinite(olo) and np.isfinite(ohi)):
            nonfinite += 1
            continue
        pts = [[(l + h) / 2 for l, h in zip(lo, hi)]] + [[rng.uniform(l, h) if h > l else l for l, h in zip(lo, hi)]
                                                       for _ in range(4)]
        for p in pts:
            v = ref(*p, *extra)
            tested += 1
            if not (mp.mpf(olo) <= v <= mp.mpf(ohi)):
                fails += 1
                if fails <= 5:
                    print('   %s FAIL: box lo %s hi %s -> [%r, %r]; reference at %s = %s' % (name, lo, hi, olo, ohi, p, mp.nstr(v, 20)))
            slack = min(v - mp.mpf(olo), mp.mpf(ohi) - v) / (abs(v) + mp.mpf('1e-300'))
            worst_slack = slack if worst_slack is None else min(worst_slack, slack)
    print('%-10s %d boxes (a quarter of them point boxes), %d reference points: containment failures %d; '
          'non-finite enclosures %d; min relative slack %s' % (name, nb, tested, fails, nonfinite, mp.nstr(worst_slack, 3)))


run('main', mixed_cert.I_iv, ref_main, [(0.2, math.pi - 1e-9), (0.0, 1.0), (-1.0, 1.0)], 40)
run('band', mixed_charts.I4_band, ref_band, [(1e-4, 0.2), (0.0, 1.0), (-1.0, 1.0)], 60)
run('corner+', lambda r, th, rh, sg: mixed_charts.J1n_corner(r, th, rh, sg), ref_corner,
    [(1e-4, 0.05), (1e-3, 0.999), (-1.0, 1.0)], 50, extra=(1,))
run('corner-', lambda r, th, rh, sg: mixed_charts.J1n_corner(r, th, rh, sg), ref_corner,
    [(1e-4, 0.05), (1e-3, 0.999), (-1.0, 1.0)], 50, extra=(-1,))
