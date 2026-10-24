"""Containment test for the general-exponent certificate mixed_certF.py: the analogue of audit/b4c_chart_containment.py
(which tests the Coulomb chart functions) with F = K/2 - 1 an interval parameter, F in [F(1/4), F(2)] = [13/12, 13/5].

The reference values come from the independent implementation in audit/b0_common.py: the boundary point is found by
bisection on the defining relation (not by the closed forms of the certificate), and I, I/e^4, J1/U2^2 are evaluated from
t at 40-60 digits with K = 2(F + 1). Random boxes (widths 1e-7 .. 1e-1 of the domain in every coordinate, F included, plus
point boxes) are evaluated by the chart functions of mixed_certF.py; the reference values at the centre and at random
points of each box must lie inside the returned interval.
usage: python check_certF_containment.py [n_boxes_per_chart]"""
import math
import os
import sys

import mpmath as mp
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'audit')))
import b0_common as C  # noqa: E402  (independent reference implementation)
import mixed_certF     # noqa: E402  (object under test)

rng = np.random.default_rng(2026)
F_DOM = (13 / 12, 13 / 5)


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


def ref_main(A12, lam, rho, F):
    A12, lam, rho, F = mp.mpf(A12), mp.mpf(lam), mp.mpf(rho), mp.mpf(F)
    A1, A2 = lam * A12, (1 - lam) * A12
    A = (A1, A2, mp.pi - A12)
    pm = mp.acos((2 - (F - 1) * mp.cos(A[2])) / (F + 1))
    psi = rho * pm
    s = (A12 + psi) / 2
    beta = ((F - 1) * mp.cos(A[2]) + (F + 1) * mp.cos(psi)) / 2
    Dl = mp.acosh(max(beta, mp.mpf(1)))
    t1, t2 = s_line_point(A, s, Dl)
    return C.I_vertex3(A, 2 * (F + 1), t1, t2)[0]


def ref_band(e, lam, rho, F):
    return ref_main(e, lam, rho, F) / mp.mpf(e) ** 4


def ref_corner(r, th, rho1, F, sig):
    r, th, rho1, F = mp.mpf(r), mp.mpf(th), mp.mpf(rho1), mp.mpf(F)
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
    return C.J1n_vertex1(A, 2 * (F + 1), t1, t2)[0]


def run(name, f, ref, domain, dps, nb, extra=()):
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
                    print('   %s FAIL: box lo %s hi %s -> [%r, %r]; reference at %s = %s'
                          % (name, lo, hi, olo, ohi, p, mp.nstr(v, 20)))
            slack = min(v - mp.mpf(olo), mp.mpf(ohi) - v) / (abs(v) + mp.mpf('1e-300'))
            worst_slack = slack if worst_slack is None else min(worst_slack, slack)
    print('%-10s %d boxes (a quarter of them point boxes), %d reference points: containment failures %d; '
          'non-finite enclosures %d; min relative slack %s' % (name, nb, tested, fails, nonfinite, mp.nstr(worst_slack, 3)),
          flush=True)


if __name__ == '__main__':
    nb = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    run('main', mixed_certF.I_main, ref_main, [(0.2, math.pi - 1e-9), (0.0, 1.0), (-1.0, 1.0), F_DOM], 40, nb)
    run('band', mixed_certF.I4_band, ref_band, [(1e-4, 0.4), (0.0, 1.0), (-1.0, 1.0), F_DOM], 60, nb)
    run('corner+', lambda r, th, rh, F, sg: mixed_certF.J1n_corner(r, th, rh, F, sg), ref_corner,
        [(1e-4, 0.05), (1e-3, 0.999), (-1.0, 1.0), F_DOM], 50, nb, extra=(1,))
    run('corner-', lambda r, th, rh, F, sg: mixed_certF.J1n_corner(r, th, rh, F, sg), ref_corner,
        [(1e-4, 0.05), (1e-3, 0.999), (-1.0, 1.0), F_DOM], 50, nb, extra=(-1,))
