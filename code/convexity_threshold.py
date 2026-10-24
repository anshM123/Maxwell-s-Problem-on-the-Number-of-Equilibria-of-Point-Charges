"""Numerical (not rigorous) location of the exponent where Claim M stops holding.

(1) scan: for each alpha, random triangles and boundary points, with extra weight near the degenerate corners where
    the minimum lives. The scale-free quantity I/(A1+A2)^4 is evaluated with the main chart of mixed_certF at point
    boxes, which gives tight values. (The band chart is not used here, because its enclosures are wide even for point
    boxes.) The smallest value found is re-evaluated with the independent reference of check_certF_containment: the
    boundary point is found by bisection, and I uses the formulas of audit/b0_common, at 40 digits.
(2) corner: at alpha = 5.6, 5.7, 5.8, the reference value at 50 digits along the degenerate corner A1 + A2 -> 0,
    A2/A1 -> 0. The minimum over the boundary parameter rho is taken on a grid.
usage: python convexity_threshold.py [alpha ...]"""
import math
import sys

import mpmath as mp
import numpy as np

from check_certF_containment import ref_main
from mixed_certF import I_main, _F_of_alpha

rng = np.random.default_rng(12)


def pt(x):
    return (x, x.copy())


def scan(F, n=1_600_000, batch=200_000):
    best = (np.inf, None)
    for _ in range(0, n, batch):
        near = rng.random(batch) < 0.5
        A = np.where(near, 10 ** rng.uniform(-2.7, -0.7, batch), rng.uniform(0.02, math.pi - 1e-6, batch))
        u = np.where(near, 10 ** rng.uniform(-7, -1, batch), 10 ** rng.uniform(-5, 0, batch))
        lam = np.where(rng.random(batch) < 0.5, u, 1 - u)
        rho = rng.uniform(-1, 1, batch)
        with np.errstate(all='ignore'):
            v = I_main(pt(A), pt(lam), pt(rho), pt(np.full(batch, F)))[0] / A ** 4
        v = np.where(np.isfinite(v), v, np.inf)
        i = int(np.argmin(v))
        if v[i] < best[0]:
            best = (v[i], (A[i], lam[i], rho[i]))
    return best


if __name__ == '__main__':
    alphas = [float(a) for a in sys.argv[1:]] or [0.25, 2, 4, 5, 5.2, 5.4, 5.5, 5.6, 5.8, 6, 8, 10]
    mp.mp.dps = 40
    print('(1) minimum of I/(A1+A2)^4 over sampled triangles and boundary points')
    for a in alphas:
        F = _F_of_alpha(a)
        v, p = scan(F)
        ref = ref_main(*p, F) / mp.mpf(p[0]) ** 4
        A1, A2 = p[0] * p[1], p[0] * (1 - p[1])
        print('alpha %5.2f: %+.4e  [40-digit reference %+.4e]  at angles (%.3f, %.3f, %.3f) deg, rho %+.4f'
              % (a, v, float(ref), math.degrees(A1), math.degrees(A2), 180 - math.degrees(p[0]), p[2]), flush=True)
    mp.mp.dps = 50
    print('(2) the degenerate corner A1 + A2 = e -> 0 with A2/A1 = 1e-5 (50 digits; min over rho in [0.80, 0.94])')
    for a in (5.6, 5.7, 5.8):
        F = _F_of_alpha(a)
        row = []
        for e in (3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4):
            row.append(min(ref_main(e, 1 - 1e-5, 0.80 + 0.005 * k, F) for k in range(29)) / mp.mpf(e) ** 4)
        print('alpha %.1f: ' % a + '  '.join('e=%.0e %+.4e' % (e, float(v)) for e, v in
                                             zip((3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4), row)), flush=True)
