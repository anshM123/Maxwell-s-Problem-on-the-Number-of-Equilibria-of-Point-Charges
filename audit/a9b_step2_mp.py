"""AUDIT: is the slightly negative float64 minimum of the Step-2 quantity in a9_collinear_step2.py real or rounding?

m(theta, theta') = sum_i (f_i(theta_i) - f_i(theta'_i)) (theta_i - theta'_i) / |theta - theta'|^2 on boundary pairs of
{E <= (b+1)^2/b}.  Per shape: float64 Nelder-Mead minimisation of m over pairs of boundary angles (as in a9), then the
minimiser is re-evaluated with 50-digit arithmetic (boundary points recomputed by 50-digit bisection), together with
the separation |theta - theta'| and the tangential value t^T diag(d) t / |t|^2 at the boundary point (the limit of m
for coincident pairs; it vanishes exactly where the null direction of the form is tangent to the boundary).
usage: python a9b_step2_mp.py [nshapes]"""
import math
import sys

import mpmath as mp
import numpy as np
from scipy.optimize import minimize

mp.mp.dps = 50
rng = np.random.default_rng(31337)
nsh = int(sys.argv[1]) if len(sys.argv) > 1 else 40


def make(A, b, lib):
    k = (b + 1) ** 2 / b
    if lib is np:
        sin, log = np.sin, np.log
    else:
        sin, log = mp.sin, mp.log

    def E(th):
        W = [sin(A[i]) / (sin(th[i] - A[i]) * sin(th[i])) for i in range(3)]
        return W[0] * W[1] + W[1] * W[2] + W[2] * W[0]

    def F(th):
        return [log(sin(th[i])) + b * log(sin(th[i] - A[i])) for i in range(3)]

    th0 = [(lib.pi + A[i]) / 2 for i in range(3)]

    def bpoint(phi, iters):
        e = (lib.cos(phi), lib.sin(phi))
        lo, hi = 0 * k, k * 0 + 1e-9

        def th_at(s):
            t1, t2 = th0[0] + s * e[0], th0[1] + s * e[1]
            return [t1, t2, 2 * lib.pi - t1 - t2]

        def inside(s):
            t = th_at(s)
            return all(A[i] < t[i] < lib.pi for i in range(3)) and E(t) <= k

        while inside(hi):
            lo, hi = hi, 2 * hi
        for _ in range(iters):
            m = (lo + hi) / 2
            if inside(m):
                lo = m
            else:
                hi = m
        return th_at(lo)

    def mono(t1, t2):
        f1, f2 = F(t1), F(t2)
        v = [t1[i] - t2[i] for i in range(3)]
        return sum((f1[i] - f2[i]) * v[i] for i in range(3)) / sum(vi * vi for vi in v), math.sqrt(float(sum(vi * vi for vi in v)))

    return k, E, F, bpoint, mono


rows = []
for s in range(nsh):
    while True:
        Ad = rng.dirichlet(np.ones(3)) * 180
        if Ad.min() > [0.5, 5, 20][s % 3]:
            break
    alpha = [0.05, 0.25, 0.5, 1.0, 2.5, 10.0][s % 6]
    b = 2 * alpha + 1
    A = np.radians(Ad)
    k, E, F, bpoint, mono = make(A, b, np)
    phis = rng.uniform(0, 2 * math.pi, 120)
    B = [bpoint(p, 70) for p in phis]
    cand = sorted(((mono(B[i], B[j])[0], phis[i], phis[j]) for i in range(120) for j in range(i + 1, 120)))[:2]
    med = np.median([mono(B[i], B[(i + 1) % 120])[0] for i in range(120)])
    best = None
    for v0, p1, p2 in cand:
        res = minimize(lambda z: mono(bpoint(z[0], 70), bpoint(z[1], 70))[0] if abs(z[0] - z[1]) > 1e-7 else 1e9,
                       [p1, p2], method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-14, 'maxiter': 500})
        if best is None or res.fun < best[0]:
            best = (res.fun, res.x)
    # 50-digit re-evaluation at the minimiser
    Am = [mp.mpf(float(a)) * 1 for a in A]
    Am = [mp.mpf(float(Ad[0])) * mp.pi / 180, mp.mpf(float(Ad[1])) * mp.pi / 180]
    Am.append(mp.pi - Am[0] - Am[1])
    km, Em, Fm, bpm, monom = make(Am, 2 * mp.mpf(alpha) + 1, mp)
    t1 = bpm(mp.mpf(float(best[1][0])), 200)
    t2 = bpm(mp.mpf(float(best[1][1])), 200)
    mm, sep = monom(t1, t2)
    # tangential limit at t1: tangent ~ t2 - t1 direction; form value with d_i = cot th + b cot(th - A)
    bm = 2 * mp.mpf(alpha) + 1
    d = [mp.cot(t1[i]) + bm * mp.cot(t1[i] - Am[i]) for i in range(3)]
    v = [t2[i] - t1[i] for i in range(3)]
    tang = sum(d[i] * v[i] ** 2 for i in range(3)) / sum(vi * vi for vi in v)
    rows.append((s, Ad.round(3).tolist(), alpha, best[0] / med, float(mm / med), sep, float(tang / med)))
    print('shape %2d angles %s alpha %5.2f: float64 min m/median %+.3e | 50-digit m/median at that pair %+.3e '
          '(separation %.1e, tangential form/median %+.3e)' % rows[-1], flush=True)
neg64 = sum(1 for r in rows if r[3] < 0)
negmp = sum(1 for r in rows if r[4] <= 0)
print('SUMMARY: %d shapes; float64 minimum negative for %d; 50-digit value at the same pair <= 0 for %d; '
      'min 50-digit m/median %.3e' % (len(rows), neg64, negmp, min(r[4] for r in rows)))
