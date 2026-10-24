"""AUDIT item 6 (approach A): direct search for two planar minima, WITHOUT using the reduction (Lemma 1-4) to count.

Charges are generated so that at least one planar minimum exists, q = Q(x0) = (lambda_i r_i^beta) with x0 in M (half of
the x0 drawn close to the fold D = c, where the image of M is thinnest and a second sheet would first appear); then
ALL critical points of V_q in T are located by root finding on the unit-scale field sum w_i d_i / sum w_i r_i
(w_i = q_i r_i^-beta) from ~450 seeds (barycentric grid + geometric clusters at the vertices), classified by the
eigenvalues of the analytic Hessian, and checked by the index count #saddles - #minima = 2 (a failure means a missed
or spurious point).  Any configuration with >= 2 minima is printed in full.
usage: python a6_direct_search.py N [alpha] [seed] [min_angle_deg]"""
import math
import sys
import time

import numpy as np
from scipy.optimize import root

N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
alpha = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 99)
amin = float(sys.argv[4]) if len(sys.argv) > 4 else 8.0
beta = 2 * alpha + 2
c = (beta - 1) / beta ** 2


def tri(A1, A2):
    A3 = math.pi - A1 - A2
    L13 = math.sin(A2) / math.sin(A3)
    return np.array([[0.0, 0.0], [1.0, 0.0], [L13 * math.cos(A1), L13 * math.sin(A1)]])


def bary(P, x):
    return np.linalg.solve(np.vstack([P.T, np.ones(3)]), np.r_[x, 1.0])


def Dval(P, x):
    lam = bary(P, x)
    z = 0j
    for i in range(3):
        d = x - P[i]
        u = complex(d[0], d[1]) / math.hypot(d[0], d[1])
        z += lam[i] * u * u
    return (1 - abs(z) ** 2) / 4


def hess(P, q, x):
    H = np.zeros((2, 2))
    for p, qi in zip(P, q):
        d = x - p
        r = math.hypot(*d)
        u = d / r
        H += -2 * alpha * qi * r ** (-beta) * (np.eye(2) - beta * np.outer(u, u))
    return H


def crit(P, q, nseed=26):
    def F(x):
        r = np.hypot(*(x - P).T)
        w = q * r ** (-beta)
        return (w[:, None] * (x - P)).sum(0) / (w * r).sum()

    seeds = []
    for i in range(1, nseed):
        for j in range(1, nseed - i):
            k = nseed - i - j
            seeds.append((i * P[0] + j * P[1] + k * P[2]) / nseed)
    for v in range(3):
        for t in np.geomspace(1e-5, 0.3, 10):
            for s in np.linspace(0.1, 0.9, 6):
                seeds.append(P[v] + t * ((1 - s) * P[(v + 1) % 3] + s * P[(v + 2) % 3] - P[v]))
    pts = []
    for s0 in seeds:
        sol = root(F, s0, method='hybr', options={'xtol': 1e-14})
        x = sol.x
        if not np.all(np.isfinite(x)) or bary(P, x).min() <= 0 or np.linalg.norm(F(x)) > 1e-12:
            continue
        if all(np.linalg.norm(x - y) > 1e-8 for y in pts):
            pts.append(x)
    return pts


t0 = time.time()
tally, parity_bad, found = {}, 0, 0
k = 0
while k < N:
    w = rng.dirichlet(np.ones(3)) * 180
    if w.min() < amin:
        continue
    P = tri(math.radians(w[0]), math.radians(w[1]))
    near_fold = k % 2 == 1
    x0 = None
    for t in range(20000):
        y = rng.dirichlet(np.ones(3)) @ P
        dv = Dval(P, y)
        if dv > c and (not near_fold or dv < c + 0.004):
            x0 = y
            break
    if x0 is None:
        continue
    k += 1
    lam0 = bary(P, x0)
    q = lam0 * np.hypot(*(x0 - P).T) ** beta
    q /= q.sum()
    pts = crit(P, q)
    kinds = []
    for x in pts:
        ev = np.linalg.eigvalsh(hess(P, q, x))
        kinds.append('min' if ev[0] > 0 else ('saddle' if ev[0] < 0 < ev[1] else 'degenerate'))
    m, s = kinds.count('min'), kinds.count('saddle')
    tally[(m, s)] = tally.get((m, s), 0) + 1
    if s - m != 2 or not any(np.linalg.norm(x - x0) < 1e-7 for x in pts):
        parity_bad += 1
    if m >= 2:
        found += 1
        print('*** m0 = %d: angles %s q %s points %s' % (m, np.round(w, 6).tolist(), q.tolist(),
                                                      [(x.tolist(), kk) for x, kk in zip(pts, kinds)]), flush=True)
    if k % 250 == 0:
        print('  %d configs, %.0fs, tally %s, parity/x0 failures %d' % (k, time.time() - t0, dict(sorted(tally.items())),
                                                                       parity_bad), flush=True)
print('alpha %.3g: %d configurations with a forced minimum (half near the fold), angles >= %.1f deg: (minima, saddles) '
      'tally %s; index-count or x0-recovery failures %d; configurations with >= 2 minima: %d  (%.0fs)'
      % (alpha, N, amin, dict(sorted(tally.items())), parity_bad, found, time.time() - t0), flush=True)
