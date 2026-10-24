"""AUDIT item 6: re-examine the configurations where a6_direct_search.py failed its index count (large alpha).

Replays exactly the random sequence of a6_direct_search.py (same seed, same code path) and, for every configuration
whose first pass fails the index count (#saddles - #minima != 2) or misses the forced minimum x0, re-solves with:
 - many more seeds: barycentric grid n = 90, geometric clusters at each vertex (distance 1e-9 .. 0.5, 24 directions),
   points close to each edge (distance 1e-9 .. 0.2 from the edge, 40 positions along it);
 - the solver in two parametrisations (cartesian hybr on the unit-scale field, and hybr in log-barycentric
   coordinates (log l1/l3, log l2/l3), which resolves points very close to edges/vertices).
Reports the corrected (minima, saddles) and whether the index count now holds.
usage: python a6d_failures.py N alpha seed min_angle_deg"""
import math
import sys

import numpy as np
from scipy.optimize import root

N, alpha, seed, amin = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
rng = np.random.default_rng(seed)
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


def Ffield(P, q):
    def F(x):
        r = np.hypot(*(x - P).T)
        w = q * r ** (-beta)
        return (w[:, None] * (x - P)).sum(0) / (w * r).sum()
    return F


def crit_basic(P, q, nseed=26):
    F = Ffield(P, q)
    seeds = []
    for i in range(1, nseed):
        for j in range(1, nseed - i):
            k = nseed - i - j
            seeds.append((i * P[0] + j * P[1] + k * P[2]) / nseed)
    for v in range(3):
        for t in np.geomspace(1e-5, 0.3, 10):
            for s in np.linspace(0.1, 0.9, 6):
                seeds.append(P[v] + t * ((1 - s) * P[(v + 1) % 3] + s * P[(v + 2) % 3] - P[v]))
    return solve_from(P, q, F, seeds)


def solve_from(P, q, F, seeds, pts=None):
    pts = [] if pts is None else pts
    for s0 in seeds:
        sol = root(F, s0, method='hybr', options={'xtol': 1e-14})
        x = sol.x
        if not np.all(np.isfinite(x)) or bary(P, x).min() <= 0 or np.linalg.norm(F(x)) > 1e-12:
            continue
        if all(np.linalg.norm(x - y) > 1e-9 for y in pts):
            pts.append(x)
    return pts


def crit_robust(P, q):
    F = Ffield(P, q)
    seeds = []
    n = 90
    for i in range(1, n):
        for j in range(1, n - i):
            seeds.append((i * P[0] + j * P[1] + (n - i - j) * P[2]) / n)
    for v in range(3):
        for t in np.geomspace(1e-9, 0.5, 30):
            for s in np.linspace(0.02, 0.98, 24):
                seeds.append(P[v] + t * ((1 - s) * P[(v + 1) % 3] + s * P[(v + 2) % 3] - P[v]))
    for e in range(3):
        a, b, opp = P[e], P[(e + 1) % 3], P[(e + 2) % 3]
        for t in np.geomspace(1e-9, 0.2, 20):
            for s in np.linspace(0.01, 0.99, 40):
                seeds.append((1 - t) * ((1 - s) * a + s * b) + t * opp)
    pts = solve_from(P, q, F, seeds)

    # log-barycentric parametrisation
    def xb(z):
        l = np.exp(np.r_[z, 0.0])
        l /= l.sum()
        return l @ P

    G = lambda z: F(xb(z))
    for s0 in seeds[::3]:
        l = bary(P, s0)
        if l.min() <= 0:
            continue
        sol = root(G, np.log(l[:2] / l[2]), method='hybr', options={'xtol': 1e-14})
        if np.all(np.isfinite(sol.x)):
            x = xb(sol.x)
            if bary(P, x).min() > 0 and np.linalg.norm(F(x)) < 1e-12 and all(np.linalg.norm(x - y) > 1e-9 for y in pts):
                pts.append(x)
    return pts


def classify(P, q, pts):
    kinds = []
    for x in pts:
        ev = np.linalg.eigvalsh(hess(P, q, x))
        kinds.append('min' if ev[0] > 0 else ('saddle' if ev[0] < 0 < ev[1] else 'degenerate'))
    return kinds.count('min'), kinds.count('saddle')


k = 0
nfail, fixed, still = 0, 0, 0
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
    pts = crit_basic(P, q)
    m, s = classify(P, q, pts)
    if s - m == 2 and any(np.linalg.norm(x - x0) < 1e-7 for x in pts):
        continue
    nfail += 1
    pts2 = crit_robust(P, q)
    m2, s2 = classify(P, q, pts2)
    ok = s2 - m2 == 2 and any(np.linalg.norm(x - x0) < 1e-7 for x in pts2)
    fixed += ok
    still += not ok
    newpts = [x for x in pts2 if all(np.linalg.norm(x - y) > 1e-7 for y in pts)]
    print('config %d angles %s: first pass (%d,%d) -> robust (%d,%d), index count %s; points added: %s'
          % (k, np.round(w, 3).tolist(), m, s, m2, s2, 'OK' if ok else 'FAIL',
             [(np.round(bary(P, x), 10).tolist()) for x in newpts]), flush=True)
print('alpha %.3g: %d configurations replayed, %d failed the first pass; robust re-solve: %d now satisfy the index count '
      '(1 min + 3 saddles expected), %d still fail' % (alpha, N, nfail, fixed, still))
