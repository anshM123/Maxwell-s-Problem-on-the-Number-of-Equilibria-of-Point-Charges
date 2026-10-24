"""AUDIT_2 item 6 (direct): search for mixed-sign Coulomb configurations with two minima or more than four critical
points, WITHOUT using the reduction to count.

Each configuration has a minimum forced: a point t of M' = {E < K} (rejection sampling in Q) is mapped to x0 in R+out
by the auditor's geometric inverse (inversion at p3), and q = (Q1, Q2, Q3)(x0) with Q_i = lambda_i r_i^3 (q3 < 0).
Shapes: uniformly random (half), and near-degenerate families (half): A1 + A2 small (down to 1e-3 rad), A1 near pi,
A3 near 0, A1 or A2 near 0.  All critical points in the plane are located by root finding on the bounded field
sum w_i d_i / sum |w_i| r_i (w_i = q_i r_i^-3) from ~1100 seeds (log-polar grid to 80 x the triangle size and clusters
at the vertices), classified by the Hessian, and checked by the index count (#saddles - #extrema = 2, or 1 if
Q_tot = 0) and by recovery of x0.  Any configuration with >= 2 extrema or > 4 critical points is printed.
usage: python b6a_direct_search.py N [seed]"""
import math
import sys
import time

import numpy as np
from scipy.optimize import root

N = int(sys.argv[1]) if len(sys.argv) > 1 else 500
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 6)
K = 4.5


def tri(A1, A2):
    A3 = math.pi - A1 - A2
    L13 = math.sin(A2) / math.sin(A3)
    return np.array([[0.0, 0.0], [1.0, 0.0], [L13 * math.cos(A1), L13 * math.sin(A1)]])


def Ufun(t, a):
    return 1 / math.tan(t) + 1 / math.tan(a - t)


def E_of_t(A1, A2, t1, t2):
    U1, U2, U3 = Ufun(t1, A1), Ufun(t2, A2), Ufun(t1 + t2, A1 + A2)
    return U1 * U2 - U3 * (U1 + U2)


def x_of_t(P, t1, t2):
    p = P[:, 0] + 1j * P[:, 1]
    d1, d2 = p[0] - p[2], p[1] - p[2]
    a = d1 * np.exp(-1j * t1)
    bb = -d2 * np.exp(1j * t2)
    rhs = d1 - d2
    r1, r2 = np.linalg.solve(np.array([[a.real, bb.real], [a.imag, bb.imag]]), np.array([rhs.real, rhs.imag]))
    Y = (1 - r1 * np.exp(-1j * t1)) / d2
    z = p[2] + 1 / Y
    return np.array([z.real, z.imag]), (r1 > 0 and r2 > 0)


def bary(P, x):
    return np.linalg.solve(np.vstack([P.T, np.ones(3)]), np.r_[x, 1.0])


def hess(P, q, x):
    H = np.zeros((2, 2))
    for p, qi in zip(P, q):
        d = x - p
        r = np.linalg.norm(d)
        u = d / r
        H += -qi * r ** -3 * (np.eye(2) - 3 * np.outer(u, u))
    return H


def crit(P, q):
    def G(x):
        r = np.linalg.norm(x[None] - P, axis=1)
        w = q * r ** -3.0
        return (w[:, None] * (x[None] - P)).sum(0) / (np.abs(w) * r).sum()

    cen = P.mean(0)
    sc = np.max(np.linalg.norm(P - cen, axis=1))
    seeds = []
    for rr in np.geomspace(0.02, 80, 26) * sc:
        for th in np.linspace(0, 2 * math.pi, 32, endpoint=False):
            seeds.append(cen + rr * np.array([math.cos(th), math.sin(th)]))
    for v in range(3):
        for rr in np.geomspace(1e-6, 0.5, 12) * sc:
            for th in np.linspace(0, 2 * math.pi, 12, endpoint=False):
                seeds.append(P[v] + rr * np.array([math.cos(th), math.sin(th)]))
    out = []
    for s0 in seeds:
        sol = root(G, s0, method='hybr', options={'xtol': 1e-15})
        x = sol.x
        if not np.all(np.isfinite(x)) or np.linalg.norm(x - cen) > 1e4 * sc:
            continue
        if np.linalg.norm(G(x)) > 1e-11 or np.min(np.linalg.norm(x[None] - P, axis=1)) < 1e-12 * sc:
            continue
        if all(np.linalg.norm(x - y) > 1e-9 * sc for y in out):
            out.append(x)
    return out


def rand_shape(k):
    fam = k % 10
    if fam < 5:
        while True:
            w = rng.dirichlet(np.ones(3)) * math.pi
            if w.min() > math.radians(0.5):
                return w[0], w[1], 'random'
    e = 10 ** rng.uniform(-3, -0.7)
    if fam == 5:                                        # A1 + A2 small (band)
        lam = rng.uniform(0.02, 0.98)
        return lam * e, (1 - lam) * e, 'A12 small'
    if fam == 6:                                        # A1 near pi (corner)
        th = rng.uniform(0.02, 0.98)
        return math.pi - e, th * e, 'A1 near pi'
    if fam == 7:                                        # A3 near 0
        lam = rng.uniform(0.05, 0.95)
        return lam * (math.pi - e), (1 - lam) * (math.pi - e), 'A3 near 0'
    if fam == 8:                                        # A1 near 0
        A2 = rng.uniform(0.2, math.pi - 0.2)
        return e, A2 - e if A2 - e > 0.05 else A2, 'A1 near 0'
    th = rng.uniform(0.02, 0.98)                         # A2 near pi
    return th * e, math.pi - e, 'A2 near pi'


import mpmath as mp  # noqa: E402
mp.mp.dps = 50


def mp_verify(P, q, pts):
    """refine every float candidate by 50-digit Newton on grad V = 0 (analytic Hessian); keep the converged ones,
    deduplicate, classify by the 50-digit Hessian eigenvalues"""
    Pm = [[mp.mpf(float(v)) for v in p] for p in P]
    qm = [mp.mpf(float(v)) for v in q]

    def gH(x):
        g = [mp.mpf(0), mp.mpf(0)]
        H = mp.matrix(2, 2)
        for p, qi in zip(Pm, qm):
            d = [x[0] - p[0], x[1] - p[1]]
            r2 = d[0] ** 2 + d[1] ** 2
            r = mp.sqrt(r2)
            for a in range(2):
                g[a] += -qi * d[a] / r ** 3
                for c in range(2):
                    H[a, c] += -qi / r ** 3 * ((1 if a == c else 0) - 3 * d[a] * d[c] / r2)
        return g, H

    good = []
    for x in pts:
        xm = [mp.mpf(float(x[0])), mp.mpf(float(x[1]))]
        ok = False
        for it in range(60):
            g, H = gH(xm)
            try:
                st = mp.lu_solve(H, mp.matrix(g))
            except ZeroDivisionError:
                break
            xm = [xm[0] - st[0], xm[1] - st[1]]
            if mp.sqrt(st[0] ** 2 + st[1] ** 2) < mp.mpf('1e-40'):
                ok = True
                break
        if not ok:
            continue
        g, H = gH(xm)
        if all(abs(xm[0] - y[0]) + abs(xm[1] - y[1]) > mp.mpf('1e-30') for y, _ in good):
            ev = mp.eigsy(H)[0]
            kind = 'ext' if (ev[0] > 0 and ev[1] > 0) or (ev[0] < 0 and ev[1] < 0) else ('sad' if ev[0] * ev[1] < 0 else 'deg')
            good.append((xm, kind))
    return good


t0 = time.time()
tally, bad, found = {}, [], 0
n_verified = 0
k = 0
while k < N:
    A1, A2, fam = rand_shape(k)
    A3 = math.pi - A1 - A2
    if min(A1, A2, A3) <= 0:
        continue
    P = tri(A1, A2)
    t = None
    for _ in range(4000):
        t1, t2 = A1 * rng.uniform(0.01, 0.99), A2 * rng.uniform(0.01, 0.99)
        if E_of_t(A1, A2, t1, t2) < K:
            t = (t1, t2)
            break
    if t is None:
        t = (A1 / 2, A2 / 2)
    x0, okb = x_of_t(P, *t)
    if not okb:
        continue
    lam = bary(P, x0)
    r = np.linalg.norm(x0[None] - P, axis=1)
    q = lam * r ** 3
    q = q / np.abs(q).sum()
    if not (q[0] > 0 and q[1] > 0 and q[2] < 0):
        continue
    k += 1
    pts = crit(P, q)
    idx = []
    for x in pts:
        ev = np.linalg.eigvalsh(hess(P, q, x))
        idx.append('ext' if (ev[0] > 0 or ev[1] < 0) else ('sad' if ev[0] < 0 < ev[1] else 'deg'))
    m, s = idx.count('ext'), idx.count('sad')
    Qt = q.sum()
    want = 1 if abs(Qt) < 1e-14 else 2
    got_x0 = any(np.linalg.norm(x - x0) < 1e-7 * np.max(np.linalg.norm(P, axis=1)) for x in pts)
    if s - m != want or not got_x0 or m >= 2 or m + s > 4:
        n_verified += 1
        good = mp_verify(P, q, pts)
        m2, s2 = sum(1 for _, kk in good if kk == 'ext'), sum(1 for _, kk in good if kk == 'sad')
        print('  float64 suspicious (family %s, angles %r, %r, q %r): float (%d ext, %d sad) -> 50-digit verified '
              '(%d ext, %d sad, %d deg)' % (fam, A1, A2, q.tolist(), m, s, m2, s2,
                                            sum(1 for _, kk in good if kk == 'deg')), flush=True)
        m, s = m2, s2
    key = (fam, m, s)
    tally[key] = tally.get(key, 0) + 1
    if s - m != want or not got_x0:
        bad.append((fam, A1, A2, q.tolist(), m, s, got_x0, abs(Qt)))
    if m >= 2 or m + s > 4:
        found += 1
        print('*** CANDIDATE (after 50-digit verification): family %s angles %r, %r q %r: %d extrema, %d saddles'
              % (fam, A1, A2, q.tolist(), m, s), flush=True)
print('%d configurations in %.0fs' % (N, time.time() - t0))
for key in sorted(tally):
    print('   family %-12s (extrema, saddles) = (%d, %d): %d' % (key[0], key[1], key[2], tally[key]))
print('configurations re-verified in 50 digits: %d; index-count or x0-recovery failures after verification: %d; '
      'candidates with >= 2 extrema or > 4 critical points after verification: %d' % (n_verified, len(bad), found))
for b_ in bad[:15]:
    print('   failure: family %s angles (%.4g, %.4g) q %s -> (%d, %d), x0 found %s, |Q_tot| %.2e' % b_)
