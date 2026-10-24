"""AUDIT_2 item 6: 50-digit re-examination of the one float64 'candidate' with 2 extrema printed by the first (crashed)
run of b6a_direct_search.py:  angles (0.0373164, 0.0262531, 3.07802) rad (printed to 6 digits),
q = [0.2921064821412595, 0.20773284746969184, -0.5001606703890485] (full precision).
All critical points are located by float root finding from a dense seed set, then every candidate is refined by
50-digit Newton on grad V = 0 and classified by its 50-digit Hessian; the index count is checked.
usage: python b6d_candidate_check.py"""
import math

import mpmath as mp
import numpy as np
from scipy.optimize import root

mp.mp.dps = 50
A1, A2 = 0.0373164, 0.0262531
q = np.array([0.2921064821412595, 0.20773284746969184, -0.5001606703890485])
A3 = math.pi - A1 - A2
L13 = math.sin(A2) / math.sin(A3)
P = np.array([[0.0, 0.0], [1.0, 0.0], [L13 * math.cos(A1), L13 * math.sin(A1)]])


def G(x):
    r = np.linalg.norm(x[None] - P, axis=1)
    w = q * r ** -3.0
    return (w[:, None] * (x[None] - P)).sum(0) / (np.abs(w) * r).sum()


cen = P.mean(0)
seeds = []
for rr in np.geomspace(1e-4, 200, 60):
    for th in np.linspace(0, 2 * math.pi, 64, endpoint=False):
        seeds.append(cen + rr * np.array([math.cos(th), math.sin(th)]))
for v in range(3):
    for rr in np.geomspace(1e-7, 0.5, 20):
        for th in np.linspace(0, 2 * math.pi, 24, endpoint=False):
            seeds.append(P[v] + rr * np.array([math.cos(th), math.sin(th)]))
cands = []
for s0 in seeds:
    sol = root(G, s0, method='hybr', options={'xtol': 1e-15})
    x = sol.x
    if np.all(np.isfinite(x)) and np.linalg.norm(x - cen) < 1e5 and np.linalg.norm(G(x)) < 1e-9:
        if np.min(np.linalg.norm(x[None] - P, axis=1)) > 1e-12 and all(np.linalg.norm(x - y) > 1e-10 for y in cands):
            cands.append(x)
print('float candidates (|G| < 1e-9): %d' % len(cands))
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
for x in cands:
    xm = [mp.mpf(float(x[0])), mp.mpf(float(x[1]))]
    ok = False
    for it in range(80):
        g, H = gH(xm)
        try:
            st = mp.lu_solve(H, mp.matrix(g))
        except ZeroDivisionError:
            break
        xm = [xm[0] - st[0], xm[1] - st[1]]
        if mp.sqrt(st[0] ** 2 + st[1] ** 2) < mp.mpf('1e-42'):
            ok = True
            break
    if ok and all(abs(xm[0] - y[0]) + abs(xm[1] - y[1]) > mp.mpf('1e-30') for y, _, _ in good):
        g, H = gH(xm)
        ev = mp.eigsy(H)[0]
        kind = 'extremum' if ev[0] * ev[1] > 0 else 'saddle'
        good.append((xm, kind, ev))
for xm, kind, ev in good:
    print('  x = (%s, %s): %s, Hessian eigenvalues %s, %s' % (mp.nstr(xm[0], 15), mp.nstr(xm[1], 15), kind,
                                                             mp.nstr(ev[0], 6), mp.nstr(ev[1], 6)))
ne = sum(1 for _, k, _ in good if k == 'extremum')
ns = sum(1 for _, k, _ in good if k == 'saddle')
print('50-digit verified: %d extrema, %d saddles (Q_tot = %.3e; index count requires saddles - extrema = %d)'
      % (ne, ns, q.sum(), 2 if abs(q.sum()) > 0 else 1))
