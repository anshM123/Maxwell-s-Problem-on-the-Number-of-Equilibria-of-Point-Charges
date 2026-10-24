"""AUDIT_2 item 5 (collinear case, improved test): on-axis equilibria are counted with a far-field scan (log-spaced up to
1e7 beyond the outer charges), and the perturbation step of MIXED.md 7(b) is tested IN THE PLANE of the perturbed
charges (part (a): all critical points of a non-collinear configuration lie in that plane) with a dense planar search.
Checks: (i) number of on-axis equilibria (claim: 0 or 2), (ii) every non-degenerate on-axis equilibrium has a nearby
critical point after moving p3 off the axis by eps (implicit function theorem), (iii) the perturbed configuration has
<= 4 critical points and satisfies the index count.
usage: python b5b_collinear_fixed.py [n]"""
import math
import sys

import numpy as np
from scipy.optimize import brentq, root

rng = np.random.default_rng(58)
n = int(sys.argv[1]) if len(sys.argv) > 1 else 150


def dV1(z, zs, q):
    return sum(-qi * np.sign(z - zi) / (z - zi) ** 2 for zi, qi in zip(zs, q))


def axis_equilibria(zs, q):
    out = []
    span = zs[2] - zs[0]
    grids = [zs[0] - span * np.geomspace(1e7, 1e-9, 20000), np.linspace(zs[0], zs[1], 20001)[1:-1],
             np.linspace(zs[1], zs[2], 20001)[1:-1], zs[2] + span * np.geomspace(1e-9, 1e7, 20000)]
    for gz in grids:
        v = np.array([dV1(z, zs, q) for z in gz])
        for j in np.nonzero(np.sign(v[1:]) != np.sign(v[:-1]))[0]:
            out.append(brentq(lambda z: dV1(z, zs, q), gz[j], gz[j + 1], xtol=1e-14, rtol=1e-15))
    return out


def crit_plane(P, q):
    def G(x):
        r = np.linalg.norm(x[None] - P, axis=1)
        w = q * r ** -3.0
        return (w[:, None] * (x[None] - P)).sum(0) / (np.abs(w) * r).sum()
    cen = P.mean(0)
    sc = np.max(np.linalg.norm(P - cen, axis=1))
    seeds = []
    for rr in np.geomspace(1e-3, 1e5, 48) * sc:
        for th in np.linspace(0, 2 * math.pi, 40, endpoint=False):
            seeds.append(cen + rr * np.array([math.cos(th), math.sin(th)]))
    for v in range(3):
        for rr in np.geomspace(1e-7, 0.5, 14) * sc:
            for th in np.linspace(0, 2 * math.pi, 16, endpoint=False):
                seeds.append(P[v] + rr * np.array([math.cos(th), math.sin(th)]))
    out = []
    for s0 in seeds:
        sol = root(G, s0, method='hybr', options={'xtol': 1e-15})
        x = sol.x
        if not np.all(np.isfinite(x)) or np.linalg.norm(G(x)) > 1e-12:
            continue
        if np.min(np.linalg.norm(x[None] - P, axis=1)) < 1e-12 * sc:
            continue
        if all(np.linalg.norm(x - y) > 1e-8 * max(1.0, np.linalg.norm(y)) for y in out):
            out.append(x)
    return out


def hess(P, q, x):
    H = np.zeros((2, 2))
    for p, qi in zip(P, q):
        d = x - p
        r = np.linalg.norm(d)
        u = d / r
        H += -qi * r ** -3 * (np.eye(2) - 3 * np.outer(u, u))
    return H


hist, unmatched, tot_hist, idx_bad = {}, 0, {}, 0
done = 0
while done < n:
    zs = np.sort(rng.uniform(-1, 1, 3))
    if np.min(np.diff(zs)) < 0.05:
        continue
    signs = [1, 1, -1]
    rng.shuffle(signs)
    q = np.array(signs) * np.exp(rng.normal(size=3) * 0.7)
    done += 1
    ax = axis_equilibria(zs, q)
    hist[len(ax)] = hist.get(len(ax), 0) + 1
    # perturb the charge with index 2 (in sorted order) off the axis, inside a plane through the axis
    P = np.c_[zs, np.zeros(3)]
    P[2, 1] = 1e-3
    pts = crit_plane(P, q)
    tot_hist[len(pts)] = tot_hist.get(len(pts), 0) + 1
    ext = sum(1 for x in pts if np.linalg.det(hess(P, q, x)) > 0)
    sad = len(pts) - ext
    want = 2 if abs(q.sum()) > 1e-12 else 1
    idx_bad += (sad - ext) != want
    for z in ax:
        f2 = sum(2 * qi / abs(z - zi) ** 3 for zi, qi in zip(zs, q))
        h = sum(-qi / abs(z - zi) ** 3 for zi, qi in zip(zs, q))
        if abs(f2) > 1e-12 and abs(h) > 1e-12:
            tol = 0.05 * max(1.0, abs(z))
            if not any(np.linalg.norm(x - np.array([z, 0.0])) < tol for x in pts):
                unmatched += 1
print('(i) collinear mixed charges, %d configurations: number of on-axis equilibria %s' % (n, dict(sorted(hist.items()))))
print('(ii) non-degenerate on-axis equilibria without a nearby critical point after moving one charge off the axis by '
      '1e-3: %d' % unmatched)
print('(iii) total critical points of the perturbed (non-collinear) configurations: %s; index-count failures %d'
      % (dict(sorted(tot_hist.items())), idx_bad))
