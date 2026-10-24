"""Exploration: three charges of MIXED signs (q1, q2 > 0 > q3), non-collinear.

All equilibria lie in the plane of the charges (off the plane, the weights w_i = q_i r_i^-beta would have to satisfy
sum w_i = 0 and sum w_i p_i = 0, impossible for non-collinear p_i).  In the plane, x is critical iff w is proportional
to the (signed) barycentric coordinates of x, so x lies in R+ = {l1, l2 > 0 > l3} (beyond side p1p2, inside the angle
at p3) or in R- = {l1, l2 < 0 < l3} (the vertical angle at p3).  Type: Hess V = 2 alpha kappa (beta S - I) with S =
sum l_i u_i u_i^T (trace 1), so x is an extremum iff D = det S > c (a minimum if kappa > 0, i.e. in R+, a maximum in R-).
If q1 + q2 + q3 != 0, Poincare-Hopf gives #saddles = #extrema + 2.
This script finds all planar critical points by damped Newton from many seeds and reports the (extrema, saddles)
histogram.  usage: python mixed_scan.py N [alpha] [seed]"""
import math
import sys

import numpy as np

N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
alpha = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
beta = 2 * alpha + 2


def grad_hess(X, P, q):
    G = np.zeros_like(X)
    H = np.zeros((len(X), 2, 2))
    for p, qi in zip(P, q):
        d = X - p
        r2 = np.einsum('ij,ij->i', d, d)
        w = qi * r2 ** (-beta / 2)
        G += w[:, None] * d
        H += w[:, None, None] * (np.eye(2)[None] - beta * d[:, :, None] * d[:, None, :] / r2[:, None, None])
    return -2 * alpha * G, -2 * alpha * H


def critical_points(P, q, iters=200):
    c0 = P.mean(axis=0)
    diam = max(np.linalg.norm(P[i] - P[j]) for i in range(3) for j in range(3))
    seeds = []
    for rad in diam * np.geomspace(0.02, 12, 40):
        for ang in np.linspace(0, 2 * math.pi, 48, endpoint=False):
            seeds.append(c0 + rad * np.array([math.cos(ang), math.sin(ang)]))
    for p in P:
        for rad in diam * np.geomspace(1e-3, 0.5, 25):
            for ang in np.linspace(0, 2 * math.pi, 24, endpoint=False):
                seeds.append(p + rad * np.array([math.cos(ang), math.sin(ang)]))
    X = np.array(seeds)
    for _ in range(iters):
        G, H = grad_hess(X, P, q)
        det = H[:, 0, 0] * H[:, 1, 1] - H[:, 0, 1] * H[:, 1, 0]
        ok = np.abs(det) > 1e-300
        st = np.zeros_like(X)
        st[ok, 0] = (H[ok, 1, 1] * G[ok, 0] - H[ok, 0, 1] * G[ok, 1]) / det[ok]
        st[ok, 1] = (-H[ok, 1, 0] * G[ok, 0] + H[ok, 0, 0] * G[ok, 1]) / det[ok]
        dmin = np.min([np.linalg.norm(X - p, axis=1) for p in P], axis=0)
        sl = np.linalg.norm(st, axis=1)
        fac = np.minimum(1.0, 0.5 * dmin / np.maximum(sl, 1e-300))
        X = X - st * fac[:, None]
        X[~np.isfinite(X).all(axis=1)] = c0
    G, H = grad_hess(X, P, q)
    scale = sum(abs(qi) * np.linalg.norm(X - p, axis=1) ** (-beta + 1) for p, qi in zip(P, q))
    good = (np.linalg.norm(G, axis=1) / scale < 1e-10) & (np.linalg.norm(X - c0, axis=1) < 1e3 * diam)
    pts = []
    for x, h in zip(X[good], H[good]):
        if all(np.linalg.norm(x - y) > 1e-7 * diam for y, _ in pts):
            ev = np.linalg.eigvalsh(h)
            kind = 'min' if ev[0] > 0 else ('max' if ev[1] < 0 else ('saddle' if ev[0] < 0 < ev[1] else 'deg'))
            pts.append((x, kind))
    return pts


if __name__ == '__main__':
    rng = np.random.default_rng(seed)
    hist, bad = {}, 0
    worst = None
    for it in range(N):
        while True:
            A = rng.dirichlet(np.ones(3)) * math.pi
            if A.min() > math.radians(3):
                break
        s = math.sin(A[1]) / math.sin(A[2])
        P = np.array([[0.0, 0.0], [1.0, 0.0], [s * math.cos(A[0]), s * math.sin(A[0])]])
        q = np.exp(rng.normal(scale=1.5, size=3))
        q[2] = -q[2]
        pts = critical_points(P, q)
        ne = sum(k in ('min', 'max') for _, k in pts)
        ns = sum(k == 'saddle' for _, k in pts)
        key = (sum(k == 'min' for _, k in pts), sum(k == 'max' for _, k in pts), ns)
        hist[key] = hist.get(key, 0) + 1
        if ns != ne + 2:
            bad += 1
        if ne >= 2:
            worst = (A, q, pts)
    print('alpha %.2f, %d mixed-sign configurations: (minima, maxima, saddles) histogram %s; parity failures %d'
          % (alpha, N, dict(sorted(hist.items())), bad))
    if worst is not None:
        A, q, pts = worst
        print('TWO OR MORE EXTREMA at angles %s q %s: %s' % (np.degrees(A).round(2).tolist(), q.round(4).tolist(),
                                                           [(x.round(4).tolist(), k) for x, k in pts]))
