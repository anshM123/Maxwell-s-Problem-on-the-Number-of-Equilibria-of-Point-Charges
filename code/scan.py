"""Numerical scan for Maxwell's problem with three positive charges (see README.md).

V(x) = sum_i q_i |x - p_i|^(-2 alpha) on the plane of the charges (alpha = 1/2: Coulomb in R^3).  Critical points are
the zeros of G(x) = sum_i q_i r_i^(-2 alpha - 2) (x - p_i) = -grad V / (2 alpha), found by damped Newton from a dense
barycentric grid of seeds in the triangle, deduplicated and classified by the Jacobian J = dG/dx = -Hess V / (2 alpha):
J negative definite <=> planar local minimum of V, det J < 0 <=> saddle.  Morse/Poincare-Hopf: #saddles - #minima = 2,
which every sample is checked against (a failure means a missed or spurious point).  A numerical exploration only.

usage: python scan.py N_CONFIGS [alpha] [seed] [logratio]"""
import math
import sys
import time

import numpy as np


def triangle_from_angles(A, B):
    """p1 = (0,0), p2 = (1,0), p3 above the x-axis with angles A at p1 and B at p2"""
    C = math.pi - A - B
    L13 = math.sin(B) / math.sin(C)          # law of sines with |p1 p2| = 1
    return np.array([[0.0, 0.0], [1.0, 0.0], [L13 * math.cos(A), L13 * math.sin(A)]])


def seeds_in(P, N=28, nrad=14, nang=9):
    """a barycentric grid, plus seeds clustered geometrically toward each vertex (an equilibrium next to a small charge
    sits at distance about sqrt(q_small / q_large) from it)"""
    pts = []
    for i in range(1, N):
        for j in range(1, N - i):
            k = N - i - j
            if k >= 1:
                pts.append((i * P[0] + j * P[1] + k * P[2]) / N)
    for v in range(3):
        a, b = P[(v + 1) % 3], P[(v + 2) % 3]
        for t in np.geomspace(1e-4, 0.4, nrad):
            for u in np.linspace(0.05, 0.95, nang):
                pts.append(P[v] + t * ((1 - u) * a + u * b - P[v]))
    return np.array(pts)


def field(X, P, q, alpha):
    """G(X) and J(X) for points X (S,2)"""
    G = np.zeros_like(X)
    J = np.zeros((X.shape[0], 2, 2))
    for p, qi in zip(P, q):
        d = X - p
        r2 = np.einsum('ij,ij->i', d, d)
        w = qi * r2 ** (-alpha - 1.0)
        G += w[:, None] * d
        w2 = -2.0 * (alpha + 1.0) * qi * r2 ** (-alpha - 2.0)
        J += w[:, None, None] * np.eye(2)[None] + w2[:, None, None] * d[:, :, None] * d[:, None, :]
    return G, J


def critical_points(P, q, alpha, N=28, iters=120):
    X = seeds_in(P, N)
    alive = np.ones(len(X), bool)
    for _ in range(iters):
        G, J = field(X, P, q, alpha)
        det = J[:, 0, 0] * J[:, 1, 1] - J[:, 0, 1] * J[:, 1, 0]
        ok = np.abs(det) > 1e-300
        step = np.zeros_like(X)
        step[ok, 0] = (J[ok, 1, 1] * G[ok, 0] - J[ok, 0, 1] * G[ok, 1]) / det[ok]
        step[ok, 1] = (-J[ok, 1, 0] * G[ok, 0] + J[ok, 0, 0] * G[ok, 1]) / det[ok]
        dmin = np.min([np.hypot(*(X - p).T) for p in P], axis=0)
        sl = np.hypot(step[:, 0], step[:, 1])
        fac = np.minimum(1.0, 0.5 * dmin / np.maximum(sl, 1e-300))
        X = X - step * fac[:, None]
        alive &= ok & np.isfinite(X).all(axis=1)
    G, J = field(X, P, q, alpha)
    scale = np.zeros(len(X))
    for p, qi in zip(P, q):
        d = X - p
        r2 = np.einsum('ij,ij->i', d, d)
        scale += qi * r2 ** (-alpha - 0.5)
    res = np.hypot(G[:, 0], G[:, 1]) / np.maximum(scale, 1e-300)
    good = alive & (res < 1e-10)
    # inside the open triangle (barycentric coordinates)
    T = np.array([P[1] - P[0], P[2] - P[0]]).T
    lam = np.linalg.solve(T, (X - P[0]).T).T
    inside = (lam[:, 0] > 0) & (lam[:, 1] > 0) & (lam.sum(axis=1) < 1)
    pts = []
    for x, j in zip(X[good & inside], J[good & inside]):
        if all(np.hypot(*(x - y)) > 1e-8 for y, _ in pts):
            ev = np.linalg.eigvalsh(j)
            kind = 'min' if ev[1] < 0 else ('saddle' if ev[0] < 0 < ev[1] else ('max' if ev[0] > 0 else 'degenerate'))
            pts.append((x, kind))
    return pts


def classify(P, q, alpha):
    pts = critical_points(P, q, alpha)
    m = sum(k == 'min' for _, k in pts)
    s = sum(k == 'saddle' for _, k in pts)
    return m, s, pts


if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    alpha = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
    rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 0)
    L = float(sys.argv[4]) if len(sys.argv) > 4 else math.log(1000.0)
    # sanity tests
    Peq = triangle_from_angles(math.pi / 3, math.pi / 3)
    m, s, _ = classify(Peq, np.ones(3), alpha)
    print('equilateral, equal charges: minima %d saddles %d (expect 1, 3)' % (m, s), flush=True)
    hist = {}
    bad = 0
    t0 = time.time()
    for k in range(n):
        w = rng.dirichlet(np.ones(3))                       # angles uniform on the simplex
        A, B = math.pi * w[0], math.pi * w[1]
        if min(A, B, math.pi - A - B) < 1e-3:
            continue
        P = triangle_from_angles(A, B)
        q = np.exp(rng.uniform(-L, L, 3))
        q /= q.sum()
        m, s, pts = classify(P, q, alpha)
        hist[(m, s)] = hist.get((m, s), 0) + 1
        if s - m != 2:
            bad += 1
            if bad <= 5:
                print('parity check failed: A=%.6f B=%.6f q=%s minima %d saddles %d' % (A, B, q.tolist(), m, s), flush=True)
        if m >= 2:
            print('*** m0 >= 2: A=%.10f B=%.10f q=%s minima %d saddles %d points %s'
                  % (A, B, q.tolist(), m, s, [(x.tolist(), kk) for x, kk in pts]), flush=True)
    print('alpha %.4g: %d configs in %.0fs; (minima, saddles) histogram %s; parity failures %d'
          % (alpha, n, time.time() - t0, dict(sorted(hist.items())), bad), flush=True)
