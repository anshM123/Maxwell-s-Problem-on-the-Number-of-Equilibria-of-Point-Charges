"""AUDIT of the rewritten PROOF.md: (C) the collinear case, and (M) an adversarial test of Section 7 Step 2.

 (C) collinear charges in R^3 (random positions on a line, random q, alpha): all zeros of V' on the line are located by
     sign changes on a fine grid of each gap plus bisection (V is strictly convex on each gap, so one zero per gap is
     predicted), and a 3D root search from off-line seeds checks that no critical point lies off the line.
 (M) for random shapes and alpha, the quantity  m(theta, theta') = sum_i (f_i(theta_i) - f_i(theta'_i)) (theta_i - theta'_i)
     / |theta - theta'|^2, f_i(t) = log sin t + b log sin(t - A_i), must be > 0 for all theta != theta' in the closed
     sublevel set {E <= (b+1)^2/b} of H_A (this is exactly what Step 2 needs).  Pairs are drawn (i) at random inside,
     (ii) both on the boundary {E = (b+1)^2/b}, (iii) nearly coincident on the boundary; then Nelder-Mead minimises m
     over boundary pairs from the worst samples.  Reported: min m (normalised by the median).
usage: python a9_collinear_step2.py"""
import math

import numpy as np
from scipy.optimize import minimize, root

rng = np.random.default_rng(777)

# ------------------------------------------------------------------ (C)
counts, offline = {}, 0
for k in range(300):
    alpha = [0.05, 0.25, 0.5, 1.0, 2.5, 10.0][k % 6]
    s = np.sort(rng.uniform(-1, 1, 3))
    if min(np.diff(s)) < 0.02:
        continue
    d = rng.normal(size=3); d /= np.linalg.norm(d)
    o = rng.normal(size=3)
    P = o[None] + s[:, None] * d[None]
    q = np.exp(rng.normal(size=3))
    dV = lambda t: sum(-2 * alpha * qi * np.sign(t - si) * np.abs(t - si) ** (-2 * alpha - 1) for qi, si in zip(q, s))
    zeros = 0
    for a, b in ((s[0], s[1]), (s[1], s[2])):
        t = np.linspace(a, b, 20001)[1:-1]
        v = dV(t)
        zeros += int(np.sum(np.sign(v[1:]) != np.sign(v[:-1])))
    # outside the hull: V' has one sign (no zero) -- check on a grid
    for t in (np.linspace(s[0] - 5, s[0], 2001)[:-1], np.linspace(s[2], s[2] + 5, 2001)[1:]):
        v = dV(t)
        zeros += int(np.sum(np.sign(v[1:]) != np.sign(v[:-1])))
    counts[zeros] = counts.get(zeros, 0) + 1
    beta = 2 * alpha + 2

    def G(x):
        r = np.linalg.norm(x[None] - P, axis=1)
        w = q * r ** (-beta)
        return (w[:, None] * (x[None] - P)).sum(0) / (w * r).sum()

    for t in range(30):
        x0 = o + rng.uniform(-1.2, 1.2) * d + rng.normal(size=3) * 0.5
        sol = root(G, x0, method='hybr', options={'xtol': 1e-14})
        if np.all(np.isfinite(sol.x)) and np.linalg.norm(G(sol.x)) < 1e-12:
            off = sol.x - o - np.dot(sol.x - o, d) * d
            if np.linalg.norm(off) > 1e-8:
                offline += 1
print('(C) collinear: number of zeros of V\' on the line (both gaps + outside the hull) per configuration: %s;  '
      'critical points found off the line by 3D root finding: %d' % (dict(sorted(counts.items())), offline))


# ------------------------------------------------------------------ (M)
def setup(A, b):
    k = (b + 1) ** 2 / b

    def E(th):
        de = th - A
        W = np.sin(A) / (np.sin(de) * np.sin(th))
        return W[..., 0] * W[..., 1] + W[..., 1] * W[..., 2] + W[..., 2] * W[..., 0]

    def F(th):
        return np.log(np.sin(th)) + b * np.log(np.sin(th - A))

    return k, E, F


def to3(t12):
    t12 = np.atleast_2d(t12)
    return np.c_[t12, 2 * math.pi - t12.sum(1)]


worst_all, nshape, npairs = np.inf, 0, 0
for k in range(150):
    while True:
        Ad = rng.dirichlet(np.ones(3)) * 180
        if Ad.min() > [0.5, 5, 20][k % 3]:
            break
    A = np.radians(Ad)
    alpha = [0.05, 0.25, 0.5, 1.0, 2.5, 10.0][k % 6]
    b = 2 * alpha + 1
    kk, E, F = setup(A, b)
    th0 = (math.pi + A) / 2                       # incentre angles: E = 4 < k
    nshape += 1

    def bpoint(phi):                               # boundary point of {E <= k} on the ray from th0 in direction phi
        e = np.array([math.cos(phi), math.sin(phi)])
        lo, hi = 0.0, 1e-9
        while True:
            t = to3(th0[:2] + hi * e)[0]
            if np.any(t <= A) or np.any(t >= math.pi) or E(t) > kk:
                break
            lo, hi = hi, hi * 2
        for _ in range(80):
            m = 0.5 * (lo + hi)
            t = to3(th0[:2] + m * e)[0]
            if np.all(t > A) and np.all(t < math.pi) and E(t) <= kk:
                lo = m
            else:
                hi = m
        return to3(th0[:2] + lo * e)[0]

    def mono(t1, t2):
        v = t1 - t2
        return np.sum((F(t1) - F(t2)) * v) / np.sum(v * v)

    vals = []
    phis = rng.uniform(0, 2 * math.pi, 300)
    B = np.array([bpoint(p) for p in phis])
    for i in range(300):                            # (ii) boundary-boundary pairs
        j = rng.integers(300)
        if j != i:
            vals.append((mono(B[i], B[j]), phis[i], phis[j]))
    for i in range(100):                            # (iii) nearly coincident boundary pairs
        p = phis[i]
        vals.append((mono(bpoint(p), bpoint(p + 1e-4)), p, p + 1e-4))
    for i in range(100):                            # (i) interior pairs
        s1, s2 = rng.uniform(0, 1, 2)
        vals.append((mono(th0 + s1 * (B[i] - th0), th0 + s2 * (B[(i + 7) % 300] - th0)), None, None))
    npairs += len(vals)
    med = np.median([v[0] for v in vals])
    best = min(vals, key=lambda v: v[0])
    worst = best[0] / med
    # local minimisation over boundary pairs from the worst boundary samples
    bb = sorted([v for v in vals if v[1] is not None], key=lambda v: v[0])[:3]
    for v0, p1, p2 in bb:
        res = minimize(lambda z: mono(bpoint(z[0]), bpoint(z[1])) if abs(z[0] - z[1]) > 1e-7 else 1e9,
                       [p1, p2], method='Nelder-Mead', options={'xatol': 1e-9, 'fatol': 1e-12, 'maxiter': 400})
        worst = min(worst, res.fun / med)
    worst_all = min(worst_all, worst)
print('(M) %d shapes (min angle 0.5/5/20 deg), alpha in {0.05,...,10}, %d sampled pairs + local minimisation over '
      'boundary pairs: min of the Step-2 quantity normalised by its median = %.4g (must be > 0)' % (nshape, npairs, worst_all))
