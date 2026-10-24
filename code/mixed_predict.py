"""Mixed signs, targeted test.  Negative charge at p3, positive at p1, p2 (any sign pattern reduces to this by relabelling
and V -> -V).  Pick x beyond side p1p2, outside the circumcircle, with D(x) > c (the region where minima live), set
q = (lambda_i r_i^beta) (signed barycentrics: (+, +, -)), and count ALL planar critical points by a global Newton search.
Prediction if Maxwell's bound holds: exactly one minimum (at x), no maximum, and 3 saddles (2 if q1 + q2 + q3 = 0).
Points are sampled uniformly in the minimum region and also near its boundary (D close to c).
usage: python mixed_predict.py [n] [alpha] [seed] [near_fraction]"""
import math
import sys

import numpy as np

from mixed_scan import critical_points

n = int(sys.argv[1]) if len(sys.argv) > 1 else 200
alpha = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
near = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
beta = 2 * alpha + 2
b = beta - 1
c = b / (b + 1) ** 2


def bary(P, x):
    T = np.array([P[:, 0], P[:, 1], np.ones(3)])
    return np.linalg.solve(T, np.array([x[0], x[1], 1.0]))


def D_of(P, x):
    lam = bary(P, x)
    S = np.zeros((2, 2))
    for i in range(3):
        u = (x - P[i]) / np.linalg.norm(x - P[i])
        S += lam[i] * np.outer(u, u)
    return np.linalg.det(S)


if __name__ == '__main__':
    rng = np.random.default_rng(seed)
    hist, bad, done = {}, 0, 0
    while done < n:
        A = rng.dirichlet(np.ones(3)) * math.pi
        if A.min() < math.radians(4):
            continue
        s = math.sin(A[1]) / math.sin(A[2])
        P = np.array([[0.0, 0.0], [1.0, 0.0], [s * math.cos(A[0]), s * math.sin(A[0])]])
        a = np.array([np.linalg.norm(P[(i + 1) % 3] - P[(i + 2) % 3]) for i in range(3)])
        I3 = (a[0] * P[0] + a[1] * P[1] - a[2] * P[2]) / (a[0] + a[1] - a[2])     # excentre opposite p3
        want_near = rng.uniform() < near
        x = None
        for _ in range(4000):
            y = I3 + rng.normal(size=2) * rng.choice([0.05, 0.2, 0.6]) * np.linalg.norm(I3 - P[2])
            lam = bary(P, y)
            if not (lam[0] > 0 and lam[1] > 0 and lam[2] < 0):
                continue
            D = D_of(P, y)
            if D > c and (not want_near or D < c * 1.02):
                x = y
                break
        if x is None:
            continue
        done += 1
        lam, r = bary(P, x), np.linalg.norm(x - P, axis=1)
        q = lam * r ** beta
        q /= np.abs(q).max()
        pts = critical_points(P, q)
        m = sum(k == 'min' for _, k in pts)
        M = sum(k == 'max' for _, k in pts)
        sd = sum(k == 'saddle' for _, k in pts)
        at_x = [k for y, k in pts if np.linalg.norm(y - x) < 1e-6 * np.linalg.norm(P[1] - P[0])]
        hist[(m, M, sd)] = hist.get((m, M, sd), 0) + 1
        if not (m == 1 and M == 0 and sd == 3 and at_x == ['min']):
            bad += 1
            print('MISMATCH angles %s q %s: (min, max, saddle) = (%d, %d, %d), at x: %s'
                  % (np.degrees(A).round(2).tolist(), q.round(5).tolist(), m, M, sd, at_x), flush=True)
    print('alpha %.2f: %d configurations with a forced minimum; histogram %s; mismatches %d'
          % (alpha, n, dict(sorted(hist.items())), bad))
