"""Numerical check of the identities behind the proof (PROOF.md, Lemmas A-D), at random triangles and points.

With theta_i the angle at x subtending the side opposite p_i, A_i the angles of the triangle, delta_i = theta_i - A_i,
b = beta - 1, P_x = R^2 - |x - O|^2 (R, O circumradius and circumcentre), W_i = cot(delta_i) - cot(theta_i):
  (A) r_j r_k sin(delta_i) = P_x sin(A_i)                       (Euler's pedal-triangle theorem)
  (A') lambda_i r_i^beta  ∝  sin(theta_i) sin(delta_i)^b / sin(A_i)^b
  (B) W_i = sin(A_i) / (sin(delta_i) sin(theta_i)),   D(x) = 1 / e2(W)
  (C) d_i = cot(theta_i) + b cot(delta_i):   e2(d) = (b + 1)^2 - b e2(W)
  (D) e2(W) is convex in theta (finite-difference Hessian on the plane sum theta = 2 pi)
  (E) the incentre has theta_i = (pi + A_i)/2, D = 1/4 and all d_i > 0
usage: python key_identities.py [n_points] [alpha]"""
import math
import sys

import numpy as np

from separable import tri, bary, angles_at, D_of

n = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
alpha = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
beta = 2 * alpha + 2
b = beta - 1


def e2(v):
    return v[0] * v[1] + v[1] * v[2] + v[2] * v[0]


def circum(P):
    (ax, ay), (bx, by), (cx, cy) = P
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
    uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
    O = np.array([ux, uy])
    return O, np.linalg.norm(P[0] - O)


def W_of(th, A):
    return 1 / np.tan(th - A) - 1 / np.tan(th)


if __name__ == '__main__':
    rng = np.random.default_rng(17)
    errA = errA2 = errB1 = errB2 = errC = 0.0
    for k in range(n):
        A = rng.dirichlet(np.ones(3) * 0.8) * math.pi
        if A.min() < 1e-3:
            continue
        P = tri(A)
        x = rng.dirichlet(np.ones(3) * 0.8) @ P
        lam, r = bary(P, x), np.linalg.norm(x - P, axis=1)
        if lam.min() < 1e-6:
            continue
        th = angles_at(P, x)
        dl = th - A
        O, R = circum(P)
        Px = R * R - np.sum((x - O) ** 2)
        for i in range(3):
            j, kk = (i + 1) % 3, (i + 2) % 3
            errA = max(errA, abs(r[j] * r[kk] * math.sin(dl[i]) - Px * math.sin(A[i])) / (Px * math.sin(A[i])))
        q1 = lam * r ** beta
        q2 = np.sin(th) * np.sin(dl) ** b / np.sin(A) ** b
        errA2 = max(errA2, np.max(np.abs((q1 / q1.sum()) / (q2 / q2.sum()) - 1)))
        W = W_of(th, A)
        errB1 = max(errB1, np.max(np.abs(W / (np.sin(A) / (np.sin(dl) * np.sin(th))) - 1)))
        errB2 = max(errB2, abs(D_of(P, x) * e2(W) - 1))
        d = 1 / np.tan(th) + b / np.tan(dl)
        errC = max(errC, abs(e2(d) - ((b + 1) ** 2 - b * e2(W))) / (1 + abs(e2(d)) + b * e2(W)))
    print('(A)  r_j r_k sin(delta_i) = P_x sin(A_i):        max rel err %.2e' % errA)
    print("(A') charge formula:                              max rel err %.2e" % errA2)
    print('(B)  W_i = sin A_i/(sin delta_i sin theta_i):      max rel err %.2e' % errB1)
    print('(B)  D * e2(W) = 1:                                max abs err %.2e' % errB2)
    print('(C)  e2(d) = (b+1)^2 - b e2(W):                    max rel err %.2e' % errC)
    # (D) convexity of e2(W) on the plane, finite differences in theta
    worst = np.inf
    E = np.array([[1, -1, 0], [1, 1, -2]], dtype=float)
    E[0] /= math.sqrt(2)
    E[1] /= math.sqrt(6)
    for k in range(3000):
        A = rng.dirichlet(np.ones(3) * 0.6) * math.pi
        if A.min() < 1e-3:
            continue
        th = A + rng.dirichlet(np.ones(3)) * math.pi
        if np.any(th >= math.pi - 1e-4) or np.any(th - A < 1e-4):
            continue
        h = 1e-4 * min(np.min(th - A), np.min(math.pi - th))
        f = lambda t: e2(W_of(t, A))
        Hm = np.zeros((2, 2))
        for a in range(2):
            for c in range(2):
                ea, ec = E[a] * h, E[c] * h
                Hm[a, c] = (f(th + ea + ec) - f(th + ea - ec) - f(th - ea + ec) + f(th - ea - ec)) / (4 * h * h)
        ev = np.linalg.eigvalsh((Hm + Hm.T) / 2)
        worst = min(worst, ev[0] / ev[1])
    print('(D)  min over samples of (smaller/larger eigenvalue) of Hess e2(W): %.3e  (>= 0: convex)' % worst)
    # (E) incentre
    errE = 0.0
    for k in range(2000):
        A = rng.dirichlet(np.ones(3)) * math.pi
        P = tri(A)
        a = np.array([np.linalg.norm(P[(i + 1) % 3] - P[(i + 2) % 3]) for i in range(3)])
        I = (a @ P) / a.sum()
        th = angles_at(P, I)
        d = 1 / np.tan(th) + b / np.tan(th - A)
        errE = max(errE, np.max(np.abs(th - (math.pi + A) / 2)), abs(D_of(P, I) - 0.25))
        assert np.all(d > 0)
    print('(E)  incentre: theta = (pi + A)/2 and D = 1/4, max err %.2e; all d_i > 0 there' % errE)
