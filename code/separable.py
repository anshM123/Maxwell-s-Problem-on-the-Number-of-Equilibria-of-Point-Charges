"""The reduction is separable in angle coordinates.

For x in the open triangle, with theta_i the angle at x subtending the side opposite p_i, A_i the triangle's angles
and a_i its sides, the charges that make x critical are
    q_i  ∝  lambda_i r_i^beta  ∝  sin(theta_i) * sin(theta_i - A_i)^(beta-1) / a_i^(beta-1).
(theta_i - A_i are the angles of the pedal triangle of x; r_i ∝ sin(theta_i - A_i)/a_i by inversion at x.)
So Psi = grad_theta H with H = sum_i F_i(theta_i) (+ linear), F_i' = f_i = log sin t + (beta-1) log sin(t - A_i), and
Hess_theta H = diag(d) restricted to {sum v = 0}, with d_i = f_i'(theta_i).
Checks:
  1. the charge formula against lambda_i r_i^beta at random points;
  2. sign(D - c) = sign(e2(d)), e2 = d1 d2 + d2 d3 + d3 d1;
  3. N = {e2 > 0, sum d < 0} (negative-definite region) does not meet the hexagon of angle coordinates;
  4. P = {e2 > 0, sum d > 0} is star-shaped from the centroid of T0 = {theta_i < t_i*}, and its boundary curvature,
     scaled by the size of P, is positive (P convex).
usage: python separable.py [alpha] [n_shapes] [min_angle_deg]"""
import math
import sys

import numpy as np

alpha = float(sys.argv[1]) if __name__ == '__main__' and len(sys.argv) > 1 else 0.5
nsh = int(sys.argv[2]) if __name__ == '__main__' and len(sys.argv) > 2 else 200
amin = math.radians(float(sys.argv[3])) if __name__ == '__main__' and len(sys.argv) > 3 else math.radians(1.0)
beta = 2 * alpha + 2
b = beta - 1
c = (beta - 1) / beta ** 2


def tri(A):
    s = math.sin(A[1]) / math.sin(A[2])          # |p1 p3| when |p1 p2| = 1 (law of sines)
    return np.array([[0.0, 0.0], [1.0, 0.0], [s * math.cos(A[0]), s * math.sin(A[0])]])


def bary(P, x):
    T = np.array([P[:, 0], P[:, 1], np.ones(3)])
    return np.linalg.solve(T, np.array([x[0], x[1], 1.0]))


def angles_at(P, x):
    th = np.empty(3)
    for i in range(3):
        u, v = P[(i + 1) % 3] - x, P[(i + 2) % 3] - x
        th[i] = math.atan2(abs(u[0] * v[1] - u[1] * v[0]), u @ v)
    return th


def D_of(P, x):
    lam = bary(P, x)
    S = np.zeros((2, 2))
    for i in range(3):
        u = (x - P[i]) / np.linalg.norm(x - P[i])
        S += lam[i] * np.outer(u, u)
    return np.linalg.det(S)


def g1(th, A):
    return b / np.tan(th - A) + 1 / np.tan(th)


def g2(th, A):
    return -b / np.sin(th - A) ** 2 - 1 / np.sin(th) ** 2


def g3(th, A):
    return 2 * b * np.cos(th - A) / np.sin(th - A) ** 3 + 2 * np.cos(th) / np.sin(th) ** 3


def e2(d):
    return d[0] * d[1] + d[1] * d[2] + d[2] * d[0]


def tstar(A):
    k = 1 / np.tan(A)
    cd = (np.sqrt((b + 1) ** 2 * k ** 2 + 4 * b) - (b + 1) * k) / (2 * b)
    return A + np.arctan2(1.0, cd)


def F_of(th, A):
    return e2(g1(th, A))


def boundary_curvature(th, A):
    d = g1(th, A)
    e = d.sum() - d
    s2, s3 = g2(th, A), g3(th, A)
    G = e * s2
    Gp = G - G.mean()
    v = np.array([G[2] - G[1], G[0] - G[2], G[1] - G[0]])
    Q = np.sum(e * s3 * v ** 2) + 2 * (s2[0] * s2[1] * v[0] * v[1] + s2[1] * s2[2] * v[1] * v[2] + s2[0] * s2[2] * v[0] * v[2])
    return -Q / (v @ v * np.linalg.norm(Gp)), np.linalg.norm(Gp)


def smax(th0, v, A):
    """largest s with A_i < th0 + s v < pi"""
    best = np.inf
    for i in range(3):
        if v[i] > 0:
            best = min(best, (math.pi - th0[i]) / v[i])
        elif v[i] < 0:
            best = min(best, (A[i] - th0[i]) / v[i])
    return best


if __name__ == '__main__':
    rng = np.random.default_rng(1)
    err1, bad2, n12 = 0.0, 0, 0
    while n12 < 3000:
        A = rng.dirichlet(np.ones(3)) * math.pi
        if A.min() < amin:
            continue
        n12 += 1
        P = tri(A)
        a = np.array([np.linalg.norm(P[(i + 1) % 3] - P[(i + 2) % 3]) for i in range(3)])
        x = rng.dirichlet(np.ones(3)) @ P
        lam, r = bary(P, x), np.linalg.norm(x - P, axis=1)
        q1 = lam * r ** beta
        th = angles_at(P, x)
        q2 = np.sin(th) * np.sin(th - A) ** b / a ** b
        err1 = max(err1, np.max(np.abs((q1 / q1.sum()) / (q2 / q2.sum()) - 1)))
        Dm, E = D_of(P, x) - c, F_of(th, A)
        if Dm * E < 0 and abs(Dm) > 1e-9:
            bad2 += 1
    print('1. charge formula: max relative error %.2e over %d points' % (err1, n12))
    print('2. sign(D - c) != sign(e2(d)) at %d of %d points' % (bad2, n12))

    worstN, worst_kR, worst_star, sum_tstar_min = -np.inf, np.inf, 0, np.inf
    worst_shape = None
    for k in range(nsh):
        while True:
            A = rng.dirichlet(np.ones(3) * 0.6) * math.pi
            if A.min() >= amin:
                break
        ts = tstar(A)
        sum_tstar_min = min(sum_tstar_min, ts.sum() - 2 * math.pi)
        # 3. grid over the hexagon in (theta_1, theta_2)
        n = 400
        t1 = np.linspace(A[0], math.pi, n + 2)[1:-1]
        t2 = np.linspace(A[1], math.pi, n + 2)[1:-1]
        T1, T2 = np.meshgrid(t1, t2)
        T3 = 2 * math.pi - T1 - T2
        ok = (T3 > A[2]) & (T3 < math.pi)
        th = np.array([T1[ok], T2[ok], T3[ok]])
        d = np.array([g1(th[i], A[i]) for i in range(3)])
        E, S = e2(d), d.sum(axis=0)
        neg = E > 0
        if np.any(neg & (S < 0)):
            worstN = max(worstN, np.max(-S[neg & (S < 0)]))
        # 4. star-shapedness and boundary curvature of P from the centroid of T0
        V = np.array([[ts[0], ts[1], 2 * math.pi - ts[0] - ts[1]],
                      [ts[0], 2 * math.pi - ts[0] - ts[2], ts[2]],
                      [2 * math.pi - ts[1] - ts[2], ts[1], ts[2]]])
        th0 = V.mean(axis=0)
        assert F_of(th0, A) > 0
        e1, e2v = np.array([1, -1, 0]) / math.sqrt(2), np.array([1, 1, -2]) / math.sqrt(6)
        rad, kap = [], []
        for phi in np.linspace(0, 2 * math.pi, 720, endpoint=False):
            v = math.cos(phi) * e1 + math.sin(phi) * e2v
            sm = smax(th0, v, A)
            ss = np.linspace(0, sm, 4001)[1:-1]
            pts = th0[:, None] + ss[None, :] * v[:, None]
            Fv = F_of(pts, A[:, None])
            sgn = Fv > 0
            first = np.argmin(sgn)                # first index with F <= 0
            if sgn[first] or np.any(sgn[first:]):
                worst_star += 1                   # not star-shaped along this ray (or no exit)
                continue
            lo, hi = ss[first - 1] if first > 0 else 0.0, ss[first]
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if F_of(th0 + mid * v, A) > 0:
                    lo = mid
                else:
                    hi = mid
            tb = th0 + lo * v
            kp, gn = boundary_curvature(tb, A)
            rad.append(lo)
            kap.append(kp)
        R = max(rad)
        kR = min(kap) * R
        if kR < worst_kR:
            worst_kR, worst_shape = kR, A.copy()
    print('3. Sigma t* - 2 pi >= %.4f;  negative-definite region met: %s'
          % (sum_tstar_min, 'NO' if worstN == -np.inf else 'YES (max -sum d = %.3g)' % worstN))
    print('4. rays not star-shaped: %d;  min over shapes of (boundary curvature x radius) = %.4f at angles %s deg'
          % (worst_star, worst_kR, np.round(np.degrees(worst_shape), 2).tolist()))
