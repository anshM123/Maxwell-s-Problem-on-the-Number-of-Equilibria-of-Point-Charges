"""Mixed signs (q3 < 0 < q1, q2): curvature of the boundary of the minimum region in angle coordinates.

In the region beyond side p1p2 and outside the circumcircle, with t_i = angle at x subtending p_j p_k (i = 1, 2),
t1 in (0, A1), t2 in (0, A2), s = t1 + t2:
    U1 = cot t1 + cot(A1 - t1),  U2 = cot t2 + cot(A2 - t2),  U3 = cot s + cot(A1 + A2 - s),
    E = U1 U2 - U3 (U1 + U2) = 1/D,   minimum <=> E < k = (b + 1)^2 / b.
E >= 4 with equality at the excentre point (A1/2, A2/2).  For each shape, trace the level curve {E = k} along rays from
that point, check star-shapedness, and compute the signed curvature (exact derivatives) scaled by the max radius.
usage: python mixed_boundary.py [alpha] [n_shapes] [min_angle_deg]"""
import math
import sys

import numpy as np

alpha = float(sys.argv[1]) if __name__ == '__main__' and len(sys.argv) > 1 else 0.5
nsh = int(sys.argv[2]) if __name__ == '__main__' and len(sys.argv) > 2 else 400
amin = math.radians(float(sys.argv[3])) if __name__ == '__main__' and len(sys.argv) > 3 else math.radians(0.5)
b = 2 * alpha + 1
k = (b + 1) ** 2 / b


def U(t, A):
    a, c = t, A - t
    u = 1 / np.tan(a) + 1 / np.tan(c)
    u1 = -1 / np.sin(a) ** 2 + 1 / np.sin(c) ** 2
    u2 = 2 / np.tan(a) / np.sin(a) ** 2 + 2 / np.tan(c) / np.sin(c) ** 2
    return u, u1, u2


def E_derivs(t1, t2, A1, A2):
    a, a1, a2 = U(t1, A1)
    c, c1, c2 = U(t2, A2)
    w, w1, w2 = U(t1 + t2, A1 + A2)
    E = a * c - w * (a + c)
    E1 = a1 * c - w1 * (a + c) - w * a1
    E2 = a * c1 - w1 * (a + c) - w * c1
    E11 = a2 * c - w2 * (a + c) - 2 * w1 * a1 - w * a2
    E22 = a * c2 - w2 * (a + c) - 2 * w1 * c1 - w * c2
    E12 = a1 * c1 - w2 * (a + c) - w1 * c1 - w1 * a1
    return E, E1, E2, E11, E12, E22


if __name__ == '__main__':
    rng = np.random.default_rng(4)
    worst, worst_shape, nonstar = np.inf, None, 0
    for it in range(nsh):
        while True:
            A = rng.dirichlet(np.ones(3) * 0.6) * math.pi
            if A.min() >= amin:
                break
        A1, A2 = A[0], A[1]
        c0 = np.array([A1 / 2, A2 / 2])
        rad, kap = [], []
        for phi in np.linspace(0, 2 * math.pi, 360, endpoint=False):
            v = np.array([math.cos(phi), math.sin(phi)])
            # largest step inside the rectangle (0, A1) x (0, A2)
            smax = min((A1 - c0[0]) / v[0] if v[0] > 0 else (-c0[0] / v[0] if v[0] < 0 else np.inf),
                       (A2 - c0[1]) / v[1] if v[1] > 0 else (-c0[1] / v[1] if v[1] < 0 else np.inf))
            ss = np.linspace(0, smax, 3001)[1:-1]
            Ev = E_derivs(c0[0] + ss * v[0], c0[1] + ss * v[1], A1, A2)[0]
            inside = Ev < k
            first = np.argmin(inside)
            if inside[first] or np.any(inside[first:]):
                nonstar += 1
                continue
            lo, hi = (ss[first - 1] if first > 0 else 0.0), ss[first]
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if E_derivs(c0[0] + mid * v[0], c0[1] + mid * v[1], A1, A2)[0] < k:
                    lo = mid
                else:
                    hi = mid
            t = c0 + lo * v
            E, E1, E2, E11, E12, E22 = E_derivs(t[0], t[1], A1, A2)
            Q = E11 * E2 ** 2 - 2 * E12 * E1 * E2 + E22 * E1 ** 2
            kap.append(Q / (E1 ** 2 + E2 ** 2) ** 1.5)      # curvature of {E = k}, > 0 = curving towards {E < k}
            rad.append(lo)
        if kap:
            val = min(kap) * max(rad)
            if val < worst:
                worst, worst_shape = val, np.degrees(A).round(3).tolist()
    print('alpha %.3f (k = %.4f): rays not star-shaped %d; min over shapes of curvature x radius = %.5f at angles %s'
          % (alpha, k, nonstar, worst, worst_shape))
