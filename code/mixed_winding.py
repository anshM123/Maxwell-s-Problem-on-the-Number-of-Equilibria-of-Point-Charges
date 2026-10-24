"""Mixed signs: is the charge map psi injective on the minimum region M' = {E < k}?

psi is a local diffeomorphism of constant orientation on M' (Hess H is definite there).  So the number of preimages in
M' of a point y equals |winding number of psi(boundary of M') around y|.  We trace the boundary {E = k} (star-shaped
from the excentre point (A1/2, A2/2)) with many rays, map it by psi, and compute the maximum |winding number| over a
grid of points y (and at the midpoints of the image polygon's self-intersections, if any).  Max winding 1 means at most
one minimum for every choice of charges at this shape.
usage: python mixed_winding.py alpha [n_shapes] [min_angle_deg] [seed]   or   python mixed_winding.py alpha A1 A2 (deg)"""
import math
import sys

import numpy as np

alpha = float(sys.argv[1])
b = 2 * alpha + 1
k = (b + 1) ** 2 / b


def make(A1, A2):
    A3 = math.pi - A1 - A2

    def U(t, A):
        return 1 / np.tan(t) + 1 / np.tan(A - t)

    def E(t1, t2):
        u1, u2, u3 = U(t1, A1), U(t2, A2), U(t1 + t2, A1 + A2)
        return u1 * u2 - u3 * (u1 + u2)

    def psi(t1, t2):
        s = t1 + t2
        p1 = np.log(np.sin(t1)) + b * np.log(np.sin(A1 - t1))
        p2 = np.log(np.sin(t2)) + b * np.log(np.sin(A2 - t2))
        p3 = np.log(np.sin(s)) + b * np.log(np.sin(s + A3))
        return np.stack([p1 - p3, p2 - p3], axis=-1)
    return E, psi


def boundary(A1, A2, E, nray=4000):
    c0 = np.array([A1 / 2, A2 / 2])
    pts = []
    for phi in np.linspace(0, 2 * math.pi, nray, endpoint=False):
        v = np.array([math.cos(phi), math.sin(phi)])
        smax = min((A1 - c0[0]) / v[0] if v[0] > 0 else (-c0[0] / v[0] if v[0] < 0 else np.inf),
                   (A2 - c0[1]) / v[1] if v[1] > 0 else (-c0[1] / v[1] if v[1] < 0 else np.inf))
        lo, hi = 0.0, smax * (1 - 1e-9)
        for _ in range(70):
            mid = 0.5 * (lo + hi)
            if E(c0[0] + mid * v[0], c0[1] + mid * v[1]) < k:
                lo = mid
            else:
                hi = mid
        pts.append(c0 + lo * v)
    return np.array(pts)


def winding(poly, y):
    d = poly - y
    ang = np.arctan2(d[:, 1], d[:, 0])
    dif = np.diff(np.concatenate([ang, ang[:1]]))
    dif = (dif + math.pi) % (2 * math.pi) - math.pi
    return int(round(dif.sum() / (2 * math.pi)))


def max_winding(A1, A2):
    E, psi = make(A1, A2)
    bd = boundary(A1, A2, E)
    img = psi(bd[:, 0], bd[:, 1])
    lo, hi = img.min(axis=0), img.max(axis=0)
    best = 0
    g = 120
    for yx in np.linspace(lo[0], hi[0], g)[1:-1]:
        for yy in np.linspace(lo[1], hi[1], g)[1:-1]:
            best = max(best, abs(winding(img, np.array([yx, yy]))))
    return best, img


if __name__ == '__main__':
    if len(sys.argv) == 4:
        A1, A2 = math.radians(float(sys.argv[2])), math.radians(float(sys.argv[3]))
        w, _ = max_winding(A1, A2)
        print('alpha %.2f shape (%.3f, %.3f, %.3f) deg: max |winding| of psi(boundary) = %d'
              % (alpha, math.degrees(A1), math.degrees(A2), 180 - math.degrees(A1) - math.degrees(A2), w))
    else:
        nsh = int(sys.argv[2]) if len(sys.argv) > 2 else 100
        amin = math.radians(float(sys.argv[3])) if len(sys.argv) > 3 else math.radians(2)
        rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 0)
        hist = {}
        bad = []
        for it in range(nsh):
            while True:
                A = rng.dirichlet(np.ones(3) * 0.7) * math.pi
                if A.min() >= amin:
                    break
            w, _ = max_winding(A[0], A[1])
            hist[w] = hist.get(w, 0) + 1
            if w > 1:
                bad.append(np.degrees(A).round(2).tolist())
        print('alpha %.2f: %d shapes, histogram of max |winding| %s; shapes with winding > 1: %s'
              % (alpha, nsh, dict(sorted(hist.items())), bad[:5]))
