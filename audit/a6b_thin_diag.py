"""AUDIT items 6/7, diagnostic for thin shapes flagged by a67_shape_scan.py (fold image self-intersections and
near-zero normalised curvature of the boundary of theta(M)).

For a shape: high-resolution boundary of theta(M) (rays from theta*), its image Gamma under Psi; every crossing of
Gamma with itself is reported with the two ray parameters, the distance between the two boundary points in theta,
and the winding numbers of Gamma at points sampled around the crossing.  The curvature of the boundary is printed
along the boundary (mpmath-free float64 central differences; step 1e-3 rho).
usage: python a6b_thin_diag.py A1_deg A2_deg [alpha] [npsi]"""
import math
import sys

import numpy as np

sys.argv += [] if len(sys.argv) > 3 else []
A1d, A2d = float(sys.argv[1]), float(sys.argv[2])
alpha = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
npsi = int(sys.argv[4]) if len(sys.argv) > 4 else 8000
import a67_shape_scan as m   # noqa: E402  (our own audit module; argv parsing there uses defaults)

m.beta, m.c = 2 * alpha + 2, (2 * alpha + 1) / (2 * alpha + 2) ** 2
S = m.Shape(math.radians(A1d), math.radians(A2d))
from scipy.optimize import root  # noqa: E402

zs = []
for s0 in m.seeds(S):
    sol = root(lambda x: [S.Z(x[None])[0].real, S.Z(x[None])[0].imag], s0, method='hybr', options={'xtol': 1e-15})
    if np.all(np.isfinite(sol.x)) and S.bary(sol.x).min() > 0 and abs(S.Z(sol.x[None])[0]) < 1e-10:
        zs.append(sol.x)
xs = zs[0]
ts = S.theta_of_x(xs[None])[0]
print('shape', A1d, A2d, 180 - A1d - A2d, ' x* =', xs, ' bary(x*) =', S.bary(xs), ' theta* =', ts)
psi = np.linspace(0, 2 * math.pi, npsi, endpoint=False)
E = np.c_[np.cos(psi), np.sin(psi)]
lo, hi = np.zeros(npsi), np.full(npsi, 1e-6)
for _ in range(80):
    grow = S.Dth(ts[None] + hi[:, None] * E) > m.c
    if not grow.any():
        break
    lo = np.where(grow, hi, lo)
    hi = np.where(grow, 2 * hi, hi)
for _ in range(70):
    mid = 0.5 * (lo + hi)
    inside = S.Dth(ts[None] + mid[:, None] * E) > m.c
    lo, hi = np.where(inside, mid, lo), np.where(inside, hi, mid)
rho = 0.5 * (lo + hi)
B = ts[None] + rho[:, None] * E
Xb, lb = S.x_of_theta(B)
Gam = S.Psi(Xb, lb)
print('rho range %.3e .. %.3e (aspect %.1f); Gamma extent %s' % (rho.min(), rho.max(), rho.max() / rho.min(),
                                                                  np.ptp(Gam, 0)))
# curvature along the boundary
hh = 1e-3 * rho[:, None]
F = lambda dt: S.Dth(B + dt)
e1, e2 = np.array([1.0, 0.0]), np.array([0.0, 1.0])
F0 = F(0 * B)
F1 = (F(hh * e1) - F(-hh * e1)) / (2 * hh[:, 0]); F2 = (F(hh * e2) - F(-hh * e2)) / (2 * hh[:, 0])
F11 = (F(hh * e1) - 2 * F0 + F(-hh * e1)) / hh[:, 0] ** 2; F22 = (F(hh * e2) - 2 * F0 + F(-hh * e2)) / hh[:, 0] ** 2
F12 = (F(hh * (e1 + e2)) - F(hh * (e1 - e2)) - F(hh * (e2 - e1)) + F(-hh * (e1 + e2))) / (4 * hh[:, 0] ** 2)
kap = -(F11 * F2 ** 2 - 2 * F12 * F1 * F2 + F22 * F1 ** 2) / (F1 ** 2 + F2 ** 2) ** 1.5
print('curvature * mean rho: min %.3e at psi=%.4f, max %.3e;  quantiles %s' % (
    kap.min() * rho.mean(), psi[kap.argmin()], kap.max() * rho.mean(),
    np.round(np.quantile(kap * rho.mean(), [0, .01, .1, .5, .9, 1]), 6).tolist()))
# self-intersections of Gamma
n = len(Gam)
a1, a2 = Gam, np.roll(Gam, -1, 0)
o = lambda p, q, r: (q[..., 0] - p[..., 0]) * (r[..., 1] - p[..., 1]) - (q[..., 1] - p[..., 1]) * (r[..., 0] - p[..., 0])
hits = []
for i in range(n):
    j = np.arange(i + 2, n) if i > 0 else np.arange(2, n - 1)
    b1, b2 = a1[j], a2[j]
    s1 = o(a1[i], a2[i], b1) * o(a1[i], a2[i], b2)
    s2 = o(b1, b2, a1[i][None]) * o(b1, b2, a2[i][None])
    for jj in j[(s1 < 0) & (s2 < 0)]:
        hits.append((i, int(jj)))
print('self-intersections of Gamma: %d' % len(hits))
for i, j in hits[:10]:
    # intersection point
    p, r = a1[i], a2[i] - a1[i]
    q_, s = a1[j], a2[j] - a1[j]
    t = ((q_ - p)[0] * s[1] - (q_ - p)[1] * s[0]) / (r[0] * s[1] - r[1] * s[0])
    X0 = p + t * r
    # winding numbers at points around the crossing
    scale = max(np.linalg.norm(r), np.linalg.norm(s))
    ws = []
    for ang in np.linspace(0, 2 * math.pi, 8, endpoint=False):
        w = X0 + 0.5 * scale * np.array([math.cos(ang), math.sin(ang)])
        dv = Gam - w
        a = np.arctan2(dv[:, 1], dv[:, 0])
        da = np.diff(np.r_[a, a[0]])
        da = (da + math.pi) % (2 * math.pi) - math.pi
        ws.append(round(da.sum() / (2 * math.pi), 3))
    # curvature sign of Gamma near i and j (cusp detection: turning of Gamma)
    print('  segments %d (psi=%.5f) and %d (psi=%.5f): ray-index gap %d, theta distance of the two boundary points %.3e, '
          'segment length %.2e; winding numbers around the crossing (radius half a segment): %s'
          % (i, psi[i], j, psi[j], min(j - i, n - (j - i)), np.linalg.norm(B[i] - B[j]), scale, ws))
