"""AUDIT items 6 and 7: shape-space scan in angle coordinates (independent of code).

For each triangle shape (angles A1, A2, A3 at p1, p2, p3; ccw labelling) and alpha:
 (i)   zeros of Z (= global maximisers of D, D = 1/4) by multi-start root finding, and local maximisers of D with
       D > c by multi-start maximisation: each component of M = {D > c} contains a local maximiser of D (D -> 0 at the
       boundary of T), so a single maximiser with D > c means M is connected;
 (ii)  theta* = theta(x*), x* the zero of Z.  Along 1440 rays theta* + t e(psi) in the theta-hexagon
       {A_i < theta_i < pi} (closed-form inverse x(theta), barycentrics prop. to sin A_i sin th_i / sin(th_i - A_i)),
       D is sampled on a geometric grid up to the hexagon boundary: checks that D > c exactly on an initial interval
       [0, rho(psi)) (theta(M) star-shaped w.r.t. theta*, no other piece of M on the ray) and records whether D is
       monotone on [0, rho);
 (iii) convexity of theta(M): discrete turning signs of the boundary polygon AND the curvature of the level curve
       {D = c} in theta-coordinates, kappa = -(F11 F2^2 - 2 F12 F1 F2 + F22 F1^2)/|grad F|^3 (F = D o x(theta),
       central differences with step 1e-3 rho), normalised by the mean radius: kappa_min * R_mean > 0 <=> strictly convex;
 (iv)  item 6, all charges at once: Psi = (log Q1/Q3, log Q2/Q3) is orientation preserving on M (audit a4), so for a
       charge vector q with log-ratios w not on the fold image Gamma = Psi(boundary of M), the number of planar minima
       of V_q equals the winding number of Gamma around w.  The winding number is evaluated at the images of 30 x 360
       interior points of theta(M) (each has winding >= 1); max winding >= 2 would be a counterexample (m0 = 2).
       Self-intersections of Gamma are also counted.
usage: python a67_shape_scan.py alpha step_deg [min_deg] [out_tag]"""
import math
import sys
import time

import numpy as np
from scipy.optimize import minimize, root

alpha = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
step = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
min_deg = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
beta = 2 * alpha + 2
c = (beta - 1) / beta ** 2
NPSI, NT = 1440, 160


def tri(A1, A2):
    A3 = math.pi - A1 - A2
    L13 = math.sin(A2) / math.sin(A3)
    return np.array([[0.0, 0.0], [1.0, 0.0], [L13 * math.cos(A1), L13 * math.sin(A1)]])


class Shape:
    def __init__(self, A1, A2):
        self.A = np.array([A1, A2, math.pi - A1 - A2])
        self.P = tri(A1, A2)
        self.Minv = np.linalg.inv(np.vstack([self.P.T, np.ones(3)]))

    def x_of_theta(self, T12):
        """T12 (..., 2) -> x (..., 2) and barycentrics (..., 3)"""
        th = np.concatenate([T12, (2 * math.pi - T12.sum(-1))[..., None]], -1)
        w = np.sin(self.A) * np.sin(th) / np.sin(th - self.A)
        lam = w / w.sum(-1, keepdims=True)
        return lam @ self.P, lam

    def theta_of_x(self, X):
        th = []
        for i in range(3):
            a, b = self.P[(i + 1) % 3] - X, self.P[(i + 2) % 3] - X
            cs = (a * b).sum(-1) / np.linalg.norm(a, axis=-1) / np.linalg.norm(b, axis=-1)
            sn = (a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]) / np.linalg.norm(a, axis=-1) / np.linalg.norm(b, axis=-1)
            th.append(np.arctan2(np.abs(sn), cs))
        return np.stack(th[:2], -1)

    def Z(self, X, lam=None):
        if lam is None:
            lam = np.concatenate([X, np.ones(X.shape[:-1] + (1,))], -1) @ self.Minv.T
        z = np.zeros(X.shape[:-1], complex)
        for i in range(3):
            d = X - self.P[i]
            u = (d[..., 0] + 1j * d[..., 1]) / np.linalg.norm(d, axis=-1)
            z = z + lam[..., i] * u * u
        return z

    def D(self, X, lam=None):
        return (1 - np.abs(self.Z(X, lam)) ** 2) / 4

    def Dth(self, T12):
        X, lam = self.x_of_theta(T12)
        return self.D(X, lam)

    def Psi(self, X, lam=None):
        if lam is None:
            lam = np.concatenate([X, np.ones(X.shape[:-1] + (1,))], -1) @ self.Minv.T
        L = [np.log(lam[..., i]) + beta * np.log(np.linalg.norm(X - self.P[i], axis=-1)) for i in range(3)]
        return np.stack([L[0] - L[2], L[1] - L[2]], -1)

    def bary(self, x):
        return self.Minv @ np.r_[x, 1.0]


def seeds(S, n=10):
    out = []
    for i in range(1, n):
        for j in range(1, n - i):
            k = n - i - j
            out.append(np.array([i, j, k]) / n @ S.P)
    for v in range(3):
        for t in (1e-3, 1e-2, 0.05, 0.15):
            for s in (0.25, 0.5, 0.75):
                out.append(S.P[v] + t * ((1 - s) * S.P[(v + 1) % 3] + s * S.P[(v + 2) % 3] - S.P[v]))
    return out


def analyse(A1, A2):
    S = Shape(A1, A2)
    side = 1.0
    # (i) zeros of Z
    zeros = []
    for s0 in seeds(S):
        f = lambda x: [S.Z(x[None])[0].real, S.Z(x[None])[0].imag]
        sol = root(f, s0, method='hybr', options={'xtol': 1e-14})
        x = sol.x
        if not np.all(np.isfinite(x)) or S.bary(x).min() <= 0 or abs(S.Z(x[None])[0]) > 1e-10:
            continue
        if all(np.linalg.norm(x - y) > 1e-7 for y in zeros):
            zeros.append(x)
    maxima = []
    for s0 in seeds(S, 8):
        res = minimize(lambda x: -S.D(x[None])[0] if S.bary(x).min() > 0 else 1.0, s0, method='Nelder-Mead',
                       options={'xatol': 1e-12, 'fatol': 1e-15, 'maxiter': 4000})
        x = res.x
        if S.bary(x).min() > 0 and -res.fun > c and all(np.linalg.norm(x - y) > 1e-5 for y in maxima):
            maxima.append(x)
    out = dict(nzeros=len(zeros), nmax_above_c=len(maxima))
    if len(zeros) == 0:
        out['fail'] = 'no zero of Z found'
        return out
    xs = zeros[0]
    ts = S.theta_of_x(xs[None])[0]
    # (ii) rays in theta space
    psi = np.linspace(0, 2 * math.pi, NPSI, endpoint=False)
    E = np.c_[np.cos(psi), np.sin(psi)]
    # hexagon: A1 < t1 < pi, A2 < t2 < pi, A3 < 2pi - t1 - t2 < pi
    G = np.array([[1, 0], [-1, 0], [0, 1], [0, -1], [-1, -1], [1, 1]], float)
    h = np.array([S.A[0], -math.pi, S.A[1], -math.pi, S.A[2] - 2 * math.pi, math.pi])   # G t >= h
    slack = G @ ts - h
    rate = E @ G.T                                                                # (NPSI, 6)
    with np.errstate(divide='ignore'):
        tm = np.where(rate < 0, slack[None] / -rate, np.inf).min(1)
    frac = np.geomspace(1e-7, 1 - 1e-9, NT)
    Tg = ts[None, None] + (tm[:, None, None] * frac[None, :, None]) * E[:, None, :]
    Dg = S.Dth(Tg)
    inM = Dg > c
    first_out = np.argmax(~inM, axis=1)
    if np.any(inM[:, 0] == False):
        out['fail'] = 'theta* not in M on the grid'
        return out
    reentry = int(sum(np.any(inM[k, first_out[k]:]) for k in range(NPSI)))
    nonmono = int(sum(np.any(np.diff(Dg[k, :first_out[k] + 1]) > 1e-15) for k in range(NPSI)))
    lo = tm * np.where(first_out > 0, frac[np.maximum(first_out - 1, 0)], 0.0)
    hi = tm * frac[first_out]
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        inside = S.Dth(ts[None] + mid[:, None] * E) > c
        lo = np.where(inside, mid, lo)
        hi = np.where(inside, hi, mid)
    rho = 0.5 * (lo + hi)
    B = ts[None] + rho[:, None] * E
    # (iii) discrete turning + curvature
    d1 = np.roll(B, -1, 0) - B
    cr = d1[:, 0] * np.roll(d1, -1, 0)[:, 1] - d1[:, 1] * np.roll(d1, -1, 0)[:, 0]
    turn_pos, turn_neg = int((cr > 0).sum()), int((cr < 0).sum())
    hh = 1e-3 * rho[:, None]
    F = lambda dt: S.Dth(B + dt)
    e1, e2 = np.array([1.0, 0.0]), np.array([0.0, 1.0])
    F0 = F(0 * B)
    F1 = (F(hh * e1) - F(-hh * e1)) / (2 * hh[:, 0])
    F2 = (F(hh * e2) - F(-hh * e2)) / (2 * hh[:, 0])
    F11 = (F(hh * e1) - 2 * F0 + F(-hh * e1)) / hh[:, 0] ** 2
    F22 = (F(hh * e2) - 2 * F0 + F(-hh * e2)) / hh[:, 0] ** 2
    F12 = (F(hh * (e1 + e2)) - F(hh * (e1 - e2)) - F(hh * (e2 - e1)) + F(-hh * (e1 + e2))) / (4 * hh[:, 0] ** 2)
    kappa = -(F11 * F2 ** 2 - 2 * F12 * F1 * F2 + F22 * F1 ** 2) / (F1 ** 2 + F2 ** 2) ** 1.5
    Rm = rho.mean()
    # affine-invariant version: whiten theta by the covariance of the boundary polygon (convexity is affine invariant;
    # the raw normalised curvature mostly measures elongation)
    Cv = np.cov((B - ts).T)
    ev, U = np.linalg.eigh(Cv)
    Wh = U @ np.diag(ev ** 0.5) @ U.T                       # theta - ts = Wh theta'
    g1 = Wh[0, 0] * F1 + Wh[1, 0] * F2
    g2 = Wh[0, 1] * F1 + Wh[1, 1] * F2
    Hs = np.stack([np.stack([F11, F12], -1), np.stack([F12, F22], -1)], -2)
    Hw = np.einsum('ab,nbc,cd->nad', Wh.T, Hs, Wh)
    kw = -(Hw[:, 0, 0] * g2 ** 2 - 2 * Hw[:, 0, 1] * g1 * g2 + Hw[:, 1, 1] * g1 ** 2) / (g1 ** 2 + g2 ** 2) ** 1.5
    Rw = np.mean(np.linalg.norm(np.linalg.solve(Wh, (B - ts).T), axis=0))
    # (iv) winding numbers of the fold image
    Xb, lb = S.x_of_theta(B)
    Gam = S.Psi(Xb, lb)
    sgrid = np.linspace(0.03, 0.985, 30)
    pidx = np.arange(0, NPSI, NPSI // 360)
    Ti = ts[None, None] + (rho[pidx, None, None] * sgrid[None, :, None]) * E[pidx, None, :]
    Xi, li = S.x_of_theta(Ti.reshape(-1, 2))
    W = S.Psi(Xi, li)
    wind = np.zeros(len(W))
    for a in range(0, len(W), 1500):
        dv = Gam[None, :, :] - W[a:a + 1500, None, :]
        ang = np.arctan2(dv[..., 1], dv[..., 0])
        da = np.diff(np.concatenate([ang, ang[:, :1]], 1), axis=1)
        da = (da + math.pi) % (2 * math.pi) - math.pi
        wind[a:a + 1500] = da.sum(1) / (2 * math.pi)
    wr = np.round(wind).astype(int)
    # self-intersections of Gamma (O(N^2) segment test, vectorised by rows)
    n = len(Gam)
    a1, a2 = Gam, np.roll(Gam, -1, 0)
    inter, inter_rel = 0, 0.0
    ext = np.ptp(Gam, 0).max()
    seglen = np.linalg.norm(a2 - a1, axis=1)
    for i in range(n):
        j = np.arange(i + 2, n) if i > 0 else np.arange(2, n - 1)
        if len(j) == 0:
            continue
        b1, b2 = a1[j], a2[j]
        o = lambda p, q, r: (q[..., 0] - p[..., 0]) * (r[..., 1] - p[..., 1]) - (q[..., 1] - p[..., 1]) * (r[..., 0] - p[..., 0])
        s1 = o(a1[i], a2[i], b1) * o(a1[i], a2[i], b2)
        s2 = o(b1, b2, a1[i][None]) * o(b1, b2, a2[i][None])
        hit = (s1 < 0) & (s2 < 0)
        if hit.any():
            inter += int(hit.sum())
            inter_rel = max(inter_rel, float(max(seglen[i], seglen[j[hit]].max()) / ext))
    out.update(Rx=float(np.max(np.linalg.norm(Xb - xs, axis=1))), Rtheta=float(Rm), reentry=reentry, nonmono=nonmono,
               turn=(turn_pos, turn_neg), kmin=float(kappa.min() * Rm), kmax=float(kappa.max() * Rm),
               kw_min=float(kw.min() * Rw), kw_max=float(kw.max() * Rw),
               wind_min=int(wr.min()), wind_max=int(wr.max()), wind_dev=float(np.max(np.abs(wind - wr))),
               gamma_self_int=inter, self_int_rel_seglen=inter_rel)
    return out


if __name__ == '__main__':
    tag = sys.argv[4] if len(sys.argv) > 4 else ''
    shapes = []
    g = np.arange(min_deg, 180, step)
    for A1 in g:
        for A2 in g:
            A3 = 180 - A1 - A2
            if A3 >= min_deg - 1e-9 and A1 <= A2 <= A3 + 1e-9:
                shapes.append((A1, A2))
    extra = [(0.2, 0.2), (0.2, 1.0), (0.2, 89.9), (0.5, 0.5), (1, 1), (0.5, 60), (1, 89.5), (2, 2), (0.3, 30),
             (60, 60), (59, 61), (45, 45), (30, 60), (20, 20), (10, 85)]
    shapes += [s for s in extra if s not in shapes]
    t0 = time.time()
    worst = dict(kmin=1e9, kw=1e9, kw_at=None, wind_max=0, nzeros=set(), nmax=set(), reentry=0, nonconvex=0, selfint=0,
                 selfint_big=0, wind_min=9, nonmono=0)
    for A1, A2 in shapes:
        r = analyse(math.radians(A1), math.radians(A2))
        line = 'A=(%6.2f,%6.2f,%6.2f) %s' % (A1, A2, 180 - A1 - A2, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
        print(line, flush=True)
        if 'fail' in r:
            continue
        worst['kmin'] = min(worst['kmin'], r['kmin'])
        if r['kw_min'] < worst['kw']:
            worst['kw'], worst['kw_at'] = r['kw_min'], (A1, A2)
        worst['nonmono'] += int(r['nonmono'] > 0)
        worst['selfint_big'] += int(r['self_int_rel_seglen'] > 1e-6)
        worst['wind_max'] = max(worst['wind_max'], r['wind_max'])
        worst['wind_min'] = min(worst['wind_min'], r['wind_min'])
        worst['nzeros'].add(r['nzeros'])
        worst['nmax'].add(r['nmax_above_c'])
        worst['reentry'] += r['reentry']
        worst['nonconvex'] += int(r['turn'][1] > 0 and r['turn'][0] > 0)
        worst['selfint'] += int(r['gamma_self_int'] > 0)
    print('SUMMARY alpha %.3g, %d shapes (A1<=A2<=A3, step %.2f deg, min %.2f deg, + extras) in %.0fs: #zeros of Z %s; '
          '#local maxima of D above c %s; rays re-entering M %d; shapes with D non-monotone on a ray in M %d; shapes with '
          'mixed turning %d; min normalised curvature %.2e; min WHITENED normalised curvature %.4f at %s; winding numbers '
          'of the fold image at interior points in [%d, %d]; shapes with self-intersecting fold image %d (of which with '
          'crossing segments longer than 1e-6 x extent: %d)'
          % (alpha, len(shapes), step, min_deg, time.time() - t0, sorted(worst['nzeros']), sorted(worst['nmax']),
             worst['reentry'], worst['nonmono'], worst['nonconvex'], worst['kmin'], worst['kw'], worst['kw_at'],
             worst['wind_min'], worst['wind_max'], worst['selfint'], worst['selfint_big']), flush=True)
