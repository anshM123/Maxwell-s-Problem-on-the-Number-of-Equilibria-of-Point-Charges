"""AUDIT_3, item 7: direct counterexample search at alpha = 1/4, 1, 2 (no part of the reduction is used to count).

For each alpha, mixed-sign configurations in the plane (q1, q2 > 0 > q3 after relabelling; triangles generic and
near-degenerate), of two kinds:
  'forced'  a minimum is forced at x(t) for t drawn in M' = {E < K} (E from its definition), with x(t) computed by the
            auditor's inverse map (inversion at p3 turns both angle loci into rays), q_i = lambda_i(x) r_i(x)^beta;
            half of the t are drawn near the boundary of M' (K - E small), where a second extremum would appear first;
  'random'  random charges (magnitudes over four orders), no forced point.
All critical points of V|plane are searched by damped Newton iterations from ~3,700 seeds (log-polar grids around the
centroid out to 1e6 triangle sizes, plus the forced point and a cloud around it and around each charge), deduplicated and classified by the Hessian. Flags: more than
one extremum, more than four critical points, a maximum, or an index count that fails (#ext - #sad must be -2, or -1 when
the total charge is 0). Flagged configurations are re-solved at 50 digits (mpmath Newton from the float points and the forced point;
only points with |grad V| < 1e-35 x scale and a final Newton step < 1e-30 count) before being reported.
Output: audit/c6_direct_search.out (details of still-flagged cases: audit/c6b_flagged_detail.out)"""
import math
import sys
import time

import mpmath as mp
import numpy as np

rng = np.random.default_rng(424242)


def triangle(A1, A2):
    A3 = math.pi - A1 - A2
    L = math.sin(A2) / math.sin(A3)
    return np.array([[0.0, 0.0], [1.0, 0.0], [L * math.cos(A1), L * math.sin(A1)]])


def E_t(A1, A2, t1, t2):
    cot = lambda x: np.cos(x) / np.sin(x)
    U1 = cot(t1) + cot(A1 - t1)
    U2 = cot(t2) + cot(A2 - t2)
    U3 = cot(t1 + t2) + cot(A1 + A2 - t1 - t2)
    return U1 * U2 - U3 * (U1 + U2)


def x_of_t(P, t1, t2):
    """x with sigma1 = arg((p3 - x)/(p2 - x)) = t1 and sigma2 = arg((p1 - x)/(p3 - x)) = t2, via w = 1/(x - p3):
    1 - d2 w = rho e^{-i t1}, 1 - d1 w = rho' e^{i t2}, d_i = p_i - p3"""
    p = [complex(*P[i]) for i in range(3)]
    d1, d2 = p[0] - p[2], p[1] - p[2]
    # (1 - rho e^{-it1})/d2 = (1 - rho' e^{it2})/d1  ->  rho e^{-it1}/d2 - rho' e^{it2}/d1 = 1/d2 - 1/d1
    a = np.exp(-1j * t1) / d2
    b = -np.exp(1j * t2) / d1
    c = 1 / d2 - 1 / d1
    M = np.array([[a.real, b.real], [a.imag, b.imag]])
    rho, rhop = np.linalg.solve(M, np.array([c.real, c.imag]))
    w = (1 - rho * np.exp(-1j * t1)) / d2
    x = p[2] + 1 / w
    return np.array([x.real, x.imag]), rho, rhop


def bary(P, x):
    def cr(u, v):
        return u[0] * v[1] - u[1] * v[0]
    D2 = cr(P[1] - P[0], P[2] - P[0])
    return np.array([cr(P[(i + 1) % 3] - x, P[(i + 2) % 3] - x) / D2 for i in range(3)])


def grad_hess(P, q, alpha, X):
    """X: (n, 2) -> grad (n, 2), Hess (n, 2, 2) of V = sum q_i |x - p_i|^(-2 alpha)"""
    beta = 2 * alpha + 2
    g = np.zeros_like(X)
    H = np.zeros((X.shape[0], 2, 2))
    for i in range(3):
        d = X - P[i]
        r2 = (d * d).sum(1)
        w = q[i] * r2 ** (-beta / 2)
        g += -2 * alpha * w[:, None] * d
        H += -2 * alpha * w[:, None, None] * (np.eye(2)[None] - beta * d[:, :, None] * d[:, None, :] / r2[:, None, None])
    return g, H


def find_critical(P, q, alpha, extra=None):
    c = P.mean(0)
    size = max(np.linalg.norm(P[i] - P[j]) for i in range(3) for j in range(i))
    seeds = []
    th = np.linspace(0, 2 * np.pi, 48, endpoint=False)
    for R in size * 10 ** np.linspace(-2, 6, 52):
        seeds.append(c + R * np.stack([np.cos(th), np.sin(th)], 1))
    dmin = min(np.linalg.norm(P[i] - P[j]) for i in range(3) for j in range(i))
    th2 = np.linspace(0, 2 * np.pi, 24, endpoint=False)
    for i in range(3):
        for R in dmin * 10 ** np.linspace(-3, -0.1, 14):
            seeds.append(P[i] + R * np.stack([np.cos(th2), np.sin(th2)], 1))
    if extra is not None:                      # the forced point and a small cloud around it
        seeds.append(extra[None, :] + np.linalg.norm(extra - c) * 1e-3 * rng.normal(size=(20, 2)))
        seeds.append(extra[None, :])
    X = np.vstack(seeds)
    for it in range(150):
        g, H = grad_hess(P, q, alpha, X)
        det = H[:, 0, 0] * H[:, 1, 1] - H[:, 0, 1] * H[:, 1, 0]
        det = np.where(np.abs(det) < 1e-300, 1e-300, det)
        step = np.stack([(H[:, 1, 1] * g[:, 0] - H[:, 0, 1] * g[:, 1]) / det,
                         (-H[:, 1, 0] * g[:, 0] + H[:, 0, 0] * g[:, 1]) / det], 1)
        dist = np.min(np.stack([np.linalg.norm(X - P[i], axis=1) for i in range(3)], 1), 1)
        lim = 0.3 * np.maximum(dist, 1e-12)
        sn = np.linalg.norm(step, axis=1)
        fac = np.where(sn > lim, lim / np.maximum(sn, 1e-300), 1.0)
        X = X - step * fac[:, None]
        X = np.where(np.isfinite(X), X, c)
    g, H = grad_hess(P, q, alpha, X)
    dist = np.stack([np.linalg.norm(X - P[i], axis=1) for i in range(3)], 1)
    scale = sum(np.abs(q[i]) * dist[:, i] ** (-2 * alpha - 1) for i in range(3))
    ok = ((np.linalg.norm(g, axis=1) < 1e-10 * scale) & (dist.min(1) > 1e-9 * size)
          & (np.linalg.norm(X - c, axis=1) < 1e8 * size))
    pts = []
    for x in X[ok]:
        if not any(np.linalg.norm(x - y) < 1e-7 * (size + np.linalg.norm(y - c)) for y in pts):
            pts.append(x)
    pts = np.array(pts).reshape(-1, 2)
    kinds = []
    if len(pts):
        _, H = grad_hess(P, q, alpha, pts)
        for h in H:
            ev = np.linalg.eigvalsh(h)
            kinds.append('min' if ev[0] > 0 else 'max' if ev[1] < 0 else 'saddle' if ev[0] * ev[1] < 0 else 'degenerate')
    return pts, kinds


def mp_verify(P, q, alpha, pts):
    """50-digit Newton from the float points; returns the set of distinct refined critical points and their types"""
    mp.mp.dps = 50
    Pm = [mp.matrix([mp.mpf(P[i][0]), mp.mpf(P[i][1])]) for i in range(3)]
    qm = [mp.mpf(v) for v in q]
    a = mp.mpf(alpha)
    beta = 2 * a + 2

    def gh(x):
        g = mp.matrix([0, 0])
        H = mp.matrix(2, 2)
        for i in range(3):
            d = x - Pm[i]
            r2 = d[0] ** 2 + d[1] ** 2
            w = qm[i] * r2 ** (-beta / 2)
            g += -2 * a * w * d
            for j in range(2):
                for k in range(2):
                    H[j, k] += -2 * a * w * ((1 if j == k else 0) - beta * d[j] * d[k] / r2)
        return g, H
    out = []
    for p0 in pts:
        x = mp.matrix([mp.mpf(p0[0]), mp.mpf(p0[1])])
        step = mp.mpf(1)
        try:
            for _ in range(80):
                g, H = gh(x)
                dx = mp.lu_solve(H, g)
                step = mp.norm(dx)
                x = x - dx
        except (ZeroDivisionError, ValueError):
            continue
        g, H = gh(x)
        sc = sum(abs(qm[i]) * mp.norm(x - Pm[i]) ** (-2 * a - 1) for i in range(3))
        # count only genuine (converged) critical points
        if not (mp.norm(g) < mp.mpf(10) ** -35 * sc and step < mp.mpf(10) ** -30 * (1 + mp.norm(x))):
            continue
        ev = mp.eigsy(H)[0]
        kind = 'min' if ev[0] > 0 else 'max' if ev[1] < 0 else 'saddle'
        if not any(mp.norm(x - y) < mp.mpf(10) ** -25 * (1 + mp.norm(y)) for y, _ in out):
            out.append((x, kind))
    return out


def shape_sample(fam):
    if fam == 'generic':
        A1 = rng.uniform(0.1, 2.9)
        A2 = rng.uniform(0.05, math.pi - A1 - 0.05)
    elif fam == 'thinA12':
        e = 10 ** rng.uniform(-3, -1)
        lam = rng.uniform(0.05, 0.95)
        A1, A2 = lam * e, (1 - lam) * e
    elif fam == 'flatA3':
        A3 = 10 ** rng.uniform(-3, -1)
        lam = rng.uniform(0.05, 0.95)
        A1, A2 = lam * (math.pi - A3), (1 - lam) * (math.pi - A3)
    else:  # 'cornerA1'
        r = 10 ** rng.uniform(-3, -1)
        th = rng.uniform(0.05, 0.95)
        A1, A2 = math.pi - r, th * r
    return A1, A2


def main():
    n_per = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    lines = []
    alphas = [float(a) for a in sys.argv[2].split(',')] if len(sys.argv) > 2 else (0.25, 1.0, 2.0)
    for alpha in alphas:
        b = 2 * alpha + 1
        K = (b + 1) ** 2 / b
        beta = b + 1
        tally, flagged, t0 = {}, 0, time.time()
        for kind in ('forced', 'random'):
            for fam in ('generic', 'thinA12', 'flatA3', 'cornerA1'):
                for _ in range(n_per):
                    A1, A2 = shape_sample(fam)
                    P = triangle(A1, A2)
                    if kind == 'forced':
                        # t in M' (half of them close to the boundary)
                        near = rng.random() < 0.5
                        found = False
                        for _ in range(20000):
                            sc = 10 ** rng.uniform(-4, 0)
                            t1 = A1 / 2 * (1 + sc * rng.uniform(-1, 1))
                            t2 = A2 / 2 * (1 + sc * rng.uniform(-1, 1))
                            E = E_t(A1, A2, t1, t2)
                            if E < K and (not near or K - E < 0.05 * (K - 4)):
                                found = True
                                break
                        if not found:
                            tally['no t found in M-prime'] = tally.get('no t found in M-prime', 0) + 1
                            continue
                        x, _, _ = x_of_t(P, t1, t2)
                        lam = bary(P, x)
                        r = np.linalg.norm(x - P, axis=1)
                        q = lam * r ** beta
                    else:
                        q = np.array([10 ** rng.uniform(-2, 2), 10 ** rng.uniform(-2, 2), -10 ** rng.uniform(-2, 2)])
                    if kind == 'forced':
                        g0, H0 = grad_hess(P, q, alpha, x[None, :])
                        if np.linalg.eigvalsh(H0[0])[0] <= 0:
                            tally['forced point not a minimum in float'] = tally.get('forced point not a minimum in float', 0) + 1
                            continue
                    if not (q[0] > 0 and q[1] > 0 and q[2] < 0):
                        tally['bad construction'] = tally.get('bad construction', 0) + 1
                        continue
                    pts, kinds = find_critical(P, q, alpha, x if kind == 'forced' else None)
                    n_ext = sum(k in ('min', 'max') for k in kinds)
                    n_sad = kinds.count('saddle')
                    Qt = q.sum()
                    want = -1 if abs(Qt) < 1e-12 * np.abs(q).max() else -2
                    flag = (n_ext > 1 or len(pts) > 4 or 'max' in kinds or n_ext - n_sad != want or
                            'degenerate' in kinds or (kind == 'forced' and n_ext != 1))
                    if flag:
                        ref = mp_verify(P, q, alpha, list(pts) + ([x] if kind == 'forced' else []))
                        rk = [k for _, k in ref]
                        n_ext2 = sum(k in ('min', 'max') for k in rk)
                        still = n_ext2 > 1 or len(ref) > 4 or 'max' in rk or (kind == 'forced' and n_ext2 != 1)
                        key = '%s/%s flagged in float (%d ext, %d sad) -> 50 digits: %d ext, %d sad%s' % (
                            kind, fam, n_ext, n_sad, n_ext2, rk.count('saddle'), ' STILL FLAGGED' if still else '')
                        tally[key] = tally.get(key, 0) + 1
                        flagged += still
                        if still:
                            sc = sum(abs(q[i]) * np.linalg.norm(x - P[i]) ** (-2 * alpha - 1) for i in range(3)) if kind == 'forced' else 0
                            g0, H0 = grad_hess(P, q, alpha, x[None, :]) if kind == 'forced' else (None, None)
                            print('   DETAIL %s/%s: A1 %r A2 %r; q %r; ' % (kind, fam, A1, A2, list(q))
                                  + ('t %r %r, E %r, K %r, (K-E)/(K-4) %.3e; forced x %r, |grad|/scale %.2e, Hess eig %r; '
                                     % (t1, t2, E, K, (K - E) / (K - 4), list(x), np.linalg.norm(g0) / sc,
                                        list(np.linalg.eigvalsh(H0[0]))) if kind == 'forced' else '')
                                  + '50-digit points: %s' % [(mp.nstr(y[0], 12), mp.nstr(y[1], 12), k) for y, k in ref],
                                  flush=True)
                    else:
                        key = '%s (%d ext, %d sad)' % (kind, n_ext, n_sad)
                        tally[key] = tally.get(key, 0) + 1
        lines.append('alpha %.2f (K = %.4f): %d configurations per kind (4 shape families), %.0fs; outcomes:'
                     % (alpha, K, 4 * n_per, time.time() - t0))
        for k in sorted(tally):
            lines.append('    %-60s %d' % (k, tally[k]))
        lines.append('    configurations still flagged after 50-digit refinement: %d' % flagged)
        print('\n'.join(lines[-(len(tally) + 2):]), flush=True)


if __name__ == '__main__':
    main()
