"""Hyperbolic form of the minimum region.

X = cot(delta), Y = cot(gamma) (delta_i = theta_i - A_i pedal angles, gamma_i = pi - theta_i) lie on the hyperboloid
e2 = 1 (the cot-model of triangle shapes = Poincare half-plane), and with d = bX - Y (b = beta - 1),
    e2(d) = b^2 + 1 - 2b <X, Y>,   <X, Y> = ((sum X)(sum Y) - sum X_i Y_i)/2 = cosh rho(delta, gamma).
So M = {rho(delta, gamma) < log b}.  With m = (pi - A)/2 (the incentre: delta = gamma = m) and delta = m + u,
gamma = m - u, the region is {Xi(u) < (b + 1/b)/2}, Xi(u) = <cot(m + u), cot(m - u)>.
Checks: 1. e2(d) = b^2 + 1 - 2b Xi at random points; 2. is Xi convex on the hexagon |u_i| < m_i (sum u = 0)?
usage: python hyperbolic.py [n_shapes] [min_angle_deg] [b]"""
import math
import sys

import numpy as np

nsh = int(sys.argv[1]) if len(sys.argv) > 1 else 200
amin = math.radians(float(sys.argv[2])) if len(sys.argv) > 2 else math.radians(1.0)
b = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0


def Xi_parts(u, m):
    """Xi and its Hessian (in R^3) at u (shape (3, ...)), m (3,)"""
    mm = m.reshape((3,) + (1,) * (u.ndim - 1))
    f, h = 1 / np.tan(mm + u), 1 / np.tan(mm - u)
    f1, h1 = -1 / np.sin(mm + u) ** 2, 1 / np.sin(mm - u) ** 2
    f2, h2 = 2 * f / np.sin(mm + u) ** 2, 2 * h / np.sin(mm - u) ** 2
    Xi = 0.5 * (f.sum(0) * h.sum(0) - (f * h).sum(0))
    Hs = np.empty((3, 3) + u.shape[1:])
    for i in range(3):
        oth = [j for j in range(3) if j != i]
        Hs[i, i] = 0.5 * (f2[i] * h[oth].sum(0) + h2[i] * f[oth].sum(0))
        for j in oth:
            Hs[i, j] = 0.5 * (f1[i] * h1[j] + f1[j] * h1[i])
    return Xi, Hs


if __name__ == '__main__':
    rng = np.random.default_rng(3)
    # 1. identity check
    err = 0.0
    for k in range(2000):
        A = rng.dirichlet(np.ones(3)) * math.pi
        th = A + rng.dirichlet(np.ones(3)) * math.pi         # delta = th - A in the simplex
        if np.any(th >= math.pi):
            continue
        d = b / np.tan(th - A) + 1 / np.tan(th)
        X, Y = 1 / np.tan(th - A), 1 / np.tan(math.pi - th)
        xi = 0.5 * (X.sum() * Y.sum() - (X * Y).sum())
        e2 = d[0] * d[1] + d[1] * d[2] + d[2] * d[0]
        err = max(err, abs(e2 - (b * b + 1 - 2 * b * xi)) / (1 + abs(e2)))
    print('1. e2(d) = b^2 + 1 - 2b Xi: max relative error %.2e' % err)
    # 2. convexity of Xi on the hexagon (and on the region Xi < (b + 1/b)/2)
    E = np.array([[1, -1, 0], [1, 1, -2]], dtype=float)
    E[0] /= math.sqrt(2)
    E[1] /= math.sqrt(6)
    kb = (b + 1 / b) / 2
    worst_all, worst_P, worst_shape = np.inf, np.inf, None
    for k in range(nsh):
        while True:
            A = rng.dirichlet(np.ones(3) * 0.6) * math.pi
            if A.min() >= amin:
                break
        m = (math.pi - A) / 2
        n = 300
        s1 = np.linspace(-1, 1, n)
        S1, S2 = np.meshgrid(s1, s1)
        u = S1[None] * E[0][:, None, None] * 2 * m.max() + S2[None] * E[1][:, None, None] * 2 * m.max()
        ok = np.all(np.abs(u) < 0.999 * m[:, None, None], axis=0)
        u = u[:, ok]
        Xi, Hs = Xi_parts(u, m)
        # restricted Hessian, 2x2, normalised by Xi - 1 + small (scale) -> smallest eigenvalue
        R = np.einsum('ai,ijn,bj->abn', E, Hs, E)
        tr = R[0, 0] + R[1, 1]
        det = R[0, 0] * R[1, 1] - R[0, 1] * R[1, 0]
        lmin = tr / 2 - np.sqrt(np.maximum(tr * tr / 4 - det, 0))
        rel = lmin / (np.abs(tr) + 1e-300)
        worst_all = min(worst_all, rel.min())
        inP = Xi < kb
        if inP.any() and rel[inP].min() < worst_P:
            worst_P, worst_shape = rel[inP].min(), A.copy()
    print('2. min over hexagon of (smallest eigenvalue / trace) of Hess Xi: %.4f' % worst_all)
    print('   restricted to {Xi < %.4f}: %.4f at angles %s deg' % (kb, worst_P,
          np.round(np.degrees(worst_shape), 2).tolist() if worst_shape is not None else None))
