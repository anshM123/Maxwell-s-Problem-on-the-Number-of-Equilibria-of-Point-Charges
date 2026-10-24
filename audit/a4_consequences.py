"""AUDIT item 4: the angle map and the definiteness of Hess_theta H = DPsi A^{-1} on M.

 (4a) injectivity of x -> theta and its image.  Closed-form inverse (classical): the point with subtended angles
      theta_i has barycentric coordinates proportional to 1/(cot A_i - cot theta_i), A_i the triangle angle at p_i.
      Checked: (i) the formula reproduces x from theta(x) for random x in T (so theta is injective);
      (ii) for random theta in the hexagon {A_i < theta_i < pi, sum = 2 pi} the formula gives a point of T whose
      angles are theta (so theta(T) is exactly this convex hexagon).
 (4b) det A (A = D_x theta, unsigned angles) never vanishes; its sign is the sign of the labelling (ccw +, cw -), and
      det A = cross(e1,e2) + cross(e2,e3) + cross(e3,e1) (ccw), e_a = (x - p_a)/r_a^2, i.e. twice the area of the
      triangle of the inverted vertices: positive for ccw.  Points near edges and vertices included.
 (4c) THE SIGN THAT PROOF.md DOES NOT ESTABLISH: det Hess_theta H = det DPsi / det A has the sign of (D - c) (not of
      -(D - c)), for both orientations; on M both eigenvalues of the (symmetric) Hess_theta H are positive.
usage: python a4_consequences.py"""
import math

import numpy as np

rng = np.random.default_rng(12)


def tri_from_angles(A, B, orient=1):
    C = math.pi - A - B
    L13 = math.sin(B) / math.sin(C)
    P = np.array([[0.0, 0.0], [1.0, 0.0], [L13 * math.cos(A), L13 * math.sin(A)]])
    if orient == -1:
        P[:, 1] *= -1
    return P


def angles_of_triangle(P):
    out = []
    for i in range(3):
        a, b = P[(i + 1) % 3] - P[i], P[(i + 2) % 3] - P[i]
        out.append(math.acos(np.clip(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b), -1, 1)))
    return np.array(out)


def thetas(P, x):
    out = []
    for i in range(3):
        a, b = P[(i + 1) % 3] - x, P[(i + 2) % 3] - x
        out.append(math.acos(np.clip(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b), -1, 1)))
    return np.array(out)


def bary(P, x):
    return np.linalg.solve(np.vstack([P.T, np.ones(3)]), np.r_[x, 1.0])


def x_of_theta(P, th):
    Ai = angles_of_triangle(P)
    w = 1.0 / (1.0 / np.tan(Ai) - 1.0 / np.tan(th))
    return (w / w.sum()) @ P, w


def jac_fd(f, x, h):
    J = np.zeros((2, 2))
    for k in range(2):
        e = np.zeros(2); e[k] = h
        J[:, k] = (f(x + e) - f(x - e)) / (2 * h)
    return J


def DPsi(P, x, beta):
    Minv = np.linalg.inv(np.vstack([P.T, np.ones(3)]))
    lam = Minv @ np.r_[x, 1.0]
    gl = Minv[:, :2]                       # grad lambda_i = row i
    G = [gl[i] / lam[i] + beta * (x - P[i]) / np.dot(x - P[i], x - P[i]) for i in range(3)]
    return np.array([G[0] - G[2], G[1] - G[2]])


def Dval(P, x):
    lam = bary(P, x)
    U = [(x - p) / np.linalg.norm(x - p) for p in P]
    S = sum(lam[i] * np.outer(U[i], U[i]) for i in range(3))
    return np.linalg.det(S)


def rand_shape(rng, min_deg):
    while True:
        w = rng.dirichlet(np.ones(3)) * math.pi
        if w.min() > math.radians(min_deg):
            return w


def rand_point(P, rng, kind):
    if kind == 'any':
        w = rng.dirichlet(np.ones(3))
    elif kind == 'edge':
        w = rng.dirichlet(np.ones(3)); w[rng.integers(3)] = 10 ** rng.uniform(-6, -2); w /= w.sum()
    else:
        w = np.zeros(3); v = rng.integers(3); w[v] = 1.0
        o = rng.dirichlet(np.ones(2)) * 10 ** rng.uniform(-6, -2)
        w[(v + 1) % 3], w[(v + 2) % 3] = o
        w /= w.sum()
    return w @ P


# ------------------------------------------------------------------ (4a)
worst_inv, worst_img, outside = {'any': 0.0, 'edge': 0.0, 'vertex': 0.0}, 0.0, 0
for k in range(3000):
    sh = rand_shape(rng, 0.5)
    P = tri_from_angles(sh[0], sh[1], 1 if k % 2 else -1)
    kind = ['any', 'edge', 'vertex'][k % 3]
    x = rand_point(P, rng, kind)
    xr, _ = x_of_theta(P, thetas(P, x))
    worst_inv[kind] = max(worst_inv[kind], np.linalg.norm(xr - x) / np.linalg.norm(P[1] - P[0]))
    # random theta in the hexagon
    Ai = angles_of_triangle(P)
    while True:
        t12 = rng.uniform(Ai[:2], math.pi)
        th = np.r_[t12, 2 * math.pi - t12.sum()]
        if Ai[2] < th[2] < math.pi:
            break
    xt, w = x_of_theta(P, th)
    if bary(P, xt).min() <= 0:
        outside += 1
    else:
        worst_img = max(worst_img, np.max(np.abs(thetas(P, xt) - th)))
print('(4a) inverse formula reproduces x from theta(x): max error by point type %s (float64; near edges/vertices the cot difference cancels); '
      'random theta in the hexagon {A_i < theta_i < pi}: %d of 3000 mapped outside T, max |theta(x(theta)) - theta| %.2e'
      % ({k: '%.1e' % v for k, v in worst_inv.items()}, outside, worst_img))

# ------------------------------------------------------------------ (4b), (4c)
signs = {1: set(), -1: set()}
worst_formula, worst_fd = 0.0, 0.0
mism, onM, offM, posdef, tr_signs = 0, 0, 0, 0, set()
min_ratio = 1e9
for k in range(20000):
    orient = 1 if k % 2 else -1
    sh = rand_shape(rng, 0.5 if k % 4 < 2 else 20)
    P = tri_from_angles(sh[0], sh[1], orient)
    x = rand_point(P, rng, ['any', 'any', 'edge', 'vertex'][k % 4] if k % 5 else 'any')
    e = [(x - p) / np.dot(x - p, x - p) for p in P]
    cr = lambda a, b: a[0] * b[1] - a[1] * b[0]
    # analytic A from the unsigned-angle gradient: grad theta_i = s J (e_k - e_j), s = orientation sign
    Jr = lambda v: np.array([-v[1], v[0]])
    A = np.array([orient * Jr(e[2] - e[1]), orient * Jr(e[0] - e[2])])
    detA = np.linalg.det(A)
    signs[orient].add(int(np.sign(detA)))
    worst_formula = max(worst_formula, abs(detA - (cr(e[0], e[1]) + cr(e[1], e[2]) + cr(e[2], e[0]))) / abs(detA))
    if k % 50 == 0 and k % 4 < 2:      # finite-difference cross-check of A (well inside T)
        h = 1e-6 * np.linalg.norm(P[1] - P[0])
        Afd = jac_fd(lambda y: thetas(P, y)[:2], x, h)
        worst_fd = max(worst_fd, np.max(np.abs(Afd - A)) / np.max(np.abs(A)))
    alpha = [0.5, 0.5, 0.25, 1.0, 2.5][k % 5]
    beta = 2 * alpha + 2
    c = (beta - 1) / beta ** 2
    Dv = Dval(P, x)
    if abs(Dv - c) < 1e-9:
        continue
    Ht = DPsi(P, x, beta) @ np.linalg.inv(A)
    dH = np.linalg.det(Ht)
    if np.sign(dH) != np.sign(Dv - c):
        mism += 1
    if Dv > c:
        onM += 1
        ev = np.linalg.eigvalsh((Ht + Ht.T) / 2)
        tr_signs.add(int(np.sign(np.trace(Ht))))
        if ev.min() > 0:
            posdef += 1
            min_ratio = min(min_ratio, ev.min() / ev.max())
    else:
        offM += 1
print('(4b) sign(det A): ccw-labelled %s, cw-labelled %s (never 0);  |det A - sum_cyc cross(e_i,e_j)| / |det A| max %.1e; '
      'finite-difference check of A max rel err %.1e' % (sorted(signs[1]), sorted(signs[-1]), worst_formula, worst_fd))
print('(4c) sign(det Hess_theta H) != sign(D - c) at %d of %d points (both orientations, alpha in {.25,.5,1,2.5});'
      ' on M (%d points): Hess_theta H positive definite at %d, trace signs %s, min eigenvalue ratio %.3g'
      % (mism, onM + offM, onM, posdef, sorted(tr_signs), min_ratio))
