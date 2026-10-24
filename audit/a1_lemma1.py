"""AUDIT item 1 (independent code, nothing imported from code).

Checks Lemma 1 of PROOF.md and the quoted standard facts:
 (1a) gradient and Hessian of V(x) = sum q_i |x-p_i|^(-2 alpha) in R^3 against 40-digit finite differences of V;
 (1b) barycentric step: for x in T and q = Q(x) = (lambda_i r_i^beta), grad V_q(x) = 0 (to 1e-35);
      converse: every numerically found critical point x of a random V_q has Q(x) proportional to q;
 (1c) S = sum lambda_i u_i u_i^T: trace S = 1, det S = sum_{i<j} lambda_i lambda_j (u_i x u_j)^2 = (1-|Z|^2)/4,
      det(I - beta S) = beta^2 (D - c), trace(I - beta S) = 2 - beta, and the type (min/saddle) from the Hessian
      eigenvalues equals sign(D - c); M = {D > c} = {|Z| < alpha/(alpha+1)};
 (1d) standard facts: (i) 3D Newton from off-plane seeds only finds critical points in the plane inside T, and the
      3D Hessian at a planar critical point has normal eigenvalue -2 alpha sum_i q_i r_i^-beta < 0;
      (ii) the planar Laplacian of V equals 4 alpha^2 sum q_i r_i^(-2 alpha - 2) > 0;
      (iii) Poincare-Hopf bookkeeping: winding number of grad V on a large circle (+1) minus the windings on small
      circles around the charges (+1 each) gives total index -2 for the critical points; checked against the
      minima/saddles found, i.e. #saddles - #minima = 2.
usage: python a1_lemma1.py [seed]"""
import math
import sys

import mpmath as mp
import numpy as np
from scipy.optimize import root

mp.mp.dps = 40
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261024
rng = np.random.default_rng(seed)


# ---------------------------------------------------------------- basic objects (written from scratch)
def rand_triangle_2d(rng, min_angle_deg=4.0):
    while True:
        P = rng.normal(size=(3, 2))
        ang = []
        for i in range(3):
            a, b = P[(i + 1) % 3] - P[i], P[(i + 2) % 3] - P[i]
            ang.append(math.degrees(math.acos(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b))))
        if min(ang) > min_angle_deg:
            return P


def bary2(P, x):
    A = np.vstack([P.T, np.ones(3)])
    return np.linalg.solve(A, np.array([x[0], x[1], 1.0]))


def V_mp(x, P, q, alpha):
    s = mp.mpf(0)
    for p, qi in zip(P, q):
        r2 = sum((xk - pk) ** 2 for xk, pk in zip(x, p))
        s += qi * r2 ** (-alpha)
    return s


def grad_formula(x, P, q, alpha):
    beta = 2 * alpha + 2
    g = np.zeros(len(x))
    for p, qi in zip(P, q):
        d = x - p
        r = np.linalg.norm(d)
        g += -2 * alpha * qi * r ** (-beta) * d
    return g


def hess_formula(x, P, q, alpha):
    beta = 2 * alpha + 2
    n = len(x)
    H = np.zeros((n, n))
    for p, qi in zip(P, q):
        d = x - p
        r = np.linalg.norm(d)
        u = d / r
        H += -2 * alpha * qi * r ** (-beta) * (np.eye(n) - beta * np.outer(u, u))
    return H


def fd_grad_hess_mp(x, P, q, alpha, h=mp.mpf('1e-12')):
    n = len(x)
    xm = [mp.mpf(float(v)) for v in x]
    Pm = [[mp.mpf(float(v)) for v in p] for p in P]
    qm = [mp.mpf(float(v)) for v in q]
    am = mp.mpf(alpha)
    g = []
    for k in range(n):
        xp, xn = list(xm), list(xm)
        xp[k] += h
        xn[k] -= h
        g.append((V_mp(xp, Pm, qm, am) - V_mp(xn, Pm, qm, am)) / (2 * h))
    hh = mp.mpf('1e-10')
    H = [[None] * n for _ in range(n)]
    for k in range(n):
        for l in range(n):
            vals = []
            for sk, sl in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                y = list(xm)
                y[k] += sk * hh
                y[l] += sl * hh
                vals.append(V_mp(y, Pm, qm, am))
            H[k][l] = (vals[0] - vals[1] - vals[2] + vals[3]) / (4 * hh * hh)
    return np.array([float(v) for v in g]), np.array([[float(v) for v in row] for row in H])


def S_D_Z(P, x):
    lam = bary2(P, x)
    U = [(x - p) / np.linalg.norm(x - p) for p in P]
    S = sum(lam[i] * np.outer(U[i], U[i]) for i in range(3))
    cr = lambda a, b: a[0] * b[1] - a[1] * b[0]
    Dcb = lam[0] * lam[1] * cr(U[0], U[1]) ** 2 + lam[1] * lam[2] * cr(U[1], U[2]) ** 2 + lam[2] * lam[0] * cr(U[2], U[0]) ** 2
    Z = sum(lam[i] * complex(U[i][0], U[i][1]) ** 2 for i in range(3))
    return lam, S, Dcb, Z


def critical_points_2d(P, q, alpha, nseed=24):
    """zeros of F(x) = x - sum_i w_i p_i / sum_i w_i, w_i = q_i r_i^-beta (a bounded, well-scaled form of grad V = 0)
    from seeds on a barycentric grid and geometric clusters at the vertices; MINPACK hybrd; deduplicated."""
    beta = 2 * alpha + 2

    def F(x):
        # unit-scale form of grad V: sum w_i d_i / sum w_i r_i (|F| -> 1 at the vertices, no spurious zeros there)
        r = np.array([np.linalg.norm(x - p) for p in P])
        w = q * r ** (-beta)
        return (w[:, None] * (x - P)).sum(0) / (w * r).sum()

    seeds = []
    for i in range(1, nseed):
        for j in range(1, nseed - i):
            k = nseed - i - j
            seeds.append((i * P[0] + j * P[1] + k * P[2]) / nseed)
    for v in range(3):
        for t in np.geomspace(1e-5, 0.3, 12):
            for s in np.linspace(0.1, 0.9, 7):
                seeds.append(P[v] + t * ((1 - s) * P[(v + 1) % 3] + s * P[(v + 2) % 3] - P[v]))
    found = []
    for s0 in seeds:
        sol = root(F, s0, method='hybr', options={'xtol': 1e-15})
        x = sol.x
        if not np.all(np.isfinite(x)):
            continue
        lam = bary2(P, x)
        if lam.min() <= 0:
            continue
        if np.linalg.norm(F(x)) > 1e-12:
            continue
        if all(np.linalg.norm(x - y) > 1e-7 for y in found):
            found.append(x)
    return found


def winding(fun, curve):
    """winding number of the planar vector field fun along the closed polygon curve (N,2)"""
    v = np.array([fun(c) for c in curve])
    a = np.arctan2(v[:, 1], v[:, 0])
    da = np.diff(np.concatenate([a, a[:1]]))
    da = (da + math.pi) % (2 * math.pi) - math.pi
    assert np.max(np.abs(da)) < 0.5, 'curve sampling too coarse for the winding number'
    return da.sum() / (2 * math.pi)


# ---------------------------------------------------------------- (1a) gradient and Hessian in R^3
print('(1a) gradient/Hessian formulas vs 40-digit finite differences of V, random 3D points (off-plane too)')
worst_g, worst_h = 0.0, 0.0
for k in range(60):
    alpha = [0.5, 0.25, 1.0, 2.5, 0.1][k % 5]
    P3 = rng.normal(size=(3, 3))
    q = np.exp(rng.normal(size=3))
    x = rng.normal(size=3)
    if min(np.linalg.norm(x - p) for p in P3) < 0.2:
        continue
    g_fd, H_fd = fd_grad_hess_mp(x, P3, q, alpha)
    g, H = grad_formula(x, P3, q, alpha), hess_formula(x, P3, q, alpha)
    worst_g = max(worst_g, np.max(np.abs(g - g_fd)) / np.max(np.abs(g_fd)))
    worst_h = max(worst_h, np.max(np.abs(H - H_fd)) / np.max(np.abs(H_fd)))
print('     max relative error: gradient %.2e, Hessian %.2e' % (worst_g, worst_h))

# ---------------------------------------------------------------- (1b) barycentric step
print('(1b) q = Q(x) = lambda_i r_i^beta makes x critical (40-digit residual), random triangles/points/alpha')
worst = mp.mpf(0)
for k in range(200):
    alpha = [0.5, 0.25, 1.0, 2.5, 0.1][k % 5]
    beta = 2 * alpha + 2
    P = rand_triangle_2d(rng)
    lam = rng.dirichlet(np.ones(3))
    x = lam @ P
    Pm = [[mp.mpf(float(v)) for v in p] for p in P]
    xm = [sum(mp.mpf(float(lam[i])) * Pm[i][c] for i in range(3)) for c in range(2)]
    lm = [mp.mpf(float(l)) for l in lam]
    am = mp.mpf(alpha)
    bm = 2 * am + 2
    r = [mp.sqrt(sum((xm[c] - Pm[i][c]) ** 2 for c in range(2))) for i in range(3)]
    qm = [lm[i] * r[i] ** bm for i in range(3)]
    gm = [sum(-2 * am * qm[i] * r[i] ** (-bm) * (xm[c] - Pm[i][c]) for i in range(3)) for c in range(2)]
    scale = sum(abs(2 * am * qm[i] * r[i] ** (1 - bm)) for i in range(3))
    worst = max(worst, max(abs(gc) for gc in gm) / scale)
print('     max |grad V_Q(x)(x)| / scale = %s' % mp.nstr(worst, 3))

print('(1b) converse + (1c) type: critical points of random V_q (hybrd from seeds) vs Q(x) and sign(D - c)')
n_pts, worst_q, mism, n_min, n_sad = 0, 0.0, 0, 0, 0
for k in range(120):
    alpha = [0.5, 0.5, 0.25, 1.0, 2.5][k % 5]
    beta = 2 * alpha + 2
    c = (beta - 1) / beta ** 2
    P = rand_triangle_2d(rng, 20.0)
    if k % 2 == 0:
        q = np.exp(rng.normal(scale=0.3, size=3))          # near balanced
    else:                                                 # charges for which a chosen point of M is critical
        q = np.exp(rng.normal(scale=1.0, size=3))
    for x in critical_points_2d(P, q, alpha):
        lam, S, Dcb, Z = S_D_Z(P, x)
        r = np.array([np.linalg.norm(x - p) for p in P])
        Q = lam * r ** beta
        worst_q = max(worst_q, np.max(np.abs(Q / Q.sum() - q / q.sum())))
        ev = np.linalg.eigvalsh(hess_formula(x, P, q, alpha))
        kind = 'min' if ev[0] > 0 else ('saddle' if ev[0] < 0 < ev[1] else 'other')
        n_pts += 1
        n_min += kind == 'min'
        n_sad += kind == 'saddle'
        if (kind == 'min') != (Dcb > c) or kind == 'other':
            mism += 1
print('     %d critical points (%d minima, %d saddles): max |Q(x)/|Q| - q/|q|| = %.2e; type/sign(D-c) mismatches %d'
      % (n_pts, n_min, n_sad, worst_q, mism))

print('(1c) identities for S, D, Z at random points of random triangles')
w = dict(tr=0.0, cb=0.0, z=0.0, det=0.0, M=0)
for k in range(2000):
    alpha = [0.5, 0.25, 1.0, 2.5, 0.1][k % 5]
    beta = 2 * alpha + 2
    c = (beta - 1) / beta ** 2
    P = rand_triangle_2d(rng, 2.0)
    x = rng.dirichlet(np.ones(3)) @ P
    lam, S, Dcb, Z = S_D_Z(P, x)
    w['tr'] = max(w['tr'], abs(np.trace(S) - 1))
    w['cb'] = max(w['cb'], abs(np.linalg.det(S) - Dcb))
    w['z'] = max(w['z'], abs(Dcb - (1 - abs(Z) ** 2) / 4))
    w['det'] = max(w['det'], abs(np.linalg.det(np.eye(2) - beta * S) - beta ** 2 * (Dcb - c)))
    if (Dcb > c) != (abs(Z) < alpha / (alpha + 1)) and abs(Dcb - c) > 1e-12:
        w['M'] += 1
print('     max |tr S - 1| %.1e; max |det S - CauchyBinet| %.1e; max |D - (1-|Z|^2)/4| %.1e; '
      'max |det(I - beta S) - beta^2 (D-c)| %.1e; {D>c} vs {|Z|<a/(a+1)} disagreements %d'
      % (w['tr'], w['cb'], w['z'], w['det'], w['M']))
# the step "both eigenvalues > 1/beta <=> s1 s2 > c" needs 1/beta <= 1/2, i.e. beta >= 2: true for every alpha > 0.

# ---------------------------------------------------------------- (1d)(i) 3D: critical points lie in the plane, inside T
print('(1d-i) 3D Newton (hybr) from off-plane seeds: where do critical points of V in R^3 lie?')
off, outside, inplane, normal_ok, nfound = 0, 0, 0, 0, 0
for k in range(30):
    alpha = [0.5, 0.25, 1.0][k % 3]
    beta = 2 * alpha + 2
    P2 = rand_triangle_2d(rng, 15.0)
    P3 = np.hstack([P2, np.zeros((3, 1))])
    q = np.exp(rng.normal(scale=0.4, size=3))

    def G3(x):
        r = np.array([np.linalg.norm(x - p) for p in P3])
        w = q * r ** (-beta)
        return (w[:, None] * (x - P3)).sum(0) / (w * r).sum()

    pts = []
    for s in range(150):
        s0 = rng.normal(size=3) * 1.5 + np.r_[P2.mean(0), 0.0]
        sol = root(G3, s0, method='hybr', options={'xtol': 1e-14})
        if np.all(np.isfinite(sol.x)) and np.linalg.norm(G3(sol.x)) < 1e-12:
            if all(np.linalg.norm(sol.x - y) > 1e-7 for y in pts):
                pts.append(sol.x)
    for x in pts:
        nfound += 1
        if abs(x[2]) > 1e-9:
            off += 1
            continue
        inplane += 1
        if bary2(P2, x[:2]).min() <= 0:
            outside += 1
        H = hess_formula(x, P3, q, alpha)
        kappa = sum(qi * np.linalg.norm(x - p) ** (-beta) for p, qi in zip(P3, q))
        if abs(H[2, 2] + 2 * alpha * kappa) < 1e-9 * abs(H[2, 2]) and abs(H[0, 2]) + abs(H[1, 2]) < 1e-9 * abs(H[2, 2]):
            normal_ok += 1
print('     %d distinct critical points found: %d off the plane, %d in the plane, %d of these outside T; '
      'normal Hessian eigenvalue = -2 alpha kappa (block diagonal) at %d/%d' % (nfound, off, inplane, outside, normal_ok, inplane))

# ---------------------------------------------------------------- (1d)(ii) subharmonicity
print('(1d-ii) planar Laplacian of V vs 4 alpha^2 sum q_i r_i^(-2alpha-2)')
worst = 0.0
for k in range(50):
    alpha = [0.5, 0.25, 1.0, 2.5, 0.1][k % 5]
    P = rand_triangle_2d(rng)
    q = np.exp(rng.normal(size=3))
    x = rng.normal(size=2)
    if min(np.linalg.norm(x - p) for p in P) < 0.2:
        continue
    _, H_fd = fd_grad_hess_mp(x, P, q, alpha)
    lap = np.trace(H_fd)
    pred = 4 * alpha ** 2 * sum(qi * np.linalg.norm(x - p) ** (-2 * alpha - 2) for p, qi in zip(P, q))
    worst = max(worst, abs(lap - pred) / pred)
print('     max relative error %.2e (the formula is positive, so no planar local maximum)' % worst)

# ---------------------------------------------------------------- (1d)(iii) Poincare-Hopf bookkeeping
print('(1d-iii) index bookkeeping: winding(big circle) - sum winding(small circles) vs found critical points')
bad, tally = 0, {}
for k in range(60):
    alpha = [0.5, 0.5, 0.25, 1.0][k % 4]
    beta = 2 * alpha + 2
    P = rand_triangle_2d(rng, 20.0)
    q = np.exp(rng.normal(scale=0.25 if k % 2 else 1.0, size=3))
    gf = lambda y: grad_formula(y, P, q, alpha)
    cen = P.mean(0)
    Rbig = 50 * np.max(np.linalg.norm(P - cen, axis=1))
    t = np.linspace(0, 2 * math.pi, 4000, endpoint=False)
    circ = np.c_[np.cos(t), np.sin(t)]
    wb = winding(gf, cen + Rbig * circ)
    side = min(np.linalg.norm(P[i] - P[j]) for i in range(3) for j in range(i))
    ws = [winding(gf, p + 1e-3 * side * circ) for p in P]
    total = wb - sum(ws)
    pts = critical_points_2d(P, q, alpha)
    idx = [int(np.sign(np.linalg.det(hess_formula(x, P, q, alpha)))) for x in pts]
    m, s = idx.count(1), idx.count(-1)
    tally[(m, s)] = tally.get((m, s), 0) + 1
    if round(total) != -2 or abs(total - round(total)) > 1e-6 or s - m != 2:
        bad += 1
        print('     mismatch: windings big %.3f small %s; minima %d saddles %d' % (wb, np.round(ws, 3), m, s))
print('     60 configurations: (minima, saddles) tally %s; mismatches %d (big circle winding +1, each charge +1 => '
      'total index -2 => #saddles - #minima = 2, #equilibria = 2 m0 + 2)' % (dict(sorted(tally.items())), bad))

# ---------------------------------------------------------------- (1d)(iv) same bookkeeping with a forced minimum
print('(1d-iv) configurations with a planar minimum forced: q = Q(x0), x0 in M (near-equilateral triangles)')
bad, tally, n = 0, {}, 0
while n < 60:
    alpha = [0.5, 0.5, 0.25, 1.0][n % 4]
    beta = 2 * alpha + 2
    c = (beta - 1) / beta ** 2
    A = np.radians(rng.uniform(48, 72, 2))
    C = math.pi - A.sum()
    if not np.radians(48) <= C <= np.radians(72):
        continue
    L13 = math.sin(A[1]) / math.sin(C)
    P = np.array([[0.0, 0.0], [1.0, 0.0], [L13 * math.cos(A[0]), L13 * math.sin(A[0])]])
    x0 = None
    for t in range(4000):
        y = rng.dirichlet(np.ones(3)) @ P
        if S_D_Z(P, y)[2] > c:
            x0 = y
            break
    if x0 is None:
        continue
    n += 1
    lam0 = bary2(P, x0)
    q = lam0 * np.array([np.linalg.norm(x0 - p) for p in P]) ** beta
    gf = lambda y: grad_formula(y, P, q, alpha)
    cen = P.mean(0)
    t = np.linspace(0, 2 * math.pi, 4000, endpoint=False)
    circ = np.c_[np.cos(t), np.sin(t)]
    wb = winding(gf, cen + 50 * circ)
    ws = [winding(gf, p + 1e-3 * circ) for p in P]
    pts = critical_points_2d(P, q, alpha)
    idx = [int(np.sign(np.linalg.det(hess_formula(x, P, q, alpha)))) for x in pts]
    m, s = idx.count(1), idx.count(-1)
    has_x0 = any(np.linalg.norm(x - x0) < 1e-8 for x in pts)
    tally[(m, s)] = tally.get((m, s), 0) + 1
    if round(wb - sum(ws)) != -2 or s - m != 2 or not has_x0:
        bad += 1
        print('     mismatch: m %d s %d, x0 found %s' % (m, s, has_x0))
print('     60 configurations: (minima, saddles) tally %s; mismatches %d' % (dict(sorted(tally.items())), bad))
