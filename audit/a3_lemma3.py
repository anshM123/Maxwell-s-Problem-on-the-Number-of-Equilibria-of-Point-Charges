"""AUDIT item 3: Lemma 3 (Psi is a gradient in the angle coordinates theta = (theta_1, theta_2)).

All angles here are computed as UNSIGNED angles theta_i = acos(<p_j - x, p_k - x> / (r_j r_k)) in (0, pi), i.e. exactly
the definition in PROOF.md, independent of any orientation convention.  30-digit mpmath throughout.
 (3a) orientation: theta_i = phi_k - phi_j (mod 2 pi) for cyclic (i,j,k) when p1 p2 p3 is counterclockwise
      (phi_a = arg(x - p_a)); for clockwise labelling it is phi_j - phi_k.  Also grad theta_i = J(e_k - e_j)
      (ccw), J = rotation by +90 deg, e_a = (x - p_a)/r_a^2.
 (3b) d Re f ^ d Im g = Re(f' conj(g')) dx^dy for holomorphic f, g (f = log(z - p_a), g = log(z - p_b), and a generic
      pair f = z^3 + e^z, g = sin(z) + z^2).
 (3c) omega = sum_i log r_i dtheta_i is closed: curl computed by numerical differentiation (normalised), and the loop
      integral of omega around random circles inside T (trapezoid rule, exponentially convergent) is 0.  For contrast
      the same loop integrals for the non-closed form sum_i r_i dtheta_i are reported.
 (3d) sum_i log sin theta_i dtheta_i = d(-sum_i Lob(theta_i)), Lob(t) = Cl_2(2t)/2 (Lobachevsky), checked by comparing
      the x-gradient of -sum Lob(theta_i) with sum log sin theta_i grad theta_i.
 (3e) D_theta Psi = DPsi A^{-1} (A = D_x theta) is symmetric: random shapes (both orientations), random points
      (including points near edges and vertices), random beta.
 (3f) explicit form (not used in the proof): d(-sum over the nine angles of the triangles x p_j p_k of Lob) - omega
      = sum_i c_i dtheta_i with constant c_i?  (residual of a constant fit over many points)
usage: python a3_lemma3.py"""
import math

import mpmath as mp
import numpy as np

mp.mp.dps = 30
rng = np.random.default_rng(11)
J90 = lambda v: (-v[1], v[0])


def rand_tri(rng, min_deg=3.0, orient=None):
    while True:
        P = rng.normal(size=(3, 2))
        cr = (P[1, 0] - P[0, 0]) * (P[2, 1] - P[0, 1]) - (P[2, 0] - P[0, 0]) * (P[1, 1] - P[0, 1])
        if orient is not None and np.sign(cr) != orient:
            P = P[[0, 2, 1]]
        ang = []
        for i in range(3):
            a, b = P[(i + 1) % 3] - P[i], P[(i + 2) % 3] - P[i]
            ang.append(math.degrees(math.acos(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b))))
        if min(ang) > min_deg:
            return [[mp.mpf(float(v)) for v in p] for p in P]


def theta(i, x, y, P):
    j, k = (i + 1) % 3, (i + 2) % 3
    ux, uy = P[j][0] - x, P[j][1] - y
    vx, vy = P[k][0] - x, P[k][1] - y
    return mp.acos((ux * vx + uy * vy) / mp.sqrt((ux * ux + uy * uy) * (vx * vx + vy * vy)))


def rdist(i, x, y, P):
    return mp.sqrt((x - P[i][0]) ** 2 + (y - P[i][1]) ** 2)


def signed2(P):
    return (P[1][0] - P[0][0]) * (P[2][1] - P[0][1]) - (P[2][0] - P[0][0]) * (P[1][1] - P[0][1])


def lam(i, x, y, P):
    j, k = (i + 1) % 3, (i + 2) % 3
    return ((P[j][0] - x) * (P[k][1] - y) - (P[j][1] - y) * (P[k][0] - x)) / signed2(P)


def rand_point(P, rng, near=None):
    if near is None:
        w = rng.dirichlet(np.ones(3))
    elif near == 'edge':
        w = rng.dirichlet(np.ones(3)); w[rng.integers(3)] = 1e-4; w /= w.sum()
    else:
        w = np.full(3, 1e-4); w[rng.integers(3)] = 1.0; w = rng.dirichlet(w * 50 + 1e-3); w = np.maximum(w, 1e-5); w /= w.sum()
    x = sum(mp.mpf(float(w[i])) * P[i][0] for i in range(3))
    y = sum(mp.mpf(float(w[i])) * P[i][1] for i in range(3))
    return x, y


def grad(f, x, y):
    return (mp.diff(lambda u: f(u, y), x), mp.diff(lambda v: f(x, v), y))


# ------------------------------------------------------------------ (3a)
bad_ccw, bad_cw, worst_g = 0, 0, mp.mpf(0)
for k in range(200):
    orient = 1 if k % 2 == 0 else -1
    P = rand_tri(rng, orient=orient)
    x, y = rand_point(P, rng)
    phi = [mp.atan2(y - P[a][1], x - P[a][0]) for a in range(3)]
    for i in range(3):
        j, k2 = (i + 1) % 3, (i + 2) % 3
        th = theta(i, x, y, P)
        cand = (phi[k2] - phi[j]) if orient == 1 else (phi[j] - phi[k2])
        if abs(mp.fmod(cand - th + 10 * mp.pi, 2 * mp.pi)) > mp.mpf('1e-20') and \
                abs(mp.fmod(cand - th + 10 * mp.pi, 2 * mp.pi) - 2 * mp.pi) > mp.mpf('1e-20'):
            if orient == 1:
                bad_ccw += 1
            else:
                bad_cw += 1
        if orient == 1:
            g = grad(lambda u, v: theta(i, u, v, P), x, y)
            e = [((x - P[a][0]) / rdist(a, x, y, P) ** 2, (y - P[a][1]) / rdist(a, x, y, P) ** 2) for a in range(3)]
            pred = J90((e[k2][0] - e[j][0], e[k2][1] - e[j][1]))
            worst_g = max(worst_g, max(abs(g[0] - pred[0]), abs(g[1] - pred[1])) / max(abs(g[0]), abs(g[1])))
print('(3a) theta_i = phi_k - phi_j mod 2pi (ccw): failures %d;  theta_i = phi_j - phi_k (cw): failures %d;  '
      'grad theta_i = J(e_k - e_j) (ccw): max rel err %s' % (bad_ccw, bad_cw, mp.nstr(worst_g, 3)))

# ------------------------------------------------------------------ (3b)
worst = mp.mpf(0)
for k in range(60):
    pa = mp.mpc(*rng.normal(size=2)); pb = mp.mpc(*rng.normal(size=2))
    if k % 2 == 0:
        f = lambda z: mp.log(z - pa); fp = lambda z: 1 / (z - pa)
        g = lambda z: mp.log(z - pb); gp = lambda z: 1 / (z - pb)
    else:
        f = lambda z: z ** 3 + mp.exp(z); fp = lambda z: 3 * z ** 2 + mp.exp(z)
        g = lambda z: mp.sin(z) + z ** 2; gp = lambda z: mp.cos(z) + 2 * z
    z0 = mp.mpc(*rng.normal(size=2))
    if min(abs(z0 - pa), abs(z0 - pb)) < 0.3:
        continue
    # make Im g continuous near z0 for the log case: log(z - pb) - (branch value at z0) is fine locally
    ReF = lambda u, v: mp.re(f(mp.mpc(u, v)))
    ImG = lambda u, v: mp.im(g(mp.mpc(u, v)) / (mp.mpc(u, v) - pb) * (mp.mpc(u, v) - pb)) if False else mp.im(g(mp.mpc(u, v)))
    u0, v0 = mp.re(z0), mp.im(z0)
    gu = grad(ReF, u0, v0)
    gw = grad(ImG, u0, v0)
    wedge = gu[0] * gw[1] - gu[1] * gw[0]
    pred = mp.re(fp(z0) * mp.conj(gp(z0)))
    worst = max(worst, abs(wedge - pred) / abs(pred))
print('(3b) d Re f ^ d Im g = Re(f\' conj g\') dx^dy: max rel err %s' % mp.nstr(worst, 3))

# ------------------------------------------------------------------ (3c) closedness of omega
worst_curl = mp.mpf(0)
for k in range(40):
    P = rand_tri(rng, orient=1 if k % 2 else -1)
    x0, y0 = rand_point(P, rng, None if k % 4 < 2 else 'edge')

    def om(u, v, comp):
        s = mp.mpf(0)
        for i in range(3):
            gth = grad(lambda a, b: theta(i, a, b, P), u, v)
            s += mp.log(rdist(i, u, v, P)) * gth[comp]
        return s

    d_omy_dx = mp.diff(lambda u: om(u, y0, 1), x0)
    d_omx_dy = mp.diff(lambda v: om(x0, v, 0), y0)
    scale = abs(d_omy_dx) + abs(d_omx_dy) + mp.mpf('1e-30')
    worst_curl = max(worst_curl, abs(d_omy_dx - d_omx_dy) / scale)
print('(3c) curl of omega = sum log r_i dtheta_i (numerical differentiation, normalised): max %s' % mp.nstr(worst_curl, 3))


def loop_integral(Pf, center, rad, kind, N=3000):
    """trapezoid rule (float64) for the line integral of sum_i F(r_i) dtheta_i around a circle, with
    grad theta_i = J(e_k - e_j) for ccw labelling (verified in 3a)"""
    t = np.linspace(0, 2 * math.pi, N, endpoint=False)
    X = center[None, :] + rad * np.c_[np.cos(t), np.sin(t)]
    dX = rad * np.c_[-np.sin(t), np.cos(t)]
    tot = np.zeros(N)
    e = [(X - p) / np.sum((X - p) ** 2, axis=1)[:, None] for p in Pf]
    r = [np.sqrt(np.sum((X - p) ** 2, axis=1)) for p in Pf]
    for i in range(3):
        j, k2 = (i + 1) % 3, (i + 2) % 3
        gth = np.c_[-(e[k2][:, 1] - e[j][:, 1]), e[k2][:, 0] - e[j][:, 0]]
        F = np.log(r[i]) if kind == 'log' else (r[i] if kind == 'r' else np.log(np.sin(np.arccos(np.clip(
            np.sum((Pf[j] - X) * (Pf[k2] - X), axis=1) / r[j] / r[k2], -1, 1)))))
        tot += F * np.sum(gth * dX, axis=1)
    return tot.mean() * 2 * math.pi


worst_loop, contrast, worst_loop_sin = 0.0, 0.0, 0.0
for k in range(200):
    P = rand_tri(rng, orient=1)
    Pf = np.array([[float(v) for v in p] for p in P])
    w = rng.dirichlet(np.ones(3))
    cen = w @ Pf
    # largest circle around cen inside T: distance to the three side lines
    dist = []
    for i in range(3):
        a, b = Pf[(i + 1) % 3], Pf[(i + 2) % 3]
        n = np.array([-(b - a)[1], (b - a)[0]]); n /= np.linalg.norm(n)
        dist.append(abs(np.dot(cen - a, n)))
    rad = 0.98 * min(dist)
    worst_loop = max(worst_loop, abs(loop_integral(Pf, cen, rad, 'log')))
    worst_loop_sin = max(worst_loop_sin, abs(loop_integral(Pf, cen, rad, 'sin')))
    contrast = max(contrast, abs(loop_integral(Pf, cen, rad, 'r')))
print('(3c) loop integrals around 200 random maximal circles in T: |oint sum log r_i dtheta_i| max %.2e;  '
      '|oint sum log sin theta_i dtheta_i| max %.2e;  contrast |oint sum r_i dtheta_i| max %.2e'
      % (worst_loop, worst_loop_sin, contrast))

# ------------------------------------------------------------------ (3d) exactness of sum log sin theta_i dtheta_i
lob = lambda t: mp.clsin(2, 2 * t) / 2
worst = mp.mpf(0)
for k in range(40):
    P = rand_tri(rng, orient=1 if k % 2 else -1)
    x0, y0 = rand_point(P, rng)
    g1 = grad(lambda u, v: -sum(lob(theta(i, u, v, P)) for i in range(3)), x0, y0)
    g2 = [mp.mpf(0), mp.mpf(0)]
    for i in range(3):
        gth = grad(lambda a, b: theta(i, a, b, P), x0, y0)
        for c in range(2):
            g2[c] += mp.log(mp.sin(theta(i, x0, y0, P))) * gth[c]
    worst = max(worst, max(abs(g1[c] - g2[c]) for c in range(2)) / max(abs(g2[0]), abs(g2[1])))
print('(3d) grad(-sum Lob(theta_i)) vs sum log sin theta_i grad theta_i: max rel err %s' % mp.nstr(worst, 3))

# ------------------------------------------------------------------ (3e) symmetry of DPsi A^{-1}
worst = mp.mpf(0)
nchk = 0
for k in range(150):
    P = rand_tri(rng, min_deg=1.0, orient=1 if k % 2 else -1)
    x0, y0 = rand_point(P, rng, [None, None, 'edge', 'vertex'][k % 4])
    bet = mp.mpf(2) + mp.mpf(float(rng.uniform(0.02, 8)))

    def L(i, u, v):
        return mp.log(lam(i, u, v, P)) + bet * mp.log(rdist(i, u, v, P))

    try:
        DP = mp.matrix([list(grad(lambda u, v: L(0, u, v) - L(2, u, v), x0, y0)),
                        list(grad(lambda u, v: L(1, u, v) - L(2, u, v), x0, y0))])
        A = mp.matrix([list(grad(lambda u, v: theta(0, u, v, P), x0, y0)),
                       list(grad(lambda u, v: theta(1, u, v, P), x0, y0))])
    except Exception as exc:          # numerical differentiation step leaving T very near a vertex
        continue
    Jt = DP * A ** -1
    nchk += 1
    worst = max(worst, abs(Jt[0, 1] - Jt[1, 0]) / max(abs(Jt[i, j]) for i in range(2) for j in range(2)))
print('(3e) DPsi A^-1 symmetric: %d points (both orientations; a quarter near edges, a quarter near vertices; '
      'beta in (2, 10)): max relative asymmetry %s' % (nchk, mp.nstr(worst, 3)))

# ------------------------------------------------------------------ (3f) explicit form (NOTE only)
mp.mp.dps = 20
P = rand_tri(rng, min_deg=20, orient=1)


def angle_at(o, a, b):
    ux, uy = a[0] - o[0], a[1] - o[1]
    vx, vy = b[0] - o[0], b[1] - o[1]
    return mp.acos((ux * vx + uy * vy) / mp.sqrt((ux * ux + uy * uy) * (vx * vx + vy * vy)))


def nine(u, v):
    s = mp.mpf(0)
    X = (u, v)
    for i in range(3):
        j, k2 = P[(i + 1) % 3], P[(i + 2) % 3]
        s += lob(angle_at(X, j, k2)) + lob(angle_at(j, k2, X)) + lob(angle_at(k2, j, X))
    return -s


rows, rhs = [], []
for k in range(25):
    x0, y0 = rand_point(P, rng)
    gG = grad(nine, x0, y0)
    gth = [grad(lambda a, b: theta(i, a, b, P), x0, y0) for i in range(3)]
    om = [sum(mp.log(rdist(i, x0, y0, P)) * gth[i][c] for i in range(3)) for c in range(2)]
    for c in range(2):                      # gG - om = c1 grad th1 + c2 grad th2 + c3 grad th3 (c3 := 0, sum grad = 0)
        rows.append([float(gth[0][c]), float(gth[1][c])])
        rhs.append(float(gG[c] - om[c]))
coef, *_ = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)
res = np.max(np.abs(np.array(rows) @ coef - np.array(rhs)))
side = [float(mp.sqrt((P[(i + 1) % 3][0] - P[(i + 2) % 3][0]) ** 2 + (P[(i + 1) % 3][1] - P[(i + 2) % 3][1]) ** 2)) for i in range(3)]
print('(3f) d(-sum_9 Lob) - omega = c1 dtheta1 + c2 dtheta2 (c3 = 0): c = %s, max residual %.2e; '
      'compare log(side_3/side_1), log(side_3/side_2) = %s'
      % (np.round(coef, 10).tolist(), res, np.round([math.log(side[2] / side[0]), math.log(side[2] / side[1])], 10).tolist()))
