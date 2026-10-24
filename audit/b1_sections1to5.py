"""AUDIT_2 item 1: MIXED.md sections 1-5, independent checks (Coulomb unless stated; 40-digit mpmath via b0_common).

 (1)  planarity: 3D root finding from off-plane seeds for random mixed charges never finds an off-plane critical point;
 (2)  critical points of random mixed V_q (found by my own 2D finder) satisfy w = kappa lambda, lie in R+ (kappa > 0) or
      R- (kappa < 0), and Hess V = 2 alpha kappa (beta S - I) there;
 (3)  Lemma D: D = P_x sin s1 sin s2 sin s3 / (2 Delta) at random points (all regions off the side lines);
 (4)  R-: prod sin sigma > 0, P_x < 0, D < 0;  R+: prod sin sigma < 0;  R+out: t_i in (0, A_i), s < A1 + A2;
 (5)  Lemma P (signed): r_j r_k sin(sigma_i - A_i) = P_x sin A_i at random points of all regions;
 (6)  injectivity / image of t on R+out: the geometric inverse reproduces x from t(x); random t in Q map into R+out with
      the same t (so t(R+out) = Q);
 (7)  separated charges lambda_i r_i^beta = C(x) sin sigma_i |sin(sigma_i - A_i)|^b / sin^b A_i, and
      grad H(t(x)) = psi_q for q = Q(x);
 (8)  W = (-U1, -U2, U3), E = e2(W) = U1U2 - U3(U1+U2) = 1/D; e2(x~) = e2(y~) = 1, E = 2 + 2<x~,y~>, E >= 4 on Q with
      equality at t0 = (A1/2, A2/2) = t(excentre opposite p3);
 (9)  d = x~ - b y~, e2(d) = (b+1)^2 - b E; Hess H (finite differences of H) = Gram matrix [[d1+d3,d3],[d3,d2+d3]];
      Hess H negative definite at t0;
 (10) index count: winding numbers on small circles (=1 for either sign) and on a large circle (=1 if Q_tot != 0,
      =2 if Q_tot = 0), and #saddles - #extrema = 2 (resp. 1) for the critical points found.
usage: python b1_sections1to5.py"""
import math

import mpmath as mp
import numpy as np
from scipy.optimize import root

import b0_common as C

rng = np.random.default_rng(20261024)
alpha = mp.mpf('0.5')
b = C.b_of_alpha(alpha)
beta_exp = b + 1
c = b / (b + 1) ** 2
worst = {}


def upd(k, v):
    worst[k] = max(worst.get(k, mp.mpf(0)), abs(v))


def cnt(k, flag):
    worst[k] = worst.get(k, 0) + (1 if flag else 0)


def rand_shape(min_deg=1.0):
    while True:
        w = rng.dirichlet(np.ones(3)) * 180
        if w.min() > min_deg:
            return [mp.mpf(float(v)) * mp.pi / 180 for v in w[:2]]


def pts_complex(P):
    return np.array([[float(p.real), float(p.imag)] for p in P])


# ---------------------------------------------------------------- (1) planarity in R^3
def crit3d(P3, q, al, nseed=200):
    beta_ = 2 * al + 2

    def G(x):
        r = np.linalg.norm(x[None] - P3, axis=1)
        w = q * r ** (-beta_)
        return (w[:, None] * (x[None] - P3)).sum(0) / (np.abs(w) * r).sum()

    out = []
    for _ in range(nseed):
        x0 = P3.mean(0) + rng.normal(size=3) * rng.choice([0.5, 2, 6])
        sol = root(G, x0, method='hybr', options={'xtol': 1e-14})
        x = sol.x
        if np.all(np.isfinite(x)) and np.linalg.norm(x) < 1e4 and np.linalg.norm(G(x)) < 1e-12:
            if min(np.linalg.norm(x - p) for p in P3) > 1e-9 and all(np.linalg.norm(x - y) > 1e-7 for y in out):
                out.append(x)
    return out


off, tot = 0, 0
for k in range(40):
    A1, A2 = rand_shape(5)
    A, P = C.shape(A1, A2)
    P2 = pts_complex(P)
    P3 = np.c_[P2, np.zeros(3)]
    q = np.exp(rng.normal(size=3)); q[2] = -q[2]
    for x in crit3d(P3, q, 0.5):
        tot += 1
        off += abs(x[2]) > 1e-9
cnt('(1) off-plane critical points found (of %s)' % 'all', off)
print('(1) planarity: %d critical points found by 3D root finding (40 configurations), %d off the plane' % (tot, off))


# ---------------------------------------------------------------- 2D critical point finder (auditor's own)
def crit2d(P2, q, al=0.5, Rbig=60.0):
    beta_ = 2 * al + 2

    def G(x):
        r = np.linalg.norm(x[None] - P2, axis=1)
        w = q * r ** (-beta_)
        return (w[:, None] * (x[None] - P2)).sum(0) / (np.abs(w) * r).sum()

    cen = P2.mean(0)
    scale = np.max(np.linalg.norm(P2 - cen, axis=1))
    seeds = []
    for rr in np.geomspace(0.05, Rbig, 36) * scale:
        for th in np.linspace(0, 2 * math.pi, 48, endpoint=False):
            seeds.append(cen + rr * np.array([math.cos(th), math.sin(th)]))
    for v in range(3):
        for rr in np.geomspace(1e-5, 0.5, 14) * scale:
            for th in np.linspace(0, 2 * math.pi, 16, endpoint=False):
                seeds.append(P2[v] + rr * np.array([math.cos(th), math.sin(th)]))
    out = []
    for s0 in seeds:
        sol = root(G, s0, method='hybr', options={'xtol': 1e-14})
        x = sol.x
        if not np.all(np.isfinite(x)) or np.linalg.norm(x - cen) > 1e3 * scale:
            continue
        if np.linalg.norm(G(x)) > 1e-12 or min(np.linalg.norm(x - p) for p in P2) < 1e-10:
            continue
        if all(np.linalg.norm(x - y) > 1e-8 * scale for y in out):
            out.append(x)
    return out


def hessV(P2, q, x, al=0.5):
    beta_ = 2 * al + 2
    H = np.zeros((2, 2))
    for p, qi in zip(P2, q):
        d = x - p
        r = np.linalg.norm(d)
        u = d / r
        H += -2 * al * qi * r ** (-beta_) * (np.eye(2) - beta_ * np.outer(u, u))
    return H


# ---------------------------------------------------------------- (2) w = kappa lambda, regions, Hessian formula
nR = {'R+': 0, 'R-': 0, 'other': 0}
for k in range(60):
    A1, A2 = rand_shape(5)
    A, P = C.shape(A1, A2)
    P2 = pts_complex(P)
    q = np.exp(rng.normal(size=3) * 0.8); q[2] = -q[2] * rng.uniform(0.2, 3)
    for x in crit2d(P2, q):
        xm = mp.mpc(x[0], x[1])
        g = C.geom(P, xm)
        lam = np.array([float(v) for v in g['lam']])
        r = np.array([float(v) for v in g['r']])
        w = q * r ** (-3.0)
        kap = w / lam
        upd('(2) max rel spread of w_i/lambda_i', (kap.max() - kap.min()) / abs(kap).max())
        reg = 'R+' if (lam[0] > 0 and lam[1] > 0 and lam[2] < 0) else ('R-' if (lam[0] < 0 and lam[1] < 0 and lam[2] > 0) else 'other')
        nR[reg] += 1
        cnt('(2) kappa sign inconsistent with region', (reg == 'R+' and kap[0] < 0) or (reg == 'R-' and kap[0] > 0))
        S = np.array([[float(g['S'][i, j]) for j in range(2)] for i in range(2)])
        Hf = 2 * 0.5 * kap.mean() * (3 * S - np.eye(2))
        upd('(2) rel |Hess V - 2 alpha kappa (beta S - I)|', np.max(np.abs(hessV(P2, q, x) - Hf)) / np.max(np.abs(Hf)))
print('(2) critical points by region: %s' % nR)

# ---------------------------------------------------------------- (3)-(5), (8) on random points of the plane
for k in range(3000):
    A1, A2 = rand_shape(1)
    A, P = C.shape(A1, A2)
    cen = sum(P) / 3
    x = cen + mp.mpc(*rng.normal(size=2)) * rng.choice([0.3, 1, 3, 10])
    g = C.geom(P, x)
    if min(abs(v) for v in g['lam']) < mp.mpf('1e-6'):
        continue
    sig, lam, r, Px, D = g['sig'], g['lam'], g['r'], g['Px'], g['D']
    prodsin = mp.sin(sig[0]) * mp.sin(sig[1]) * mp.sin(sig[2])
    upd('(3) rel |D - P_x prod sin sigma / (2 Delta)|', (D - Px * prodsin / g['D2']) / max(abs(D), mp.mpf('1e-30')))
    for i in range(3):
        j, kk = (i + 1) % 3, (i + 2) % 3
        upd('(5) rel |r_j r_k sin(sigma_i - A_i) - P_x sin A_i|',
            (r[j] * r[kk] * mp.sin(sig[i] - A[i]) - Px * mp.sin(A[i])) / (abs(Px) * mp.sin(A[i]) + r[j] * r[kk] * mp.mpf('1e-30')))
    inRm = lam[0] < 0 and lam[1] < 0 and lam[2] > 0
    inRp = lam[0] > 0 and lam[1] > 0 and lam[2] < 0
    if inRm:
        cnt('(4) R-: prod sin <= 0 or P_x >= 0 or D >= 0', not (prodsin > 0 and Px < 0 and D < 0))
    if inRp:
        cnt('(4) R+: prod sin >= 0', not (prodsin < 0))
        if Px < 0:
            t1, t2, s = sig[0], sig[1], -sig[2]
            cnt('(4) R+out: t_i not in (0, A_i) or s != t1 + t2 or s >= A1 + A2',
                not (0 < t1 < A[0] and 0 < t2 < A[1] and abs(s - t1 - t2) < mp.mpf('1e-30') and s < A[0] + A[1]))
            # (6) geometric inverse
            xr, (r1, r2) = C.x_of_t(P, t1, t2)
            upd('(6) |x_of_t(t(x)) - x| / |x|', abs(xr - x) / abs(x))
            # (7) separated charges
            Cx = (r[0] * r[1] * r[2]) ** (b + 1) / (g['D2'] * abs(Px) ** b)
            for i in range(3):
                val = lam[i] * r[i] ** beta_exp
                pred = Cx * mp.sin(sig[i]) * abs(mp.sin(sig[i] - A[i])) ** b / mp.sin(A[i]) ** b
                upd('(7) rel |lambda_i r_i^beta - C sin sigma |sin(sigma - A)|^b / sin^b A|', val / pred - 1)
            Q = [lam[i] * r[i] ** beta_exp for i in range(3)]
            psi_q = (mp.log(Q[0] / abs(Q[2])) + b * mp.log(mp.sin(A[0]) / mp.sin(A[2])),
                     mp.log(Q[1] / abs(Q[2])) + b * mp.log(mp.sin(A[1]) / mp.sin(A[2])))
            gH = C.gradH(A, t1, t2, b)
            upd('(7) |grad H(t(x)) - psi_q|', max(abs(gH[0] - psi_q[0]), abs(gH[1] - psi_q[1])))
            # (8)
            U1, U2, U3 = C.U(t1, A[0]), C.U(t2, A[1]), C.U(s, A[0] + A[1])
            W = [mp.cot(sig[i] - A[i]) - mp.cot(sig[i]) for i in range(3)]
            upd('(8) |W - (-U1, -U2, U3)|', max(abs(W[0] + U1), abs(W[1] + U2), abs(W[2] - U3)) / max(U1, U2, U3))
            E = U1 * U2 - U3 * (U1 + U2)
            upd('(8) |E D - 1|', E * D - 1)
print('(3)-(8) evaluated at random points of the plane (3000 tries), see table below')

# random t in Q -> x in R+out with the same t; E >= 4 on Q; hyperbolic identities; (9)
for k in range(1500):
    A1, A2 = rand_shape(0.5)
    A, P = C.shape(A1, A2)
    t1 = A[0] * mp.mpf(float(rng.uniform(0.001, 0.999)))
    t2 = A[1] * mp.mpf(float(rng.uniform(0.001, 0.999)))
    s = t1 + t2
    x, (r1, r2) = C.x_of_t(P, t1, t2)
    cnt('(6) x_of_t: r1 or r2 <= 0 (wrong branch)', not (r1 > 0 and r2 > 0))
    g = C.geom(P, x)
    lam = g['lam']
    cnt('(6) random t in Q maps outside R+out', not (lam[0] > 0 and lam[1] > 0 and lam[2] < 0 and g['Px'] < 0))
    upd('(6) |t(x_of_t(t)) - t|', max(abs(g['sig'][0] - t1), abs(g['sig'][1] - t2)))
    xt = (mp.cot(t1), mp.cot(t2), -mp.cot(s))
    yt = (mp.cot(A[0] - t1), mp.cot(A[1] - t2), -mp.cot(A[0] + A[1] - s))
    e2 = lambda v: v[0] * v[1] + v[1] * v[2] + v[2] * v[0]
    upd('(8) |e2(x~) - 1|, |e2(y~) - 1|', max(abs(e2(xt) - 1), abs(e2(yt) - 1)))
    E = C.E_of_t(A, t1, t2)
    ip = (xt[0] * (yt[1] + yt[2]) + xt[1] * (yt[0] + yt[2]) + xt[2] * (yt[0] + yt[1])) / 2
    upd('(8) |E - (2 + 2<x~,y~>)|', (E - (2 + 2 * ip)) / E)
    cnt('(8) E < 4 somewhere in Q', E < 4 - mp.mpf('1e-25'))
    d = (xt[0] - b * yt[0], xt[1] - b * yt[1], xt[2] - b * yt[2])
    phi1p = mp.cot(t1) - b * mp.cot(A[0] - t1)
    phi2p = mp.cot(t2) - b * mp.cot(A[1] - t2)
    phi3p = mp.cot(s) + b * mp.cot(s + A[2])
    upd('(9) |d - (phi1\', phi2\', -phi3\')|', max(abs(d[0] - phi1p), abs(d[1] - phi2p), abs(d[2] + phi3p)))
    upd('(9) rel |e2(d) - ((b+1)^2 - bE)|', (e2(d) - ((b + 1) ** 2 - b * E)) / ((b + 1) ** 2 + b * E))
    if k % 10 == 0:
        Hf = lambda u, v: None
        f = lambda u, v: (lambda gg: None)
        # Hessian of H via grad H (finite differences of the explicit gradient)
        g1 = lambda u, v: C.gradH(A, u, v, b)[0]
        g2 = lambda u, v: C.gradH(A, u, v, b)[1]
        Hm = mp.matrix([[mp.diff(g1, (t1, t2), (1, 0)), mp.diff(g1, (t1, t2), (0, 1))],
                        [mp.diff(g2, (t1, t2), (1, 0)), mp.diff(g2, (t1, t2), (0, 1))]])
        Gm = mp.matrix([[d[0] + d[2], d[2]], [d[2], d[1] + d[2]]])
        upd('(9) |Hess H - Gram(d)| (rel)', max(abs(Hm[i, j] - Gm[i, j]) for i in range(2) for j in range(2)) / max(abs(Gm[i, j]) for i in range(2) for j in range(2)))
        upd('(9) |det Gram - e2(d)| (rel)', (mp.det(Gm) - e2(d)) / (abs(e2(d)) + 1))
# t0 and the excentre; sign of Hess H at t0
for k in range(200):
    A1, A2 = rand_shape(0.5)
    A, P = C.shape(A1, A2)
    a = [abs(P[(i + 1) % 3] - P[(i + 2) % 3]) for i in range(3)]
    J3 = (a[0] * P[0] + a[1] * P[1] - a[2] * P[2]) / (a[0] + a[1] - a[2])     # excentre opposite p3
    g = C.geom(P, J3)
    upd('(8) |t(excentre_3) - (A1/2, A2/2)|', max(abs(g['sig'][0] - A[0] / 2), abs(g['sig'][1] - A[1] / 2)))
    upd('(8) |E(t0) - 4|', C.E_of_t(A, A[0] / 2, A[1] / 2) - 4)
    t1, t2 = A[0] / 2, A[1] / 2
    g1 = lambda u, v: C.gradH(A, u, v, b)[0]
    g2 = lambda u, v: C.gradH(A, u, v, b)[1]
    Hm = mp.matrix([[mp.diff(g1, (t1, t2), (1, 0)), mp.diff(g1, (t1, t2), (0, 1))],
                    [mp.diff(g2, (t1, t2), (1, 0)), mp.diff(g2, (t1, t2), (0, 1))]])
    ev = mp.eig(Hm)[0]
    cnt('(9) Hess H at t0 not negative definite', not all(mp.re(v) < 0 for v in ev))

# ---------------------------------------------------------------- (10) index count incl. Q_tot = 0
def winding(fun, curve):
    v = np.array([fun(c_) for c_ in curve])
    a = np.arctan2(v[:, 1], v[:, 0])
    da = np.diff(np.concatenate([a, a[:1]]))
    da = (da + math.pi) % (2 * math.pi) - math.pi
    return da.sum() / (2 * math.pi)


tally = {}
for k in range(80):
    A1, A2 = rand_shape(8)
    A, P = C.shape(A1, A2)
    P2 = pts_complex(P)
    q = np.exp(rng.normal(size=3) * 0.6)
    dipole = k % 2 == 0
    q[2] = -(q[0] + q[1]) if dipole else -q[2] * rng.uniform(0.3, 3)
    gradV = lambda y: sum(-qi * np.linalg.norm(y - p) ** -3 * (y - p) for p, qi in zip(P2, q))
    cen = P2.mean(0)
    tt = np.linspace(0, 2 * math.pi, 20000, endpoint=False)
    circ = np.c_[np.cos(tt), np.sin(tt)]
    w_inf = winding(gradV, cen + 400 * circ)
    w_ch = [winding(gradV, p + 1e-4 * circ) for p in P2]
    pts = crit2d(P2, q)
    idx = [int(np.sign(np.linalg.det(hessV(P2, q, x)))) for x in pts]
    ext, sad = idx.count(1), idx.count(-1)
    key = ('Qtot=0' if dipole else 'Qtot!=0', round(w_inf), tuple(round(v) for v in w_ch), ext, sad)
    tally[key] = tally.get(key, 0) + 1
print('(10) index count: (case, winding at infinity, windings at charges, #extrema, #saddles): counts')
for key in sorted(tally):
    print('     %s: %d' % (key, tally[key]))

print('\nTABLE (max abs/rel errors, and counts of failures):')
for key in sorted(worst):
    v = worst[key]
    print('  %-80s %s' % (key, mp.nstr(v, 3) if not isinstance(v, int) else v))
