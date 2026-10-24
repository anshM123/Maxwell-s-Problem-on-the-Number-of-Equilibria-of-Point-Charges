"""AUDIT_2 item 5: MIXED.md section 5 (degeneracy) and section 7 (assembly), numerical side.

 (1) mixed signs: at critical points x in R+ and R- (q = Q(x) with its signs), Hess V_q = 2 alpha kappa L o DPsi with
     kappa = w_i/lambda_i (may be negative), L(e) = sum lambda_i e_i (x - p_i) of rank 2 with kernel R(1,1,1),
     DPsi = derivative of log|Q|;
 (2) collinear mixed charges in R^3 (random positions and sign orders): on-axis equilibria (zeros of V' on the four
     axis intervals), their non-degeneracy, the Hessian structure diag(f'', h, h), and off-axis equilibria (circles,
     solutions of sum w_i = 0, sum w_i z_i = 0 in the (z, rho) half-plane);
 (3) the perturbation step of 7(b): move p3 off the axis by eps; every non-degenerate on-axis equilibrium persists as a
     critical point of the perturbed (non-collinear) configuration, whose total number of critical points is <= 4.
usage: python b5_section7.py"""
import math

import mpmath as mp
import numpy as np
from scipy.optimize import brentq, root

import b0_common as C

rng = np.random.default_rng(57)
mp.mp.dps = 30

# ---------------------------------------------------------------- (1)
worst, worstL = mp.mpf(0), mp.mpf(0)
for k in range(120):
    while True:
        w = rng.dirichlet(np.ones(3)) * 180
        if w.min() > 3:
            break
    A, P = C.shape(*[mp.mpf(float(v)) * mp.pi / 180 for v in w[:2]])
    # random point of R+ (k even) or R- (k odd): x = p3 + a (p1 - p3) + b (p2 - p3), a, b > 0, a + b > 1 (R+) or a, b < 0
    if k % 2 == 0:
        a_, b_ = rng.uniform(0.1, 3, 2)
        if a_ + b_ <= 1.05:
            a_, b_ = a_ + 1, b_ + 1
    else:
        a_, b_ = -rng.uniform(0.05, 3, 2)
    x0 = P[2] + mp.mpf(float(a_)) * (P[0] - P[2]) + mp.mpf(float(b_)) * (P[1] - P[2])
    g = C.geom(P, x0)
    lam, r = g['lam'], g['r']
    Q = [lam[i] * r[i] ** 3 for i in range(3)]
    kappa = Q[0] * r[0] ** -3 / lam[0]                   # = 1 by construction; sign carried by Q
    H = mp.matrix(2, 2)
    for i in range(3):
        d = x0 - P[i]
        u = (d.real / r[i], d.imag / r[i])
        for a in range(2):
            for c in range(2):
                H[a, c] += -Q[i] * r[i] ** -3 * ((1 if a == c else 0) - 3 * u[a] * u[c])

    def logQ(xv, yv, i):
        gg = C.geom(P, mp.mpc(xv, yv))
        return mp.log(abs(gg['lam'][i])) + 3 * mp.log(gg['r'][i])

    LD = mp.matrix(2, 2)
    for i in range(3):
        gr = (mp.diff(lambda s: logQ(s, x0.imag, i), x0.real), mp.diff(lambda s: logQ(x0.real, s, i), x0.imag))
        d = x0 - P[i]
        for a in range(2):
            for c in range(2):
                LD[a, c] += lam[i] * (d.real if a == 0 else d.imag) * gr[c]
    M = 2 * mp.mpf('0.5') * kappa * LD
    worst = max(worst, max(abs(H[i, j] - M[i, j]) for i in range(2) for j in range(2)) / max(abs(H[i, j]) for i in range(2) for j in range(2)))
    Lm = mp.matrix([[lam[i] * ((x0 - P[i]).real if a == 0 else (x0 - P[i]).imag) for i in range(3)] for a in range(2)])
    sv = mp.svd_r(Lm, compute_uv=False)
    worstL = max(worstL, 1 / min(sv[0], sv[1]))
print('(1) mixed signs, 120 critical points (60 in R+, 60 in R-): max rel |Hess V - 2 alpha kappa L DPsi| = %s; '
      'L has rank 2 (largest 1/sigma_min = %s)' % (mp.nstr(worst, 3), mp.nstr(worstL, 3)))


# ---------------------------------------------------------------- (2), (3) collinear
def V1(z, zs, q):
    return sum(qi / abs(z - zi) for zi, qi in zip(zs, q))


def dV1(z, zs, q):
    return sum(-qi * np.sign(z - zi) / (z - zi) ** 2 for zi, qi in zip(zs, q))


def crit3d(P3, q, nseed=400, spread=3.0):
    def G(x):
        r = np.linalg.norm(x[None] - P3, axis=1)
        w = q * r ** -3.0
        return (w[:, None] * (x[None] - P3)).sum(0) / (np.abs(w) * r).sum()
    out = []
    cen = P3.mean(0)
    for _ in range(nseed):
        x0 = cen + rng.normal(size=3) * spread * rng.choice([0.3, 1, 3])
        sol = root(G, x0, method='hybr', options={'xtol': 1e-15})
        x = sol.x
        if np.all(np.isfinite(x)) and np.linalg.norm(x - cen) < 1e3 and np.linalg.norm(G(x)) < 1e-12:
            if min(np.linalg.norm(x - p) for p in P3) > 1e-9 and all(np.linalg.norm(x - y) > 1e-7 for y in out):
                out.append(x)
    return out


hist_axis, hist_circ, max_nondeg, persist_fail, pert_counts = {}, {}, 0, 0, {}
for k in range(150):
    zs = np.sort(rng.uniform(-1, 1, 3))
    if np.min(np.diff(zs)) < 0.05:
        continue
    signs = [1, 1, -1]
    rng.shuffle(signs)
    q = np.array(signs) * np.exp(rng.normal(size=3) * 0.7)
    # on-axis zeros of V'
    axis = []
    grids = [np.linspace(zs[0] - 60, zs[0], 6001)[:-1], np.linspace(zs[0], zs[1], 6001)[1:-1],
             np.linspace(zs[1], zs[2], 6001)[1:-1], np.linspace(zs[2], zs[2] + 60, 6001)[1:]]
    for gz in grids:
        v = np.array([dV1(z, zs, q) for z in gz])
        for j in np.nonzero(np.sign(v[1:]) != np.sign(v[:-1]))[0]:
            axis.append(brentq(lambda z: dV1(z, zs, q), gz[j], gz[j + 1], xtol=1e-15))
    nondeg = 0
    for z in axis:
        f2 = sum(2 * qi / abs(z - zi) ** 3 for zi, qi in zip(zs, q))           # V''
        h = sum(-qi / abs(z - zi) ** 3 for zi, qi in zip(zs, q))               # transverse second derivative
        if abs(f2) > 1e-9 and abs(h) > 1e-9:
            nondeg += 1
    max_nondeg = max(max_nondeg, nondeg)
    hist_axis[len(axis)] = hist_axis.get(len(axis), 0) + 1
    # off-axis equilibria in the (z, rho) half-plane
    def Fz(v):
        z, rho = v[0], abs(v[1]) + 1e-300
        r = np.sqrt((z - zs) ** 2 + rho ** 2)
        w = q * r ** -3
        return np.array([np.sum(w * (z - zs)), np.sum(w) * rho]) / np.sum(np.abs(w) * r)
    circ = []
    for _ in range(150):
        sol = root(Fz, [rng.uniform(-3, 3), rng.uniform(0.05, 3)], method='hybr', options={'xtol': 1e-14})
        if np.all(np.isfinite(sol.x)) and abs(sol.x[1]) > 1e-6 and np.linalg.norm(Fz(sol.x)) < 1e-12 and abs(sol.x[0]) < 1e3:
            p = (sol.x[0], abs(sol.x[1]))
            if all(abs(p[0] - c_[0]) + abs(p[1] - c_[1]) > 1e-7 for c_ in circ):
                circ.append(p)
    hist_circ[len(circ)] = hist_circ.get(len(circ), 0) + 1
    # (3) perturb p3 off the axis by eps and count all critical points of the perturbed configuration
    if k % 3 == 0:
        P3 = np.c_[zs, np.zeros(3), np.zeros(3)]
        P3p = P3.copy()
        P3p[2, 1] = 1e-3
        pts = crit3d(P3p, q)
        pert_counts[len(pts)] = pert_counts.get(len(pts), 0) + 1
        for z in axis:
            f2 = sum(2 * qi / abs(z - zi) ** 3 for zi, qi in zip(zs, q))
            h = sum(-qi / abs(z - zi) ** 3 for zi, qi in zip(zs, q))
            if abs(f2) > 1e-9 and abs(h) > 1e-9:
                if not any(np.linalg.norm(x - np.array([z, 0, 0])) < 0.05 for x in pts):
                    persist_fail += 1
print('(2) collinear mixed charges (%d configurations): number of on-axis equilibria %s; max non-degenerate on-axis '
      'equilibria %d; number of off-axis equilibrium circles %s' % (sum(hist_axis.values()), dict(sorted(hist_axis.items())),
                                                                   max_nondeg, dict(sorted(hist_circ.items()))))
print('(3) after moving p3 off the axis by 1e-3: total critical points of the perturbed configurations %s; '
      'non-degenerate on-axis equilibria without a nearby perturbed critical point: %d'
      % (dict(sorted(pert_counts.items())), persist_fail))
