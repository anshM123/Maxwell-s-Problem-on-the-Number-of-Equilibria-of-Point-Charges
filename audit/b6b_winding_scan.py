"""AUDIT_2 item 6 (all charges at once, per shape): number of minima = winding number of the fold image grad H(dM')
around psi_q (grad H is orientation preserving on M' since Hess H < 0 there).  Float64, auditor's own code.

For each shape (A1 <= A2 by the 1<->2 mirror; uniform grid plus near-degenerate extras):
  - the boundary of M' in (s, Delta) coordinates: s = (A12 + rho psi_max)/2, Delta = +-arccosh beta(s), rho in [-1, 1];
    t1 from Delta by vectorised bisection on the line t1 + t2 = s (no quadratic formula);
  - Gamma = grad H(boundary), H' = (phi1(t1) - phi3(s), phi2(t2) - phi3(s));
  - interior test points (s, rho' B(s)), |rho'| <= 0.97, mapped by grad H; winding numbers of Gamma around them;
  - also: M' sampled on a fine (s, Delta) grid is checked against E < K directly (consistency of the description).
Reports the range of winding numbers (|w| >= 2 anywhere would mean two minima for some charges).
usage: python b6b_winding_scan.py step_deg"""
import math
import sys
import time

import numpy as np

step = float(sys.argv[1]) if len(sys.argv) > 1 else 3.0
K, b = 4.5, 2.0
F = K / 2 - 1


def Lam(t, a):
    return np.log(np.sin(t) / np.sin(a - t))


def U(t, a):
    return 1 / np.tan(t) + 1 / np.tan(a - t)


def t_of_sD(A1, A2, s, D):
    lo = np.maximum(0.0, s - A2)
    hi = np.minimum(A1, s)
    for _ in range(80):
        m = 0.5 * (lo + hi)
        v = Lam(m, A1) - Lam(s - m, A2)
        lo, hi = np.where(v < D, m, lo), np.where(v < D, hi, m)
    t1 = 0.5 * (lo + hi)
    return t1, s - t1


def gradH(A1, A2, t1, t2):
    A3 = math.pi - A1 - A2
    s = t1 + t2
    p1 = np.log(np.sin(t1)) + b * np.log(np.sin(A1 - t1))
    p2 = np.log(np.sin(t2)) + b * np.log(np.sin(A2 - t2))
    p3 = np.log(np.sin(s)) + b * np.log(np.sin(s + A3))
    return np.stack([p1 - p3, p2 - p3], -1)


def analyse(A1, A2, nb=1500):
    A3 = math.pi - A1 - A2
    A12 = A1 + A2
    pm = math.acos((2 - (F - 1) * math.cos(A3)) / (F + 1))
    rho = np.linspace(-1, 1, nb)
    s = (A12 + rho * pm) / 2
    beta = ((F - 1) * math.cos(A3) + (F + 1) * np.cos(rho * pm)) / 2
    B = np.arccosh(np.maximum(beta, 1.0))
    tu = t_of_sD(A1, A2, s, B)
    tl = t_of_sD(A1, A2, s[::-1], -B[::-1])
    Gam = np.vstack([gradH(A1, A2, *tu), gradH(A1, A2, *tl)])
    # interior points
    ri = np.linspace(-0.97, 0.97, 40)
    rr = np.linspace(-0.97, 0.97, 41)
    S_, R_ = np.meshgrid((A12 + ri * pm) / 2, rr, indexing='ij')
    bet_i = ((F - 1) * math.cos(A3) + (F + 1) * np.cos(2 * S_ - A12)) / 2
    Bi = np.arccosh(np.maximum(bet_i, 1.0))
    ti = t_of_sD(A1, A2, S_.ravel(), (R_ * Bi).ravel())
    W = gradH(A1, A2, *ti)
    Ei = np.array([U(a, A1) * U(c, A2) - U(a + c, A12) * (U(a, A1) + U(c, A2)) for a, c in zip(*ti)])
    inc = int(np.sum(Ei >= K))
    dv = Gam[None, :, :] - W[:, None, :]
    ang = np.arctan2(dv[..., 1], dv[..., 0])
    da = np.diff(np.concatenate([ang, ang[:, :1]], 1), axis=1)
    da = (da + math.pi) % (2 * math.pi) - math.pi
    wind = da.sum(1) / (2 * math.pi)
    wr = np.round(wind)
    return int(np.abs(wr).min()), int(np.abs(wr).max()), float(np.max(np.abs(wind - wr))), inc


shapes = []
g = np.arange(step / 2, 180, step)
for a1 in g:
    for a2 in g:
        if a1 <= a2 and 180 - a1 - a2 > 0.05:
            shapes.append((math.radians(a1), math.radians(a2), 'grid'))
for e in (1e-3, 1e-2, 0.05, 0.1, 0.2):
    for lam in (0.01, 0.1, 0.3, 0.5):
        shapes.append((lam * e, (1 - lam) * e, 'A12=%g' % e))
for r in (1e-3, 1e-2, 0.05):
    for th in (0.01, 0.3, 0.7, 0.99):
        shapes.append((th * r, math.pi - r, 'A2 near pi, r=%g' % r))      # mirror of A1 near pi
for r in (1e-3, 1e-2, 0.05):
    for lam in (0.05, 0.3, 0.5):
        shapes.append((lam * (math.pi - r), (1 - lam) * (math.pi - r), 'A3=%g' % r))
t0 = time.time()
wmin, wmax, devmax, incons, flagged = 99, 0, 0.0, 0, []
for A1, A2, tag in shapes:
    a, bmax, dev, inc = analyse(A1, A2)
    wmin, wmax, devmax, incons = min(wmin, a), max(wmax, bmax), max(devmax, dev), incons + inc
    if bmax >= 2 or a == 0 or dev > 0.1 or inc:
        flagged.append((tag, math.degrees(A1), math.degrees(A2), a, bmax, dev, inc))
print('%d shapes (grid step %.1f deg with A1 <= A2, plus near-degenerate extras) in %.0fs' % (len(shapes), step, time.time() - t0))
print('|winding numbers| at interior images: min %d, max %d; max distance from an integer %.2e; interior grid points '
      'with E >= K (description inconsistency): %d' % (wmin, wmax, devmax, incons))
print('flagged shapes: %d' % len(flagged))
for f in flagged[:20]:
    print('   %s angles (%.4f, %.4f) deg: |w| in [%d, %d], dev %.2e, inconsistent %d' % f)
