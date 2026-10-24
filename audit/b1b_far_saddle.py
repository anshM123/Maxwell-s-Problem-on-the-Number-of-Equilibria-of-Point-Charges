"""AUDIT_2 item 1 (index count, small total charge): for Q_tot = eps != 0 the large-circle winding is 2 at radii below
~|m|/|Q_tot| and 1 beyond, with a saddle near the transition radius.  Confirms the one anomalous tally entry of
b1_sections1to5.py is a finite-radius effect, consistent with s0 = m0 + 2 for Q_tot != 0.
usage: python b1b_far_saddle.py"""
import math

import numpy as np
from scipy.optimize import root

rng = np.random.default_rng(5)
P = np.array([[0.0, 0.0], [1.0, 0.0], [0.3, 0.8]])


def winding(fun, curve):
    v = np.array([fun(c_) for c_ in curve])
    a = np.arctan2(v[:, 1], v[:, 0])
    da = np.diff(np.concatenate([a, a[:1]]))
    da = (da + math.pi) % (2 * math.pi) - math.pi
    return da.sum() / (2 * math.pi)


for eps in (1e-1, 1e-2, 1e-3):
    q = np.array([1.0, 0.7, -1.7 + eps])
    Qt = q.sum()
    m = (q[:, None] * P).sum(0)
    gradV = lambda y: sum(-qi * np.linalg.norm(y - p) ** -3 * (y - p) for p, qi in zip(P, q))
    G = lambda y: gradV(y) / np.linalg.norm(gradV(y)) * min(1.0, 1.0)
    tt = np.linspace(0, 2 * math.pi, 40000, endpoint=False)
    circ = np.c_[np.cos(tt), np.sin(tt)]
    ws = [(R, round(winding(gradV, P.mean(0) + R * circ), 3)) for R in (10, 0.1 * np.linalg.norm(m) / abs(Qt),
                                                                         10 * np.linalg.norm(m) / abs(Qt))]
    # far critical point: search along the line through the dipole axis at large radii
    def Fs(y):
        r = np.linalg.norm(y[None] - P, axis=1)
        w = q * r ** -3
        return (w[:, None] * (y[None] - P)).sum(0) / (np.abs(w) * r).sum()
    far = []
    R0 = np.linalg.norm(m) / abs(Qt)
    for th in np.linspace(0, 2 * math.pi, 72, endpoint=False):
        for rr in (0.3 * R0, R0, 3 * R0):
            sol = root(Fs, P.mean(0) + rr * np.array([math.cos(th), math.sin(th)]), method='hybr', options={'xtol': 1e-15})
            if np.all(np.isfinite(sol.x)) and np.linalg.norm(Fs(sol.x)) < 1e-12 and np.linalg.norm(sol.x) > 5:
                if all(np.linalg.norm(sol.x - f) > 1e-6 * np.linalg.norm(sol.x) for f in far):
                    far.append(sol.x)
    print('Q_tot = %.0e, |m|/|Q_tot| = %.3g: windings (radius, w) %s; far critical points at radii %s'
          % (Qt, R0, [(round(float(R), 3), w) for R, w in ws], [round(float(np.linalg.norm(f)), 3) for f in far]))
