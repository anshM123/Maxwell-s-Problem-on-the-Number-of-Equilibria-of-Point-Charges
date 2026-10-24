"""Mixed signs: up to which level is E convex?  For each shape, k_conv(A) = inf{E(t) : Hess E(t) not positive
semidefinite} over the rectangle (exact second derivatives), and k_qc(A) = the lowest level whose sublevel-set boundary
has negative curvature (quasi-convexity failure).  Coulomb needs E convex (or quasi-convex) on {E <= 4.5}.
Normalised coordinates u_i = t_i / A_i are used for the eigenvalue test (convexity is affine invariant).
usage: python mixed_convexity_level.py [n_shapes] [min_angle_deg] [seed]"""
import math
import sys

import numpy as np

from mixed_boundary import E_derivs

nsh = int(sys.argv[1]) if len(sys.argv) > 1 else 300
amin = math.radians(float(sys.argv[2])) if len(sys.argv) > 2 else math.radians(0.5)
rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 0)
worst_conv, worst_qc = (np.inf, None), (np.inf, None)
for it in range(nsh):
    while True:
        A = rng.dirichlet(np.ones(3) * 0.6) * math.pi
        if A.min() >= amin:
            break
    A1, A2 = A[0], A[1]
    g = np.linspace(0, 1, 402)[1:-1]
    U1, U2 = np.meshgrid(g, g)
    t1, t2 = U1.ravel() * A1, U2.ravel() * A2
    E, E1, E2, E11, E12, E22 = E_derivs(t1, t2, A1, A2)
    # Hessian in normalised coordinates
    h11, h12, h22 = E11 * A1 * A1, E12 * A1 * A2, E22 * A2 * A2
    tr, det = h11 + h22, h11 * h22 - h12 * h12
    lmin = tr / 2 - np.sqrt(np.maximum(tr * tr / 4 - det, 0))
    bad = lmin < -1e-9 * np.abs(tr)
    kc = E[bad].min() if bad.any() else np.inf
    # bordered Hessian (quasi-convexity) in normalised coordinates
    e1, e2 = E1 * A1, E2 * A2
    Q = h11 * e2 ** 2 - 2 * h12 * e1 * e2 + h22 * e1 ** 2
    badq = Q < -1e-9 * (np.abs(h11) + np.abs(h22) + 2 * np.abs(h12)) * (e1 ** 2 + e2 ** 2)
    kq = E[badq].min() if badq.any() else np.inf
    if kc < worst_conv[0]:
        worst_conv = (kc, np.degrees(A).round(3).tolist())
    if kq < worst_qc[0]:
        worst_qc = (kq, np.degrees(A).round(3).tolist())
print('over %d shapes (min angle %.2f deg): inf E-level where E stops being convex = %.4f at %s'
      % (nsh, math.degrees(amin), *worst_conv))
print('                                   inf E-level where sublevel sets stop being convex = %.4f at %s' % worst_qc)
print('(Coulomb needs 4.5; alpha needs (2 alpha + 2)^2 / (2 alpha + 1))')
