"""Mixed signs, Coulomb: the boundary-curvature inequality for M' in explicit form.

On the boundary cosh(Delta) = beta(s), with Delta = Lambda1(t1) - Lambda2(t2), s = t1 + t2,
    beta(s) = [(F-1) cos A3 - (F+1) cos(2s + A3)]/2   (F = K/2 - 1, Coulomb F = 5/4),
M' is convex at that point iff
    I := (c1 + c2)^2 (P S^2 + C beta'^2) + S [a1 (beta' c2 + S)^2 - a2 (S - beta' c1)^2]  >=  0,
with P = -beta'', C = cosh Delta, S = sinh Delta, c_i = 1/U_i, a_i = sinh Lambda_i/(cosh Lambda_i + cos A_i).
We trace the whole boundary for many shapes (rays from the excentre point), report min I / scale, and also the crude
bound I_crude = (c1+c2)^2 (P S^2 + C beta'^2) - |S| (|a1| (beta' c2 + S)^2 + |a2| (S - beta' c1)^2).
Cross-check: the sign of I must agree with the bordered Hessian of E (code/mixed_boundary.py).
usage: python mixed_boundary_ineq.py [n_shapes] [min_angle_deg] [K]"""
import math
import sys

import numpy as np

from mixed_boundary import E_derivs

nsh = int(sys.argv[1]) if len(sys.argv) > 1 else 300
amin = math.radians(float(sys.argv[2])) if len(sys.argv) > 2 else math.radians(0.5)
K = float(sys.argv[3]) if len(sys.argv) > 3 else 4.5
F = K / 2 - 1
rng = np.random.default_rng(7)
worst, worst_crude, disagree, npts = np.inf, np.inf, 0, 0
wshape = wcshape = None
for it in range(nsh):
    while True:
        A = rng.dirichlet(np.ones(3) * 0.6) * math.pi
        if A.min() >= amin:
            break
    A1, A2, A3 = A
    c0 = np.array([A1 / 2, A2 / 2])
    for phi in np.linspace(0, 2 * math.pi, 240, endpoint=False):
        v = np.array([math.cos(phi), math.sin(phi)])
        smax = min((A1 - c0[0]) / v[0] if v[0] > 0 else (-c0[0] / v[0] if v[0] < 0 else np.inf),
                   (A2 - c0[1]) / v[1] if v[1] > 0 else (-c0[1] / v[1] if v[1] < 0 else np.inf))
        lo, hi = 0.0, smax * (1 - 1e-9)
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if E_derivs(c0[0] + mid * v[0], c0[1] + mid * v[1], A1, A2)[0] < K:
                lo = mid
            else:
                hi = mid
        t1, t2 = c0 + lo * v
        s = t1 + t2
        L1 = math.log(math.sin(t1) / math.sin(A1 - t1))
        L2 = math.log(math.sin(t2) / math.sin(A2 - t2))
        D = L1 - L2
        C, S = math.cosh(D), math.sinh(D)
        beta = ((F - 1) * math.cos(A3) - (F + 1) * math.cos(2 * s + A3)) / 2
        if abs(beta / C - 1) > 1e-6:
            continue                                   # not on the level set (numerical), skip
        bp = (F + 1) * math.sin(2 * s + A3)
        P = -2 * (F + 1) * math.cos(2 * s + A3)
        c1 = math.sin(A1) / (2 * (math.cosh(L1) + math.cos(A1)))
        c2 = math.sin(A2) / (2 * (math.cosh(L2) + math.cos(A2)))
        a1 = math.sinh(L1) / (math.cosh(L1) + math.cos(A1))
        a2 = math.sinh(L2) / (math.cosh(L2) + math.cos(A2))
        main = (c1 + c2) ** 2 * (P * S * S + C * bp * bp)
        I = main + S * (a1 * (bp * c2 + S) ** 2 - a2 * (S - bp * c1) ** 2)
        Icr = main - abs(S) * (abs(a1) * (bp * c2 + S) ** 2 + abs(a2) * (S - bp * c1) ** 2)
        scale = main + abs(S) * (abs(a1) * (bp * c2 + S) ** 2 + abs(a2) * (S - bp * c1) ** 2)
        E, E1, E2, E11, E12, E22 = E_derivs(t1, t2, A1, A2)
        Q = E11 * E2 ** 2 - 2 * E12 * E1 * E2 + E22 * E1 ** 2
        npts += 1
        if (Q > 0) != (I > 0) and abs(I) / scale > 1e-6:
            disagree += 1
        if I / scale < worst:
            worst, wshape = I / scale, np.degrees(A).round(2).tolist()
        if Icr / scale < worst_crude:
            worst_crude, wcshape = Icr / scale, np.degrees(A).round(2).tolist()
print('K = %.2f, %d boundary points: min I/scale = %.4f at %s;  min I_crude/scale = %.4f at %s;  sign disagreements '
      'with the bordered Hessian of E: %d' % (K, npts, worst, wshape, worst_crude, wcshape, disagree))
