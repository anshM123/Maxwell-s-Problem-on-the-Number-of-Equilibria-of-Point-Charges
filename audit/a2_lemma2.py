"""AUDIT item 2: Lemma 2, det DPsi = beta^2 (D - c) / (2 Delta lambda1 lambda2 lambda3), c = (beta-1)/beta^2.

Independent of code/symbolic_I2.py:
 (2a) symbolic, in a different parametrization: arbitrary vertices p_i = (a_i, b_i) (6 free symbols), barycentric
      coordinates by Cramer's rule, D computed as det S (S = sum lambda_i d_i d_i^T / r_i^2, not by Cauchy-Binet), the
      Jacobian of Psi = (L1 - L3, L2 - L3), L_i = log lambda_i + (beta/2) log r_i^2, taken by sympy and then
      evaluated at x = 0 (both sides are translation invariant, so x = 0 with arbitrary vertices is the general case).
      Delta here is the SIGNED area (half the Cramer determinant).  The numerator of lhs - rhs is expanded exactly.
 (2b) 50-digit numerical check (mpmath, numerical differentiation of Psi), random triangles with BOTH orientations,
      random points and random beta > 2: ratio det DPsi / [beta^2 (D - c)/(2 Delta_signed l1 l2 l3)] and the same
      ratio with Delta = |area|.
usage: python a2_lemma2.py"""
import time

import mpmath as mp
import numpy as np
import sympy as sp

# ------------------------------------------------------------------ (2a) symbolic
t0 = time.time()
a1, b1, a2, b2, a3, b3, X, Y = sp.symbols('a1 b1 a2 b2 a3 b3 X Y', real=True)
beta = sp.symbols('beta', positive=True)
P = [(a1, b1), (a2, b2), (a3, b3)]
Mx = sp.Matrix([[a1, a2, a3], [b1, b2, b3], [1, 1, 1]])
det2 = Mx.det()                       # = 2 * signed area
lam = []
for i in range(3):
    Mi = Mx.copy()
    Mi[:, i] = sp.Matrix([X, Y, 1])
    lam.append(Mi.det() / det2)      # Cramer's rule
d = [sp.Matrix([X - p[0], Y - p[1]]) for p in P]
r2 = [di.dot(di) for di in d]
L = [sp.log(lam[i]) + beta / 2 * sp.log(r2[i]) for i in range(3)]
J = sp.Matrix([[sp.diff(L[0] - L[2], X), sp.diff(L[0] - L[2], Y)],
               [sp.diff(L[1] - L[2], X), sp.diff(L[1] - L[2], Y)]])
at0 = {X: 0, Y: 0}
lhs = sp.together(J.subs(at0).det())
S = sp.zeros(2, 2)
for i in range(3):
    S += lam[i] * d[i] * d[i].T / r2[i]
Dexpr = S.subs(at0).det()
c = (beta - 1) / beta ** 2
lam0 = [li.subs(at0) for li in lam]
rhs = beta ** 2 * (Dexpr - c) / (det2 * lam0[0] * lam0[1] * lam0[2])     # 2 Delta_signed = det2
num = sp.numer(sp.together(lhs - rhs))
num_exp = sp.expand(num)
print('(2a) general vertices, x = 0, signed Delta: expanded numerator of lhs - rhs = %s   (%.1fs)' % (num_exp, time.time() - t0))
print('     sanity: number of monomials in the expanded numerator of lhs alone: %d'
      % len(sp.Add.make_args(sp.expand(sp.numer(lhs)))))

# ------------------------------------------------------------------ (2b) 50-digit numerics, both orientations
mp.mp.dps = 50
rng = np.random.default_rng(7)


def Psi(xv, yv, Pm, bet):
    det = (Pm[1][0] - Pm[0][0]) * (Pm[2][1] - Pm[0][1]) - (Pm[2][0] - Pm[0][0]) * (Pm[1][1] - Pm[0][1])
    lamv = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        lamv.append(((Pm[j][0] - xv) * (Pm[k][1] - yv) - (Pm[j][1] - yv) * (Pm[k][0] - xv)) / det)
    Lv = [mp.log(lamv[i]) + bet / 2 * mp.log((xv - Pm[i][0]) ** 2 + (yv - Pm[i][1]) ** 2) for i in range(3)]
    return Lv[0] - Lv[2], Lv[1] - Lv[2], lamv, det


worst_signed, worst_abs_ccw, worst_abs_cw = mp.mpf(0), mp.mpf(0), mp.mpf(0)
signs_cw = set()
for k in range(300):
    Pm = [[mp.mpf(float(v)) for v in rng.normal(size=2)] for _ in range(3)]
    det = (Pm[1][0] - Pm[0][0]) * (Pm[2][1] - Pm[0][1]) - (Pm[2][0] - Pm[0][0]) * (Pm[1][1] - Pm[0][1])
    if abs(det) < 0.05:
        continue
    w = rng.dirichlet(np.ones(3))
    xv = sum(mp.mpf(float(w[i])) * Pm[i][0] for i in range(3))
    yv = sum(mp.mpf(float(w[i])) * Pm[i][1] for i in range(3))
    bet = mp.mpf(2) + mp.mpf(float(rng.uniform(0.05, 6)))
    f1 = lambda u, v: Psi(u, v, Pm, bet)[0]
    f2 = lambda u, v: Psi(u, v, Pm, bet)[1]
    Jm = mp.matrix([[mp.diff(f1, (xv, yv), (1, 0)), mp.diff(f1, (xv, yv), (0, 1))],
                    [mp.diff(f2, (xv, yv), (1, 0)), mp.diff(f2, (xv, yv), (0, 1))]])
    dJ = mp.det(Jm)
    _, _, lamv, _ = Psi(xv, yv, Pm, bet)
    U = []
    for i in range(3):
        dx, dy = xv - Pm[i][0], yv - Pm[i][1]
        rr = mp.sqrt(dx * dx + dy * dy)
        U.append((dx / rr, dy / rr))
    Sm = mp.matrix(2, 2)
    for i in range(3):
        Sm[0, 0] += lamv[i] * U[i][0] ** 2
        Sm[0, 1] += lamv[i] * U[i][0] * U[i][1]
        Sm[1, 1] += lamv[i] * U[i][1] ** 2
    Sm[1, 0] = Sm[0, 1]
    Dv = mp.det(Sm)
    cv = (bet - 1) / bet ** 2
    if abs(Dv - cv) < mp.mpf('1e-6'):
        continue
    pred_signed = bet ** 2 * (Dv - cv) / (det * lamv[0] * lamv[1] * lamv[2])
    pred_abs = bet ** 2 * (Dv - cv) / (abs(det) * lamv[0] * lamv[1] * lamv[2])
    worst_signed = max(worst_signed, abs(dJ / pred_signed - 1))
    if det > 0:
        worst_abs_ccw = max(worst_abs_ccw, abs(dJ / pred_abs - 1))
    else:
        signs_cw.add(int(mp.sign(dJ / pred_abs)))
print('(2b) 50-digit numerics: max |ratio - 1| with signed Delta: %s;  with Delta = |area|: ccw-labelled %s, '
      'cw-labelled ratio signs %s' % (mp.nstr(worst_signed, 3), mp.nstr(worst_abs_ccw, 3), sorted(signs_cw)))
