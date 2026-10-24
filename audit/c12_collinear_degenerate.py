"""AUDIT_3, item 8 (NOTE N9): an isolated equilibrium with zero Hessian for collinear charges of mixed signs (Coulomb).

Charges at t = (0, 1, 3) on the x1-axis of R^3; the point x* = (-1, 0, 0). With w = (t2 - t3, t3 - t1, t1 - t2) = (-2, 3, -1)
(so sum w = 0 and sum w t = 0) and q_i = w_i |x1* - t_i|^beta (beta = 3), x* is critical and Hess V(x*) = 0.
Off-axis critical points (circles) are the solutions with u = rho^2 > 0 of F1 = sum q_i r_i^-3 = 0,
F2 = sum q_i t_i r_i^-3 = 0, r_i^2 = (x1 - t_i)^2 + u. If the Jacobian of (F1, F2) in (x1, u) is non-singular at (x1*, 0),
no circle of critical points comes near x*; on the axis, x* is an isolated zero of the analytic f'(x1).
So x* is an isolated critical point that Theorem 1.1(b) (non-degenerate points only) does not count.
Output: audit/c12_collinear_degenerate.out"""
import mpmath as mp

mp.mp.dps = 40
t = [mp.mpf(0), mp.mpf(1), mp.mpf(3)]
x1s = mp.mpf(-1)
w = [t[1] - t[2], t[2] - t[0], t[0] - t[1]]
q = [w[i] * abs(x1s - t[i]) ** 3 for i in range(3)]


def V(x, y, z):
    return sum(q[i] / mp.sqrt((x - t[i]) ** 2 + y ** 2 + z ** 2) for i in range(3))


g = [mp.diff(V, (x1s, 0, 0), o) for o in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
H = [[mp.diff(V, (x1s, 0, 0), tuple(int(k == i) + int(k == j) for k in range(3))) for j in range(3)] for i in range(3)]


def F(x1, u):
    r3 = [((x1 - t[i]) ** 2 + u) ** mp.mpf(-1.5) for i in range(3)]
    return sum(q[i] * r3[i] for i in range(3)), sum(q[i] * t[i] * r3[i] for i in range(3))


J = mp.matrix([[mp.diff(lambda a, b: F(a, b)[k], (x1s, 0), (1, 0)), mp.diff(lambda a, b: F(a, b)[k], (x1s, 0), (0, 1))]
               for k in range(2)])
f = lambda s: V(s, 0, 0)
derivs = [mp.diff(f, x1s, n) for n in (1, 2, 3)]
print('charges q = %s (mixed signs), positions t = (0, 1, 3), x* = (-1, 0, 0)' % [mp.nstr(v, 6) for v in q])
print('grad V(x*) = %s' % [mp.nstr(v, 3) for v in g])
print('max |Hess V(x*)| entry = %s' % mp.nstr(max(abs(v) for row in H for v in row), 3))
print('Jacobian of (F1, F2) in (x1, u) at (x1*, 0): det = %s (non-zero: no circle of critical points near x*)'
      % mp.nstr(mp.det(J), 8))
print("on the axis: f'(x*), f''(x*), f'''(x*) = %s (f' vanishes to order exactly 2 at x*, so x* is an isolated zero of f')"
      % [mp.nstr(v, 6) for v in derivs])
