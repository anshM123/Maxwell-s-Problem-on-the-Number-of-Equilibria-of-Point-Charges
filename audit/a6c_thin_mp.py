"""AUDIT item 6, high-precision check of the fold-image self-intersections seen in float64 for thin isosceles shapes.

For an isosceles triangle (A1 = A2) the swap p1 <-> p2 reflects theta(M) in the line theta1 = theta2 and Gamma in the
line Psi1 = Psi2.  A real self-intersection of Gamma between mirror segments means that, at some boundary point of
theta(M) off the axis, Psi1 - Psi2 and theta1 - theta2 have opposite signs (for convex theta(M) this is impossible:
<Psi(a) - Psi(b), a - b> > 0 for a and its mirror image b).  Here the boundary points near the axis direction are
computed with 40 digits (bisection on D(x(theta)) = c) and the two signs compared.
usage: python a6c_thin_mp.py A_deg [alpha] [n]"""
import math
import sys

import mpmath as mp

mp.mp.dps = 40
Ad = mp.mpf(sys.argv[1]) if len(sys.argv) > 1 else mp.mpf('0.2')
alpha = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf('0.5')
n = int(sys.argv[3]) if len(sys.argv) > 3 else 120
beta = 2 * alpha + 2
c = (beta - 1) / beta ** 2
A = [Ad * mp.pi / 180, Ad * mp.pi / 180, mp.pi - 2 * Ad * mp.pi / 180]
L13 = mp.sin(A[1]) / mp.sin(A[2])
P = [(mp.mpf(0), mp.mpf(0)), (mp.mpf(1), mp.mpf(0)), (L13 * mp.cos(A[0]), L13 * mp.sin(A[0]))]


def x_of_theta(t1, t2):
    th = [t1, t2, 2 * mp.pi - t1 - t2]
    w = [mp.sin(A[i]) * mp.sin(th[i]) / mp.sin(th[i] - A[i]) for i in range(3)]
    s = sum(w)
    lam = [wi / s for wi in w]
    return (sum(lam[i] * P[i][0] for i in range(3)), sum(lam[i] * P[i][1] for i in range(3))), lam


def D_and_Psi(t1, t2):
    (x, y), lam = x_of_theta(t1, t2)
    Z = mp.mpc(0)
    L = []
    for i in range(3):
        dx, dy = x - P[i][0], y - P[i][1]
        r = mp.sqrt(dx * dx + dy * dy)
        u = mp.mpc(dx, dy) / r
        Z += lam[i] * u * u
        L.append(mp.log(lam[i]) + beta * mp.log(r))
    return (1 - abs(Z) ** 2) / 4, (L[0] - L[2], L[1] - L[2])


# theta*: by symmetry on the axis t1 = t2 = t; maximise D along the axis (the maximiser of D is the zero of Z)
f = lambda t: -D_and_Psi(t, t)[0]
lo, hi = (2 * mp.pi - A[2]) / 2 - mp.mpf('0.0'), mp.pi
# golden-section search on the admissible axis interval (pi - A3/2 ... ) : theta3 = 2pi - 2t in (A3, pi) => t in (pi/2, pi - A3/2)
a, b = mp.pi / 2, mp.pi - A[2] / 2
gr = (mp.sqrt(5) - 1) / 2
for _ in range(200):
    c1, c2 = b - gr * (b - a), a + gr * (b - a)
    if f(c1) < f(c2):
        b = c2
    else:
        a = c1
ts = (a + b) / 2
print('A = %s deg: theta* = (%s, %s), D(theta*) = %s' % (Ad, mp.nstr(ts, 12), mp.nstr(ts, 12), mp.nstr(-f(ts), 20)))
bad, rows = 0, []
for k in range(n):
    psi = mp.pi / 4 - mp.mpf('0.09') + mp.mpf('0.18') * k / (n - 1)
    if abs(psi - mp.pi / 4) < mp.mpf('1e-12'):
        continue
    e = (mp.cos(psi), mp.sin(psi))
    lo_, hi_ = mp.mpf(0), mp.mpf('1e-6')
    while D_and_Psi(ts + hi_ * e[0], ts + hi_ * e[1])[0] > c:
        lo_, hi_ = hi_, 2 * hi_
    for _ in range(120):
        mid = (lo_ + hi_) / 2
        if D_and_Psi(ts + mid * e[0], ts + mid * e[1])[0] > c:
            lo_ = mid
        else:
            hi_ = mid
    t1, t2 = ts + lo_ * e[0], ts + lo_ * e[1]
    _, (p1, p2) = D_and_Psi(t1, t2)
    s_theta, s_psi = mp.sign(t1 - t2), mp.sign(p1 - p2)
    if s_theta != s_psi:
        bad += 1
    rows.append((float(psi), float(t1 - t2), float(p1 - p2)))
print('boundary points with psi in pi/4 +- 0.09: %d; sign(Psi1 - Psi2) != sign(theta1 - theta2) at %d' % (len(rows), bad))
for r in rows[:: max(1, len(rows) // 12)]:
    print('  psi %.6f  theta1-theta2 %+.3e  Psi1-Psi2 %+.3e' % r)
