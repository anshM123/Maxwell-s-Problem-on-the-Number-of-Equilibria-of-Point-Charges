"""AUDIT item 7: 50-digit convexity test of theta(M) for thin triangles with two small, possibly unequal angles
(A1, A2 small, A3 = 180 - A1 - A2), generalising a7b_thin_mp.py (no symmetry assumed).

Strip coordinates s = (theta1 - theta2)/sqrt2, n = (theta1 + theta2)/sqrt2.  theta* = argmax D (nested golden search);
the s-range of theta(M) is found on both sides of s*; for each s on a grid the upper/lower boundary n_+(s), n_-(s) are
found by bisection (after a sampled unimodality check in n); theta(M) is convex iff n_+ is concave and n_- convex.
usage: python a7c_thin_general_mp.py A1_deg A2_deg [alpha] [ns] [dps]"""
import sys
import time

import mpmath as mp

A1d, A2d = sys.argv[1], sys.argv[2]
alpha = mp.mpf(sys.argv[3]) if len(sys.argv) > 3 else mp.mpf('0.5')
ns = int(sys.argv[4]) if len(sys.argv) > 4 else 50
mp.mp.dps = int(sys.argv[5]) if len(sys.argv) > 5 else 50
beta = 2 * alpha + 2
c = (beta - 1) / beta ** 2
A = [mp.mpf(A1d) * mp.pi / 180, mp.mpf(A2d) * mp.pi / 180]
A.append(mp.pi - A[0] - A[1])
L13 = mp.sin(A[1]) / mp.sin(A[2])
P = [(mp.mpf(0), mp.mpf(0)), (mp.mpf(1), mp.mpf(0)), (L13 * mp.cos(A[0]), L13 * mp.sin(A[0]))]
r2 = mp.sqrt(2)


def D_sn(s, n):
    t1, t2 = (n + s) / r2, (n - s) / r2
    th = [t1, t2, 2 * mp.pi - t1 - t2]
    if not all(A[i] < th[i] < mp.pi for i in range(3)):
        return mp.mpf(-1)
    w = [mp.sin(A[i]) * mp.sin(th[i]) / mp.sin(th[i] - A[i]) for i in range(3)]
    sw = sum(w)
    lam = [wi / sw for wi in w]
    x = sum(lam[i] * P[i][0] for i in range(3))
    y = sum(lam[i] * P[i][1] for i in range(3))
    Z = mp.mpc(0)
    for i in range(3):
        u = mp.mpc(x - P[i][0], y - P[i][1])
        Z += lam[i] * u * u / (u.real ** 2 + u.imag ** 2)
    return (1 - abs(Z) ** 2) / 4


def golden_max(f, a, b, it=150):
    g = (mp.sqrt(5) - 1) / 2
    for _ in range(it):
        c1, c2 = b - g * (b - a), a + g * (b - a)
        if f(c1) > f(c2):
            b = c2
        else:
            a = c1
    return (a + b) / 2


t0 = time.time()
# theta3 = 2pi - sqrt2 n in (A3, pi)  =>  n in (pi/sqrt2, (2pi - A3)/sqrt2); s limited by theta1, theta2 in (A_i, pi)
nlo, nhi = mp.pi / r2, (2 * mp.pi - A[2]) / r2
slo, shi = -mp.pi / r2, mp.pi / r2


def vmax(s):
    return golden_max(lambda n: D_sn(s, n), nlo, nhi, 110)


sstar = golden_max(lambda s: D_sn(s, vmax(s)), slo, shi, 110)
nstar = vmax(sstar)
print('A = (%s, %s, %s) deg, alpha %s: D(theta*) = %s, s* = %s, n* - nlo = %s, strip width %s'
      % (A1d, A2d, mp.nstr(A[2] * 180 / mp.pi, 8), alpha, mp.nstr(D_sn(sstar, nstar), 15), mp.nstr(sstar, 10),
         mp.nstr(nstar - nlo, 6), mp.nstr(nhi - nlo, 6)))


def s_edge(direction):
    a_, b_ = sstar, (shi if direction > 0 else slo)
    for _ in range(70):
        m = (a_ + b_) / 2
        if D_sn(m, vmax(m)) > c:
            a_ = m
        else:
            b_ = m
    return a_


sa, sb = s_edge(-1), s_edge(+1)
print('s-range of theta(M): [%s, %s] (t = %.1fs)' % (mp.nstr(sa, 10), mp.nstr(sb, 10), time.time() - t0))
S = [sa + (sb - sa) * (mp.mpf('0.0025') + mp.mpf('0.995') * mp.mpf(k) / ns) for k in range(ns + 1)]
up, dn, bad_uni = [], [], 0
for s in S:
    nm = vmax(s)
    lo, hi = nm, nhi
    for _ in range(170):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if D_sn(s, mid) > c else (lo, mid)
    u = lo
    lo, hi = nlo, nm
    for _ in range(170):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if D_sn(s, mid) > c else (mid, hi)
    d = hi
    vals = [D_sn(s, d + (u - d) * mp.mpf(k) / 40) for k in range(-4, 45)]
    km = max(range(len(vals)), key=lambda k: vals[k])
    if any(vals[k + 1] < vals[k] for k in range(km)) or any(vals[k + 1] > vals[k] for k in range(km, len(vals) - 1)):
        bad_uni += 1
    up.append(u)
    dn.append(d)
hs = S[1] - S[0]
next_ = max(up) - min(dn)
L = sb - sa
print('boundary at %d values of s (t = %.1fs); vertical lines with D not unimodal (sampled): %d'
      % (len(S), time.time() - t0, bad_uni))
for stride in (1, 2):
    d2u = [(up[k - stride] - 2 * up[k] + up[k + stride]) / (stride * hs) ** 2 for k in range(stride, len(S) - stride)]
    d2d = [(dn[k - stride] - 2 * dn[k] + dn[k + stride]) / (stride * hs) ** 2 for k in range(stride, len(S) - stride)]
    print('stride %d: #violations upper %d lower %d; relative margins min(-n_+\'\') (L/2)^2/n_extent = %s, '
          'min(n_-\'\') (L/2)^2/n_extent = %s'
          % (stride, sum(1 for v in d2u if v >= 0), sum(1 for v in d2d if v <= 0),
             mp.nstr(min(-v for v in d2u) * (L / 2) ** 2 / next_, 5), mp.nstr(min(d2d) * (L / 2) ** 2 / next_, 5)))
print('total time %.1fs' % (time.time() - t0))
