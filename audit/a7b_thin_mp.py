"""AUDIT item 7: is theta(M) convex for very thin isosceles triangles (A, A, 180 - 2A)?  50-digit arithmetic.

float64 (a67_thin_trend.py) gives a negative whitened curvature for A <= 0.1 deg; this checks whether that is real.
Strip coordinates: s = (theta1 - theta2)/sqrt2 (along the needle theta(M)), n = (theta1 + theta2)/sqrt2.  For each s
on a grid, n_+(s) = upper and n_-(s) = lower boundary of theta(M) on the vertical line (bisection, after checking that D
is unimodal in n on that line by sampling).  theta(M) is convex iff n_+ is concave and n_- is convex; tested by second
differences (two spacings).  Also reports the relative concavity margin min(-n_+'') * (s-extent)^2 / (n-extent).
usage: python a7b_thin_mp.py A_deg [alpha] [ns] [dps]"""
import sys
import time

import mpmath as mp

Ad = sys.argv[1] if len(sys.argv) > 1 else '0.1'
alpha = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf('0.5')
ns = int(sys.argv[3]) if len(sys.argv) > 3 else 60
mp.mp.dps = int(sys.argv[4]) if len(sys.argv) > 4 else 50
beta = 2 * alpha + 2
c = (beta - 1) / beta ** 2
A = [mp.mpf(Ad) * mp.pi / 180] * 2
A.append(mp.pi - 2 * A[0])
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


def golden_max(f, a, b, it=160):
    g = (mp.sqrt(5) - 1) / 2
    for _ in range(it):
        c1, c2 = b - g * (b - a), a + g * (b - a)
        if f(c1) > f(c2):
            b = c2
        else:
            a = c1
    return (a + b) / 2


t0 = time.time()
# admissible n-range on the axis s = 0: theta3 = 2pi - sqrt2 n in (A3, pi)  =>  n in (pi/sqrt2, (2pi - A3)/sqrt2)
nlo, nhi = mp.pi / r2, (2 * mp.pi - A[2]) / r2
nstar = golden_max(lambda n: D_sn(0, n), nlo, nhi)
print('A = %s deg, alpha = %s: D(theta*) = %s (should be 1/4), n* - pi/sqrt2 = %s, strip width %s'
      % (Ad, alpha, mp.nstr(D_sn(0, nstar), 15), mp.nstr(nstar - nlo, 6), mp.nstr(nhi - nlo, 6)))


def vertical_max(s):
    return golden_max(lambda n: D_sn(s, n), nlo, nhi, 120)


# s-extent of theta(M): largest s with max_n D(s, n) > c
a_, b_ = mp.mpf(0), mp.pi / r2
for _ in range(80):
    m = (a_ + b_) / 2
    if D_sn(m, vertical_max(m)) > c:
        a_ = m
    else:
        b_ = m
s_end = a_
print('s-extent of theta(M): %s (t = %.1fs)' % (mp.nstr(s_end, 10), time.time() - t0))


def boundary(s, nm):
    lo, hi = nm, nhi
    for _ in range(170):
        mid = (lo + hi) / 2
        if D_sn(s, mid) > c:
            lo = mid
        else:
            hi = mid
    up = lo
    lo, hi = nlo, nm
    for _ in range(170):
        mid = (lo + hi) / 2
        if D_sn(s, mid) > c:
            hi = mid
        else:
            lo = mid
    return up, hi


S = [s_end * mp.mpf(k) / ns * mp.mpf('0.995') for k in range(ns + 1)]
up, dn, unimodal_bad = [], [], 0
for s in S:
    nm = vertical_max(s)
    u, d = boundary(s, nm)
    # unimodality check along the vertical line (sampled): D increases up to nm and decreases after
    vals = [D_sn(s, d + (u - d) * mp.mpf(k) / 40) for k in range(-4, 45)]
    km = max(range(len(vals)), key=lambda k: vals[k])
    if any(vals[k + 1] < vals[k] for k in range(km)) or any(vals[k + 1] > vals[k] for k in range(km, len(vals) - 1)):
        unimodal_bad += 1
    up.append(u)
    dn.append(d)
print('boundary computed at %d values of s (t = %.1fs); vertical lines where D is not unimodal (sampled): %d'
      % (len(S), time.time() - t0, unimodal_bad))
hs = S[1] - S[0]
next_ = max(up) - min(dn)
for stride in (1, 2):
    d2u = [(up[k - stride] - 2 * up[k] + up[k + stride]) / (stride * hs) ** 2 for k in range(stride, len(S) - stride)]
    d2d = [(dn[k - stride] - 2 * dn[k] + dn[k + stride]) / (stride * hs) ** 2 for k in range(stride, len(S) - stride)]
    print('stride %d: n_+\'\' max %s (must be < 0), n_-\'\' min %s (must be > 0); #violations upper %d lower %d; '
          'relative margins min(-n_+\'\') s_end^2/n_extent = %s, min(n_-\'\') s_end^2/n_extent = %s'
          % (stride, mp.nstr(max(d2u), 5), mp.nstr(min(d2d), 5), sum(1 for v in d2u if v >= 0), sum(1 for v in d2d if v <= 0),
             mp.nstr(min(-v for v in d2u) * s_end ** 2 / next_, 5), mp.nstr(min(d2d) * s_end ** 2 / next_, 5)))
print('n_+(s) - n*, n_-(s) - n* at s/s_end = 0, .25, .5, .75, .995: %s'
      % [(mp.nstr(up[k] - nstar, 6), mp.nstr(dn[k] - nstar, 6)) for k in (0, ns // 4, ns // 2, 3 * ns // 4, ns)])
print('total time %.1fs' % (time.time() - t0))
