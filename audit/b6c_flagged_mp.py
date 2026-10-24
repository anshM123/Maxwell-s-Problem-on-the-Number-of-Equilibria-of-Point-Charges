"""AUDIT_2 item 6: high-precision re-check of shapes flagged by the float64 winding scan (winding 0 at an interior image,
impossible in exact arithmetic).  Same construction as b6b_winding_scan.py in 40-digit mpmath.
usage: python b6c_flagged_mp.py A1_deg A2_deg [nb] [ni]"""
import sys

import mpmath as mp

mp.mp.dps = 40
A1 = mp.mpf(sys.argv[1]) * mp.pi / 180
A2 = mp.mpf(sys.argv[2]) * mp.pi / 180
nb = int(sys.argv[3]) if len(sys.argv) > 3 else 400
ni = int(sys.argv[4]) if len(sys.argv) > 4 else 15
K, b = mp.mpf(9) / 2, mp.mpf(2)
F = K / 2 - 1
A3 = mp.pi - A1 - A2
A12 = A1 + A2
Lam = lambda t, a: mp.log(mp.sin(t) / mp.sin(a - t))


def t_of_sD(s, D):
    lo, hi = max(mp.mpf(0), s - A2), min(A1, s)
    for _ in range(170):
        m = (lo + hi) / 2
        if Lam(m, A1) - Lam(s - m, A2) < D:
            lo = m
        else:
            hi = m
    t1 = (lo + hi) / 2
    return t1, s - t1


def gradH(t1, t2):
    s = t1 + t2
    p3 = mp.log(mp.sin(s)) + b * mp.log(mp.sin(s + A3))
    return (mp.log(mp.sin(t1)) + b * mp.log(mp.sin(A1 - t1)) - p3, mp.log(mp.sin(t2)) + b * mp.log(mp.sin(A2 - t2)) - p3)


pm = mp.acos((2 - (F - 1) * mp.cos(A3)) / (F + 1))
beta = lambda s: ((F - 1) * mp.cos(A3) + (F + 1) * mp.cos(2 * s - A12)) / 2
Gam = []
for j in range(nb):
    rho = -1 + 2 * mp.mpf(j) / (nb - 1)
    s = (A12 + rho * pm) / 2
    Gam.append(gradH(*t_of_sD(s, mp.acosh(max(beta(s), mp.mpf(1))))))
for j in range(nb):
    rho = 1 - 2 * mp.mpf(j) / (nb - 1)
    s = (A12 + rho * pm) / 2
    Gam.append(gradH(*t_of_sD(s, -mp.acosh(max(beta(s), mp.mpf(1))))))
ws = []
for i in range(ni):
    for j in range(ni):
        ri = -mp.mpf('0.97') + mp.mpf('1.94') * i / (ni - 1)
        rr = -mp.mpf('0.97') + mp.mpf('1.94') * j / (ni - 1)
        s = (A12 + ri * pm) / 2
        w0 = gradH(*t_of_sD(s, rr * mp.acosh(max(beta(s), mp.mpf(1)))))
        tot = mp.mpf(0)
        for k in range(len(Gam)):
            a, c = Gam[k], Gam[(k + 1) % len(Gam)]
            da = mp.arg(mp.mpc(c[0] - w0[0], c[1] - w0[1]) / mp.mpc(a[0] - w0[0], a[1] - w0[1]))
            tot += da
        ws.append(tot / (2 * mp.pi))
print('shape (%s, %s) deg: 40-digit winding numbers at %d interior images: min %s, max %s'
      % (sys.argv[1], sys.argv[2], len(ws), mp.nstr(min(ws, key=abs), 6), mp.nstr(max(ws, key=abs), 6)))
