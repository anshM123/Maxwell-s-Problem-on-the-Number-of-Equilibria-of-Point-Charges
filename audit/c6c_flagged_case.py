"""AUDIT_3, item 7: the one configuration of c6_direct_search.py (alpha = 1/4, thin triangle) that stayed flagged.

Details from audit/c6b_flagged_detail.out: A1 = 1.1688e-4, A2 = 8.978e-4 (radians), forced point x = (386.1, -1894.9)
(triangle side p1p2 = 1), (K - E)/(K - 4) = 4.1e-3 (close to the fold D = c), charges with |q1 + q2 + q3|/sum|q| ~ 4e-5.
In c6 the charges were computed in binary64, so the configuration analysed at 50 digits there is not exactly the forced
one. Here the construction is redone at 50 digits: the point x (taken as the exact binary64 values printed by c6b) is made
an exact critical point by q_i = lambda_i(x) r_i(x)^beta computed at 50 digits; then all planar critical points are
searched by 50-digit damped Newton iterations from 480 seeds (log-polar grid around the centroid out to 1e5, a cloud
around x, small circles around the charges), and classified by the Hessian.
Output: audit/c6c_flagged_case.out"""
import mpmath as mp

mp.mp.dps = 50
alpha = mp.mpf(1) / 4
beta = 2 * alpha + 2
A1, A2 = mp.mpf('0.00011687833356738089'), mp.mpf('0.0008978065167066201')
A3 = mp.pi - A1 - A2
L = mp.sin(A2) / mp.sin(A3)
P = [mp.matrix([0, 0]), mp.matrix([1, 0]), mp.matrix([L * mp.cos(A1), L * mp.sin(A1)])]
x0 = mp.matrix([mp.mpf('386.115368542459'), mp.mpf('-1894.9418221908613')])


def cross(u, v):
    return u[0] * v[1] - u[1] * v[0]


D2 = cross(P[1] - P[0], P[2] - P[0])
lam = [cross(P[(i + 1) % 3] - x0, P[(i + 2) % 3] - x0) / D2 for i in range(3)]
q = [lam[i] * mp.norm(x0 - P[i]) ** beta for i in range(3)]


def gh(x):
    g = mp.matrix([0, 0])
    H = mp.matrix(2, 2)
    for i in range(3):
        d = x - P[i]
        r2 = d[0] ** 2 + d[1] ** 2
        w = q[i] * r2 ** (-beta / 2)
        g += -2 * alpha * w * d
        for j in range(2):
            for k in range(2):
                H[j, k] += -2 * alpha * w * ((1 if j == k else 0) - beta * d[j] * d[k] / r2)
    return g, H


def newton(x):
    for _ in range(200):
        g, H = gh(x)
        try:
            dx = mp.lu_solve(H, g)
        except ZeroDivisionError:
            return None
        dist = min(mp.norm(x - p) for p in P)
        n = mp.norm(dx)
        lim = max(mp.mpf('0.3') * dist, mp.mpf('0.3') * mp.norm(x) * 0.5)
        if n > lim:
            dx = dx * (lim / n)
        x = x - dx
        if n < mp.mpf(10) ** -40 * (1 + mp.norm(x)):
            break
    g, H = gh(x)
    sc = sum(abs(q[i]) * mp.norm(x - P[i]) ** (-2 * alpha - 1) for i in range(3))
    if mp.norm(g) < mp.mpf(10) ** -38 * sc:
        return x
    return None


seeds = []
c = (P[0] + P[1] + P[2]) / 3
for k in range(30):
    R = mp.mpf(10) ** (-1 + 6 * mp.mpf(k) / 29)
    for j in range(12):
        th = 2 * mp.pi * j / 12
        seeds.append(c + R * mp.matrix([mp.cos(th), mp.sin(th)]))
for j in range(60):
    th = 2 * mp.pi * j / 60
    for R in (mp.mpf(1), mp.mpf(30), mp.mpf(300)):
        seeds.append(x0 + R * mp.matrix([mp.cos(th), mp.sin(th)]))
for p in P:
    for j in range(8):
        th = 2 * mp.pi * j / 8
        seeds.append(p + mp.mpf('1e-4') * mp.matrix([mp.cos(th), mp.sin(th)]))
found = []
for s in seeds:
    x = newton(s)
    if x is None:
        continue
    if not any(mp.norm(x - y) < mp.mpf(10) ** -20 * (1 + mp.norm(y)) for y, _ in found):
        g, H = gh(x)
        ev = mp.eigsy(H)[0]
        kind = 'min' if ev[0] > 0 else 'max' if ev[1] < 0 else 'saddle'
        found.append((x, kind))
g0, H0 = gh(x0)
ev0 = mp.eigsy(H0)[0]
print('charges (50-digit construction): %s; Q_tot/sum|q| = %s' % ([mp.nstr(v, 8) for v in q],
                                                                   mp.nstr(sum(q) / sum(abs(v) for v in q), 3)))
print('forced point: |grad V| relative %s; Hessian eigenvalues %s' % (
    mp.nstr(mp.norm(g0) / sum(abs(q[i]) * mp.norm(x0 - P[i]) ** (-2 * alpha - 1) for i in range(3)), 3),
    [mp.nstr(e, 5) for e in ev0]))
print('%d seeds; distinct critical points found: %d' % (len(seeds), len(found)))
for x, k in found:
    print('   (%s, %s)  %s' % (mp.nstr(x[0], 12), mp.nstr(x[1], 12), k))
ne = sum(k in ('min', 'max') for _, k in found)
ns = sum(k == 'saddle' for _, k in found)
print('extrema %d, saddles %d, index sum %d (expected -2 for Q_tot != 0)' % (ne, ns, ne - ns))
