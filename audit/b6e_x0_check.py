"""AUDIT_2 item 6: the 11 'x0 not recovered' entries of b6a_direct_search.out (extreme band shapes).  For each, the
50-digit critical points are recomputed; the unique extremum x_min must satisfy q proportional to Q(x_min)
= (lambda_i r_i^3)_i (then it IS the forced minimum, and the float64 location test was merely too strict)."""
import ast
import math
import re

import mpmath as mp
import numpy as np
from scipy.optimize import root

mp.mp.dps = 50
lines = open('b6a_direct_search.out', encoding='utf-8').read().splitlines()
fails = [l for l in lines if l.strip().startswith('failure:')]
susp = [l for l in lines if 'float64 suspicious' in l]
done = 0
for fl in fails:
    qs = re.search(r'q (\[[^\]]*\])', fl).group(1)
    q = ast.literal_eval(qs)
    src = [l for l in susp if qs in l]
    if not src:
        continue
    m = re.search(r'angles ([0-9.e-]+), ([0-9.e-]+), q', src[0])
    A1, A2 = float(m.group(1)), float(m.group(2))
    A3 = math.pi - A1 - A2
    L13 = math.sin(A2) / math.sin(A3)
    P = np.array([[0.0, 0.0], [1.0, 0.0], [L13 * math.cos(A1), L13 * math.sin(A1)]])
    qn = np.array(q)

    def G(x):
        r = np.linalg.norm(x[None] - P, axis=1)
        w = qn * r ** -3.0
        return (w[:, None] * (x[None] - P)).sum(0) / (np.abs(w) * r).sum()

    cen = P.mean(0)
    cands = []
    for rr in np.geomspace(1e-4, 1e4, 60):
        for th in np.linspace(0, 2 * math.pi, 48, endpoint=False):
            sol = root(G, cen + rr * np.array([math.cos(th), math.sin(th)]), method='hybr', options={'xtol': 1e-15})
            if np.all(np.isfinite(sol.x)) and np.linalg.norm(G(sol.x)) < 1e-9 and all(np.linalg.norm(sol.x - c) > 1e-9 for c in cands):
                cands.append(sol.x)
    Pm = [[mp.mpf(float(v)) for v in p] for p in P]
    qm = [mp.mpf(v) for v in q]
    D2 = (Pm[1][0] - Pm[0][0]) * (Pm[2][1] - Pm[0][1]) - (Pm[2][0] - Pm[0][0]) * (Pm[1][1] - Pm[0][1])
    best = None
    for x in cands:
        xm = [mp.mpf(float(x[0])), mp.mpf(float(x[1]))]
        for it in range(80):
            g = [mp.mpf(0)] * 2
            H = mp.matrix(2, 2)
            for p, qi in zip(Pm, qm):
                d = [xm[0] - p[0], xm[1] - p[1]]
                r2 = d[0] ** 2 + d[1] ** 2
                r = mp.sqrt(r2)
                for a in range(2):
                    g[a] += -qi * d[a] / r ** 3
                    for c in range(2):
                        H[a, c] += -qi / r ** 3 * ((1 if a == c else 0) - 3 * d[a] * d[c] / r2)
            try:
                st = mp.lu_solve(H, mp.matrix(g))
            except ZeroDivisionError:
                break
            xm = [xm[0] - st[0], xm[1] - st[1]]
            if mp.sqrt(st[0] ** 2 + st[1] ** 2) < mp.mpf('1e-40'):
                ev = mp.eigsy(H)[0]
                if ev[0] * ev[1] > 0:
                    best = xm
                break
        if best is not None:
            break
    lam = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        lam.append(((Pm[j][0] - best[0]) * (Pm[k][1] - best[1]) - (Pm[j][1] - best[1]) * (Pm[k][0] - best[0])) / D2)
    Q = [lam[i] * mp.sqrt((best[0] - Pm[i][0]) ** 2 + (best[1] - Pm[i][1]) ** 2) ** 3 for i in range(3)]
    ratios = [Q[i] / qm[i] for i in range(3)]
    spread = (max(ratios) - min(ratios)) / abs(ratios[0])
    done += 1
    print('angles (%.4g, %.4g): extremum at (%s, %s), |x_min| %s; relative spread of Q(x_min)_i / q_i = %s'
          % (A1, A2, mp.nstr(best[0], 8), mp.nstr(best[1], 8), mp.nstr(mp.sqrt(best[0] ** 2 + best[1] ** 2), 4),
             mp.nstr(spread, 3)))
print('checked %d of %d entries' % (done, len(fails)))
