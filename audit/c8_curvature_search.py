"""AUDIT_3, items 7 and 8: adversarial search for negative boundary curvature of M' at alpha = 1/4, 1, 2, and an
independent spot check of the large-exponent remark (paper Remark 9.1, MIXED.md section 8).

All values come from the auditor's reference audit/c0_ref.py (curvature numerator of {E < K} from E's derivatives,
boundary point by root finding on E), not from the certificate's closed forms.
(1) For each alpha in {1/4, 1, 2}: the scale-free I/(A1+A2)^4 (vertex 3, upper arc; the lower arc is the upper arc of
    the mirrored shape, which is sampled too) at random (shape, rho) in five families -- generic, flat (A3 down to
    1e-8), thin (A1+A2 down to 1e-5), lopsided (A2/A1 or A1/A2 down to 1e-6) -- and J1/U2^2 (vertex 1, both arcs) on
    corner shapes (A1 = pi - r, r down to 1e-6); then a local random refinement around the five smallest values.
(2) Large exponents: along the degenerate corner A1 + A2 = 1e-3, A2/A1 = 1e-5 (as in convexity_threshold.py part 2),
    min over rho of I/(A1+A2)^4 at alpha = 5.6, 5.7, 5.8 (grid + golden refinement), compared with
    cert/convexity_threshold.out; and the alpha = 10 example triangle (10.6, 29.3, 140.1 degrees): min of I over both
    arcs.
Output: audit/c8_curvature_search.out"""
import math
import os
import random
import sys
import time

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import c0_ref as R  # noqa: E402

rnd = random.Random(8)


def Fof(alpha):
    b = 2 * mp.mpf(alpha) + 1
    return (b + 1) ** 2 / (2 * b) - 1


def score(A12, lam, rho, F):
    v = R.I_ref_main(A12, lam, rho, F)
    return None if v is None else v / mp.mpf(A12) ** 4


def sample(fam):
    if fam == 'generic':
        return mp.mpf(rnd.uniform(0.4, math.pi - 0.01)), mp.mpf(rnd.uniform(0.01, 0.99))
    if fam == 'flat':
        return mp.pi - mp.mpf(10) ** rnd.uniform(-8, -1), mp.mpf(rnd.uniform(0.02, 0.98))
    if fam == 'thin':
        u = mp.mpf(10) ** rnd.uniform(-6, -0.3)
        return mp.mpf(10) ** rnd.uniform(-5, -0.4), (u if rnd.random() < 0.5 else 1 - u)
    u = mp.mpf(10) ** rnd.uniform(-6, -1)      # lopsided
    return mp.mpf(rnd.uniform(0.4, 3.0)), (u if rnd.random() < 0.5 else 1 - u)


def rho_sample():
    return mp.mpf(rnd.uniform(-1, 1)) if rnd.random() < 0.6 else mp.mpf(rnd.choice((-1, 1)) * (1 - 10 ** rnd.uniform(-6, -1)))


def part1():
    lines = []
    for alpha in (mp.mpf(1) / 4, mp.mpf(1), mp.mpf(2)):
        mp.mp.dps = 40
        F = Fof(alpha)
        t0 = time.time()
        res = []
        neg = 0
        fammin = {}
        for fam, n in (('generic', 300), ('flat', 250), ('thin', 250), ('lopsided', 200)):
            for _ in range(n):
                A12, lam = sample(fam)
                rho = rho_sample()
                v = score(A12, lam, rho, F)
                if v is None:
                    continue
                neg += v <= 0
                res.append((v, fam, A12, lam, rho))
                fammin[fam] = min(fammin.get(fam, mp.inf), v)
        res.sort(key=lambda z: z[0])
        # local refinement around the five smallest
        best = list(res[:5])
        for k in range(len(best)):
            v, fam, A12, lam, rho = best[k]
            for it in range(50):
                sc = mp.mpf(10) ** rnd.uniform(-4, -1)
                A12n = A12 * (1 + sc * mp.mpf(rnd.uniform(-1, 1)))
                if fam == 'flat':
                    A3n = (mp.pi - A12) * (1 + sc * 10 * mp.mpf(rnd.uniform(-0.9, 0.9)))
                    A12n = mp.pi - A3n
                laml = min(max(lam * (1 + sc * mp.mpf(rnd.uniform(-1, 1))), mp.mpf(10) ** -12), 1 - mp.mpf(10) ** -12)
                rhon = min(max(rho + sc * mp.mpf(rnd.uniform(-1, 1)), mp.mpf(-1)), mp.mpf(1))
                if not (0 < A12n < mp.pi):
                    continue
                vn = score(A12n, laml, rhon, F)
                if vn is not None and vn < v:
                    v, A12, lam, rho = vn, A12n, laml, rhon
                    neg += v <= 0
            best[k] = (v, fam, A12, lam, rho)
        best.sort(key=lambda z: z[0])
        v, fam, A12, lam, rho = best[0]
        # corner chart (vertex 1), both arcs
        cmin, cneg = mp.inf, 0
        for _ in range(150):
            r = mp.mpf(10) ** rnd.uniform(-6, math.log10(0.05))
            th = mp.mpf(rnd.uniform(0, 1)) if rnd.random() < 0.5 else mp.mpf(10) ** rnd.uniform(-6, -1)
            if rnd.random() < 0.5:
                th = 1 - th
            rho1 = rho_sample()
            for sig in (1, -1):
                w = R.J1n_ref_corner(r, th, rho1, F, sig)
                if w is not None:
                    cmin = min(cmin, w)
                    cneg += w <= 0
        A1d, A2d = float(lam * A12) * 180 / math.pi, float((1 - lam) * A12) * 180 / math.pi
        lines.append('alpha %-4s (F = %s): %d samples + refinement, %.0fs; non-positive values: %d; minimum of I/(A1+A2)^4 '
                     '= %s (%s family) at angles (%.4f, %.4f, %.4f) deg, rho %s; family minima: %s; corner chart: '
                     'min J1/U2^2 = %s over 300 points, non-positive %d'
                     % (mp.nstr(alpha, 3), mp.nstr(F, 6), len(res), time.time() - t0, neg, mp.nstr(v, 6), fam, A1d, A2d,
                        180 - A1d - A2d, mp.nstr(rho, 6),
                        ', '.join('%s %s' % (k, mp.nstr(x, 4)) for k, x in fammin.items()), mp.nstr(cmin, 5), cneg))
        print(lines[-1], flush=True)
    return lines


def part2():
    lines = []
    mp.mp.dps = 50
    ref_table = {'5.6': '+1.0709e-02', '5.7': '+5.2544e-03', '5.8': '-3.4430e-04'}
    for a in ('5.6', '5.7', '5.8'):
        F = Fof(mp.mpf(a))
        e, lam = mp.mpf('1e-3'), 1 - mp.mpf('1e-5')
        grid = [(score(e, lam, mp.mpf('0.80') + mp.mpf('0.005') * k, F), mp.mpf('0.80') + mp.mpf('0.005') * k)
                for k in range(29)]
        gmin = min(grid, key=lambda z: z[0])
        f = lambda r: score(e, lam, r, F)
        rho_star = R._golden_min(f, gmin[1] - mp.mpf('0.005'), gmin[1] + mp.mpf('0.005'), 80)
        lines.append('alpha %s: e = 1e-3, A2/A1 = 1e-5: grid minimum over rho in [0.80, 0.94] = %s (convexity_threshold.out:'
                     ' %s); refined minimum %s at rho %s'
                     % (a, mp.nstr(gmin[0], 5), ref_table[a], mp.nstr(f(rho_star), 5), mp.nstr(rho_star, 5)))
        print(lines[-1], flush=True)
    # alpha = 10 example
    F = Fof(10)
    A1, A2 = mp.radians(mp.mpf('10.6')), mp.radians(mp.mpf('29.3'))
    A12 = A1 + A2
    worst = (mp.inf, None)
    for arc, lam in (('upper', A1 / A12), ('lower (mirror)', A2 / A12)):
        for k in range(201):
            rho = -1 + mp.mpf(k) / 100
            v = R.I_ref_main(A12, lam, rho, F)
            if v is not None and v < worst[0]:
                worst = (v, (arc, rho))
    lines.append('alpha 10, triangle (10.6, 29.3, 140.1) deg: min of I over both arcs (201 points each) = %s at %s arc, '
                 'rho %s' % (mp.nstr(worst[0], 5), worst[1][0], mp.nstr(worst[1][1], 4)))
    print(lines[-1], flush=True)
    return lines


if __name__ == '__main__':
    part1()
    part2()
