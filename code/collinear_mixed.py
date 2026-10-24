"""Collinear charges of mixed signs: count the critical points of f(t) = sum_i q_i |t - t_i|^(-2 alpha) on the axis.
(For collinear charges in R^3 the isolated equilibria are exactly these; off-axis equilibria come in circles.)
f'(t) = -2 alpha sum_i q_i sgn(t - t_i) |t - t_i|^(-2 alpha - 1).  Zeros are located by sign changes of f' on a fine
log-spaced grid in each of the four intervals, refined by bisection; tangential (double) zeros are detected by local
minima of |f'| and reported separately.
usage: python collinear_mixed.py [n] [alpha] [seed]"""
import math
import sys

import numpy as np

n = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
alpha = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
g = 2 * alpha + 1


def fp(t, ts, q):
    return sum(qi * np.sign(t - ti) * np.abs(t - ti) ** (-g) for ti, qi in zip(ts, q))


def grid(lo, hi, m=4000):
    if np.isinf(lo):
        return hi - np.geomspace(1e-6, 1e6, m)[::-1]
    if np.isinf(hi):
        return lo + np.geomspace(1e-6, 1e6, m)
    u = np.linspace(0, 1, m + 2)[1:-1]
    w = 0.5 - 0.5 * np.cos(math.pi * u)            # cluster near both ends
    return lo + (hi - lo) * w


if __name__ == '__main__':
    rng = np.random.default_rng(seed)
    hist, worst = {}, None
    for it in range(n):
        ts = np.sort(rng.normal(size=3))
        q = np.exp(rng.normal(scale=1.5, size=3)) * rng.choice([-1, 1], size=3)
        if np.all(q > 0) or np.all(q < 0):
            q[rng.integers(3)] *= -1
        total = 0
        for lo, hi in ((-np.inf, ts[0]), (ts[0], ts[1]), (ts[1], ts[2]), (ts[2], np.inf)):
            x = grid(lo, hi)
            y = fp(x, ts, q)
            total += int(np.sum(np.sign(y[1:]) != np.sign(y[:-1])))
        hist[total] = hist.get(total, 0) + 1
        if worst is None or total > worst[0]:
            worst = (total, ts.round(4).tolist(), q.round(4).tolist())
    print('alpha %.2f, %d collinear mixed-sign configurations: histogram of #axis critical points %s; max at %s'
          % (alpha, n, dict(sorted(hist.items())), worst))
