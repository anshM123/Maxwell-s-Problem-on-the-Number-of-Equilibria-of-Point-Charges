"""AUDIT_2 item 4a: soundness tests of code/ia.py (the object under test; imported read-only from code).

For every operation, random input intervals (narrow, medium, wide; including edge cases: intervals containing 0,
containing extrema of sin/cos, near the singular ends of arccos/arcsin/asinc) are evaluated by ia.py and the result is
compared with the TRUE range computed by the auditor in 50-digit mpmath (endpoint values plus every interior critical
point, plus dense sampling).  A failure is a true value outside the returned interval.
Also measured: the actual ulp error of numpy's transcendentals on this machine (ia.py assumes <= 2 ulps and widens by 4).
usage: python b4b_ia_tests.py [n]"""
import math
import os
import sys

import mpmath as mp
import numpy as np

CODE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'code'))
sys.path.insert(0, CODE)
import ia  # noqa: E402  (object under test)

mp.mp.dps = 50
rng = np.random.default_rng(99)
n = int(sys.argv[1]) if len(sys.argv) > 1 else 4000


def rand_iv(lo, hi, kind):
    a = rng.uniform(lo, hi)
    w = {'narrow': 1e-12, 'medium': 1e-3, 'wide': 1.0}[kind] * rng.uniform(0, 1) * (hi - lo)
    return a, min(a + w, hi)


def true_range(f, a, b, crit=()):
    """range of f on [a, b] (mpmath): endpoints, given interior critical points, dense sampling"""
    a, b = mp.mpf(a), mp.mpf(b)
    vals = [f(a), f(b)] + [f(c) for c in crit if a <= c <= b]
    for k in range(1, 16):
        vals.append(f(a + (b - a) * k / 16))
    return min(vals), max(vals)


def check(name, res, tr, fails):
    lo, hi = float(res[0]), float(res[1])
    if not (mp.mpf(lo) <= tr[0] and tr[1] <= mp.mpf(hi)):
        fails[name] = fails.get(name, 0) + 1
        if fails[name] <= 3:
            print('   FAIL %s: returned [%r, %r], true [%s, %s]' % (name, lo, hi, mp.nstr(tr[0], 20), mp.nstr(tr[1], 20)))


fails, counts = {}, {}
for k in range(n):
    kind = ['narrow', 'medium', 'wide'][k % 3]
    # binary ops on random signs
    a = rand_iv(-5, 5, kind); b_ = rand_iv(-5, 5, kind)
    X, Y = (np.float64(a[0]), np.float64(a[1])), (np.float64(b_[0]), np.float64(b_[1]))
    Xm = [mp.mpf(a[0]), mp.mpf(a[1])]
    Ym = [mp.mpf(b_[0]), mp.mpf(b_[1])]
    prods = [x * y for x in Xm for y in Ym]
    check('add', ia.add(X, Y), (Xm[0] + Ym[0], Xm[1] + Ym[1]), fails)
    check('sub', ia.sub(X, Y), (Xm[0] - Ym[1], Xm[1] - Ym[0]), fails)
    check('mul', ia.mul(X, Y), (min(prods), max(prods)), fails)
    sq_lo = mp.mpf(0) if Xm[0] <= 0 <= Xm[1] else min(Xm[0] ** 2, Xm[1] ** 2)
    check('sq', ia.sq(X), (sq_lo, max(Xm[0] ** 2, Xm[1] ** 2)), fails)
    if not (Ym[0] <= 0 <= Ym[1]):
        q = [x / y for x in Xm for y in Ym]
        check('div', ia.div(X, Y), (min(q), max(q)), fails)
        check('recip', ia.recip(Y), (1 / Ym[1], 1 / Ym[0]), fails)
    else:
        r = ia.div(X, Y)
        if not (r[0] == -np.inf and r[1] == np.inf):
            fails['div by interval containing 0 not (-inf, inf)'] = fails.get('div by interval containing 0 not (-inf, inf)', 0) + 1
        try:
            ia.recip(Y)
            fails['recip of interval containing 0 did not raise'] = fails.get('recip of interval containing 0 did not raise', 0) + 1
        except ZeroDivisionError:
            pass
    cval = rng.uniform(-3, 3)
    sc = ia.scal(cval, X)
    check('scal', sc, (min(mp.mpf(cval) * Xm[0], mp.mpf(cval) * Xm[1]), max(mp.mpf(cval) * Xm[0], mp.mpf(cval) * Xm[1])), fails)
    # sin / cos on [-10, 10], with extrema
    s_ = rand_iv(-10, 10, kind)
    S = (np.float64(s_[0]), np.float64(s_[1]))
    crit_s = [mp.pi / 2 + j * mp.pi for j in range(-5, 5)]
    check('sin', ia.sin(S), true_range(mp.sin, *s_, crit=crit_s), fails)
    check('cos', ia.cos(S), true_range(mp.cos, *s_, crit=[j * mp.pi for j in range(-4, 5)]), fails)
    # monotone functions
    p_ = rand_iv(1e-8, 20, kind); Pv = (np.float64(p_[0]), np.float64(p_[1]))
    check('log', ia.log(Pv), true_range(mp.log, *p_), fails)
    check('sqrt', ia.sqrt(Pv), true_range(mp.sqrt, *p_), fails)
    e_ = rand_iv(-20, 20, kind); Ev = (np.float64(e_[0]), np.float64(e_[1]))
    check('exp', ia.exp(Ev), true_range(mp.exp, *e_), fails)
    check('sinh', ia.sinh(Ev), true_range(mp.sinh, *e_), fails)
    check('cosh', ia.cosh(Ev), true_range(mp.cosh, *e_, crit=[mp.mpf(0)]), fails)
    check('arcsinh', ia.arcsinh(Ev), true_range(mp.asinh, *e_), fails)
    check('arctan', ia.arctan(Ev), true_range(mp.atan, *e_), fails)
    u_ = rand_iv(-1, 1, kind); Uv = (np.float64(u_[0]), np.float64(u_[1]))
    check('arccos', ia.arccos(Uv), true_range(mp.acos, *u_), fails)
    check('arcsin', ia.arcsin(Uv), true_range(mp.asin, *u_), fails)
    # sinc on |t| < pi, asinc on |t| < 1, asinhc anywhere (even functions)
    z_ = rand_iv(-3.1, 3.1, kind); Zv = (np.float64(z_[0]), np.float64(z_[1]))
    sincf = lambda t: mp.sin(t) / t if t != 0 else mp.mpf(1)
    check('sinc', ia.sinc(Zv), true_range(sincf, *z_, crit=[mp.mpf(0)]), fails)
    w_ = rand_iv(-0.999, 0.999, kind); Wv = (np.float64(w_[0]), np.float64(w_[1]))
    asincf = lambda t: mp.asin(t) / t if t != 0 else mp.mpf(1)
    check('asinc', ia.asinc(Wv), true_range(asincf, *w_, crit=[mp.mpf(0)]), fails)
    h_ = rand_iv(-30, 30, kind); Hv = (np.float64(h_[0]), np.float64(h_[1]))
    asinhcf = lambda t: mp.asinh(t) / t if t != 0 else mp.mpf(1)
    check('asinhc', ia.asinhc(Hv), true_range(asinhcf, *h_, crit=[mp.mpf(0)]), fails)
# edge cases
edge = []
for (a, b_) in [(0.0, 0.0), (-1e-300, 1e-300), (0.0, 1e-20), (-2.0, 0.0)]:
    r = ia.sinc((np.float64(a), np.float64(b_)))
    edge.append(('sinc', a, b_, float(r[0]), float(r[1])))
for (a, b_) in [(np.pi / 2 - 1e-17, np.pi / 2 + 1e-17), (np.pi / 2, np.pi / 2), (-np.pi / 2, -np.pi / 2), (np.pi, np.pi)]:
    r = ia.sin((np.float64(a), np.float64(b_)))
    tr = true_range(mp.sin, a, b_, crit=[mp.pi / 2, -mp.pi / 2])
    check('sin edge', r, tr, fails)
    edge.append(('sin', a, b_, float(r[0]), float(r[1])))
r = ia.cos((np.float64(0.0), np.float64(0.0)))
check('cos edge', r, (mp.mpf(1), mp.mpf(1)), fails)
r = ia.cos((np.float64(np.pi), np.float64(np.pi)))
check('cos edge pi', r, true_range(mp.cos, np.pi, np.pi), fails)
r = ia.log((np.float64(0.0), np.float64(1.0)))
edge.append(('log [0,1]', 0.0, 1.0, float(r[0]), float(r[1])))
r = ia.sqrt((np.float64(-1e-20), np.float64(4.0)))
edge.append(('sqrt [-1e-20, 4]', -1e-20, 4.0, float(r[0]), float(r[1])))
print('ia.py containment tests: %d random cases per operation (narrow/medium/wide); failures: %s' % (n, fails if fails else 'none'))
for e in edge:
    print('   edge case %s on [%r, %r] -> [%r, %r]' % e)

# libm ulp errors on this machine
def ulp_err(fn_np, fn_mp, xs):
    worst = 0.0
    y = fn_np(xs)
    for x, yv in zip(xs, y):
        t = fn_mp(mp.mpf(float(x)))
        if t == 0:
            continue
        u = np.spacing(np.float64(abs(float(t))))
        worst = max(worst, float(abs(mp.mpf(float(yv)) - t) / u))
    return worst


N = 20000
tests = [('sin', np.sin, mp.sin, rng.uniform(-10, 10, N)), ('cos', np.cos, mp.cos, rng.uniform(-10, 10, N)),
         ('exp', np.exp, mp.exp, rng.uniform(-30, 30, N)), ('log', np.log, mp.log, np.exp(rng.uniform(-30, 5, N))),
         ('sinh', np.sinh, mp.sinh, rng.uniform(-20, 20, N)), ('cosh', np.cosh, mp.cosh, rng.uniform(-20, 20, N)),
         ('arcsinh', np.arcsinh, mp.asinh, rng.uniform(-50, 50, N)), ('arccos', np.arccos, mp.acos, rng.uniform(-1, 1, N)),
         ('arcsin', np.arcsin, mp.asin, rng.uniform(-1, 1, N)), ('arctan', np.arctan, mp.atan, rng.uniform(-50, 50, N)),
         ('sqrt', np.sqrt, mp.sqrt, rng.uniform(0, 50, N))]
print('numpy %s: measured max error in ulps over %d random arguments each (assumption in ia.py: <= 2, widened by 4):'
      % (np.__version__, N))
print('   ' + ', '.join('%s %.2f' % (nm, ulp_err(fnp, fmp, xs)) for nm, fnp, fmp, xs in tests))
