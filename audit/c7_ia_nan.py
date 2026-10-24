"""AUDIT_3, item 2: non-finite handling in code/ia.py and whether it matters for mixed_certF.py.

Part A (a property of ia.py): sinc, asinc and asinhc compute f(y)/y through ia._ratio, which returns 1.0 wherever
y > 0 is False -- also for y = NaN. A NaN input interval therefore produces a finite, wrong enclosure (about [1, 1]),
instead of NaN. Every other ia function propagates NaN to the lower end, or returns a valid global bound
(-inf for log, 0 for exp, [-1, 1] for wide sin).
Part B (does it matter?): the band chart (the only chart that calls sinc/asinc/asinhc) is run through
mixed_certF.bnb on the F-boxes of two certified slices (alpha in [0.25, 0.3] and [1.998046875, 2.0], e0 = 0.4) with
sinc/asinc/asinhc wrapped to count non-finite input intervals; any NaN or inf input would be reported. Every input of
these calls is built from finite box coordinates by finite operations, so none is expected.
Also, every chart's lower bounds are checked: a certified box must have a finite lower bound (bnb requires
isfinite(lo) & lo > 0) -- this is code inspection, recorded here.
Output: audit/c7_ia_nan.out"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'code'))
import ia  # noqa: E402  (object under test)
import mixed_certF  # noqa: E402  (object under test)

nan = np.array([np.nan])
print('Part A: NaN input intervals')
for name in ('sinc', 'asinc', 'asinhc', 'sin', 'cos', 'log', 'exp', 'sinh', 'cosh', 'sqrt', 'sq', 'arccos', 'arcsin',
             'arcsinh', 'arctan'):
    f = getattr(ia, name)
    try:
        with np.errstate(all='ignore'):
            lo, hi = f((nan, nan))
        print('   ia.%-8s(NaN) -> [%r, %r]' % (name, float(np.asarray(lo).ravel()[0]), float(np.asarray(hi).ravel()[0])))
    except Exception as ex:  # noqa: BLE001
        print('   ia.%-8s(NaN) -> raises %s' % (name, type(ex).__name__))

print('Part B: non-finite inputs reaching sinc/asinc/asinhc during band-chart runs')
counts = {}
orig = {n: getattr(ia, n) for n in ('sinc', 'asinc', 'asinhc')}


def wrap(n):
    def g(x):
        bad = ~(np.isfinite(x[0]) & np.isfinite(x[1]))
        counts[n] = counts.get(n, 0) + int(bad.sum())
        counts[n + '_calls'] = counts.get(n + '_calls', 0) + int(np.size(x[0]))
        return orig[n](x)
    return g


for n in orig:
    setattr(ia, n, wrap(n))
for a_lo, a_hi in ((0.25, 0.3), (1.998046875, 2.0)):
    counts.clear()
    Fd = (ia.dn(mixed_certF._F_of_alpha(a_lo)), ia.up(mixed_certF._F_of_alpha(a_hi)))
    lo, hi = mixed_certF.bnb('band', mixed_certF.I4_band, [(0.0, 0.4), (0.0, 1.0), (-1.0, 1.0), Fd], (8, 16, 16, 2))
    print('   alpha [%r, %r]: failed %d; non-finite inputs: sinc %d of %d, asinc %d of %d, asinhc %d of %d'
          % (a_lo, a_hi, lo.shape[1], counts.get('sinc', 0), counts.get('sinc_calls', 0), counts.get('asinc', 0),
             counts.get('asinc_calls', 0), counts.get('asinhc', 0), counts.get('asinhc_calls', 0)))
for n, f in orig.items():
    setattr(ia, n, f)
