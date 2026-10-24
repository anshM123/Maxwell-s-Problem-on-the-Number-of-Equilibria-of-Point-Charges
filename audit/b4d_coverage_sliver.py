"""AUDIT_2 item 4c: the main chart of code/mixed_cert_all.py clamps the top edge of A12 to the float math.pi, which is
1.22e-16 BELOW pi, so shapes with 0 < A3 < 1.22e-16 and A1, A2 < pi - r0 are not covered by any chart.  This replays the
driver's branch-and-bound on the main chart (project's mixed_cert.I_iv = object under test, imported read-only) with
the top edge widened to ia.up(pi) > pi, to see whether the gap closes.
usage: python b4d_coverage_sliver.py [e0] [r0]"""
import math
import os
import sys

import numpy as np

CODE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'code'))
sys.path.insert(0, CODE)
import ia          # noqa: E402
import mixed_cert  # noqa: E402

e0 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.2
r0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
top = float(ia.up(np.pi))
print('ia.up(pi) = %r > pi (float pi = %r)' % (top, math.pi))
domain, grid = [(e0, top), (0.0, 1.0), (-1.0, 1.0)], (24, 16, 16)
edges = [np.linspace(lo, hi, n + 1) for (lo, hi), n in zip(domain, grid)]
lo = np.array(np.meshgrid(*[e[:-1] for e in edges], indexing='ij')).reshape(3, -1)
hi = np.array(np.meshgrid(*[e[1:] for e in edges], indexing='ij')).reshape(3, -1)
for k, (a, b) in enumerate(domain):
    lo[k][lo[k] <= a] = a
    hi[k][hi[k] >= b] = b
scale = np.array([b - a for a, b in domain])[:, None]
ncert, ndrop, minv = 0, 0, np.inf
for rnd in range(60):
    A1min, A2min = lo[1] * lo[0], (1.0 - hi[1]) * lo[0]
    d = (A1min >= math.pi - r0 + 1e-15) | (A2min >= math.pi - r0 + 1e-15)
    ndrop += d.sum()
    lo, hi = lo[:, ~d], hi[:, ~d]
    if lo.shape[1] == 0:
        break
    with np.errstate(all='ignore'):
        res = mixed_cert.I_iv((lo[0], hi[0]), (lo[1], hi[1]), (lo[2], hi[2]))
    ok = (res[0] > 0) & np.isfinite(res[0])
    if ok.any():
        minv = min(minv, res[0][ok].min())
    ncert += ok.sum()
    lo, hi = lo[:, ~ok], hi[:, ~ok]
    if lo.shape[1] == 0:
        break
    w = (hi - lo) / scale
    kk = np.argmax(w, axis=0)
    m = 0.5 * (lo + hi)
    idx = np.arange(lo.shape[1])
    lo2, hi2 = lo.copy(), hi.copy()
    hi[kk, idx] = m[kk, idx]
    lo2[kk, idx] = m[kk, idx]
    lo, hi = np.hstack([lo, lo2]), np.hstack([hi, hi2])
print('main chart with A12 in [%g, up(pi)]: certified %d, dropped %d, FAILED %d, min lower bound %.3e (%d rounds)'
      % (e0, ncert, ndrop, lo.shape[1], minv, rnd + 1))
