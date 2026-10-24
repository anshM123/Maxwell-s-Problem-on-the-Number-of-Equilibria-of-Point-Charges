"""AUDIT_3, item 2 (robustness of the libm assumption): re-certify two alpha-slices with much wider outward rounding.

The alpha-range certificate assumes that the math library is accurate to 2 ulps (ia.py widens IEEE operations by 2 ulps
and library transcendentals by 4 ulps). Here ia.dn/ia.up are replaced by 4-ulp steps and ia.dn4/ia.up4 by 16-ulp steps
(so a library function may be off by up to 14 ulps), and the four charts of mixed_certF.py are run through
mixed_certF.bnb exactly as mixed_certF.py's main block does (e0 = 0.4, r0 = 0.05, same grids, same drop rule) on the
slices alpha in [0.25, 0.3] and [1.998046875, 2.0] (the last slice, the hardest). If both still certify, the
certificate's margins are far from the rounding level there.
Output: audit/c11_wide_rounding.out"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'code'))
import ia  # noqa: E402  (object under test, patched)
import mixed_certF as M  # noqa: E402  (object under test)


def steps(x, k, d):
    for _ in range(k):
        x = np.nextafter(x, d)
    return x


ia.dn = lambda x: steps(x, 4, -np.inf)
ia.up = lambda x: steps(x, 4, np.inf)
ia.dn4 = lambda x: steps(x, 16, -np.inf)
ia.up4 = lambda x: steps(x, 16, np.inf)
ia.PI = (steps(np.float64(math.pi), 4, -np.inf), steps(np.float64(math.pi), 4, np.inf))
e0, r0 = 0.4, 0.05
for a_lo, a_hi in ((0.25, 0.3), (1.998046875, 2.0)):
    F_dom = (ia.dn(M._F_of_alpha(a_lo)), ia.up(M._F_of_alpha(a_hi)))
    print('alpha [%r, %r], F-box [%r, %r], rounding: 4 ulps per operation, 16 ulps per library function'
          % (a_lo, a_hi, float(F_dom[0]), float(F_dom[1])), flush=True)

    def drop_main(lo, hi):
        return (lo[1] * lo[0] >= math.pi - r0 + 1e-15) | ((1.0 - hi[1]) * lo[0] >= math.pi - r0 + 1e-15)
    fails = 0
    lo, hi = M.bnb('main', M.I_main, [(e0, math.pi), (0.0, 1.0), (-1.0, 1.0), F_dom], (24, 16, 16, 2), drop=drop_main)
    fails += lo.shape[1]
    lo, hi = M.bnb('band', M.I4_band, [(0.0, e0), (0.0, 1.0), (-1.0, 1.0), F_dom], (8, 16, 16, 2))
    fails += lo.shape[1]
    for sig in (+1, -1):
        lo, hi = M.bnb('corner %+d' % sig, lambda r, th, rh, F: M.J1n_corner(r, th, rh, F, sig),
                       [(0.0, r0), (0.0, 1.0), (-1.0, 1.0), F_dom], (8, 8, 16, 2))
        fails += lo.shape[1]
    print('   TOTAL failed: %d -> %s' % (fails, 'CERTIFIED with wide rounding' if fails == 0 else 'not complete'), flush=True)
