"""Coverage patch for the alpha sweep (AUDIT_2): mixed_certF.py runs the main chart on A1 + A2 in [e0, math.pi], and
math.pi < pi, so triangles with 0 < A3 < pi - math.pi (about 1.2e-16) are not in its main domain. This run certifies the
main chart on the strip A1 + A2 in [pi - d, up(pi)] for the whole F-range of alpha in [a_lo, a_hi], with the same corner
drop rule (r0) as the sweep. Together with the sweep this covers A1 + A2 in [e0, pi).
usage: python mixed_certF_sliver.py [a_lo] [a_hi] [r0] [d]
"""
import math
import sys

import numpy as np

import ia
from mixed_certF import I_main, _F_of_alpha, bnb

if __name__ == '__main__':
    a_lo = float(sys.argv[1]) if len(sys.argv) > 1 else 0.25
    a_hi = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
    r0 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05
    d = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-6
    F_dom = (ia.dn(_F_of_alpha(a_lo)), ia.up(_F_of_alpha(a_hi)))
    A_dom = (math.pi - d, float(ia.up(np.pi)))
    print('alpha in [%g, %g] -> F in [%.6f, %.6f];  A1 + A2 in [%r, %r]' % (a_lo, a_hi, F_dom[0], F_dom[1], *A_dom))

    def drop_main(lo, hi):
        return ((lo[1] * lo[0] >= math.pi - r0 + 1e-15) | ((1.0 - hi[1]) * lo[0] >= math.pi - r0 + 1e-15))
    lo, hi = bnb('main sliver', I_main, [A_dom, (0.0, 1.0), (-1.0, 1.0), F_dom], (1, 16, 16, 16), drop=drop_main)
    print('TOTAL failed: %d  ->  %s' % (lo.shape[1], 'SLIVER CERTIFIED' if lo.shape[1] == 0 else 'not complete'))
