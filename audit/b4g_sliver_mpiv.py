"""AUDIT_2 item 4c (rigorous): close the coverage sliver A12 in (fl(pi), pi] in mp.iv.  Every certified main-chart box
whose top A12 edge is the float math.pi is extended to the float just above pi (> true pi) and re-certified with the
auditor's mp.iv evaluator of b4e_mpiv_recert.py (local bisection allowed).  Together with b4e this rigorously covers
A12 in [e0, pi] for the main chart."""
import math
import os
import sys

import mpmath as mp
import numpy as np

src = open('b4e_mpiv_recert.py', encoding='utf-8').read()
ns = {'__file__': os.path.abspath('b4e_mpiv_recert.py'), '__name__': 'helpers'}
code = src.split('dat = np.load')[0]
code = code.replace('prec = int(sys.argv[1]) if len(sys.argv) > 1 else 64', 'prec = 80')
code = code.replace('stride = int(sys.argv[2]) if len(sys.argv) > 2 else 1', 'stride = 1')
exec(code, ns)
dat = np.load('b4a_main_certified_boxes.npz')
lo, hi = dat['lo'], dat['hi']
top = hi[0] == math.pi
up_pi = float(np.nextafter(np.nextafter(math.pi, np.inf), np.inf))
print('boxes with top edge = float pi: %d; extending to %r (> pi: %s)' % (top.sum(), up_pi, mp.mpf(up_pi) > mp.pi))
fails, leaves = 0, 0
for j in np.nonzero(top)[0]:
    box = ((float(lo[0, j]), up_pi), (float(lo[1, j]), float(hi[1, j])), (float(lo[2, j]), float(hi[2, j])))
    ok, nl, mn = ns['certify'](box)
    leaves += nl
    fails += not ok
print('mp.iv re-certification of the extended top layer: NOT certified %d (leaf boxes %d)' % (fails, leaves))
