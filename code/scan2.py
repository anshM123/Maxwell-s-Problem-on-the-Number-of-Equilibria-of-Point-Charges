"""Zoomed scan of the region where planar minima exist (Coulomb unless alpha is given): triangle angles in
[lo, hi] degrees and charge log-ratios in [-L, L].  Reports the (minima, saddles) histogram, parity failures, and every
configuration with two or more minima.

usage: python scan2.py N lo_deg hi_deg L [alpha] [seed]"""
import math
import sys
import time

import numpy as np

from scan import classify, triangle_from_angles

n = int(sys.argv[1])
lo, hi = math.radians(float(sys.argv[2])), math.radians(float(sys.argv[3]))
L = float(sys.argv[4])
alpha = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5
rng = np.random.default_rng(int(sys.argv[6]) if len(sys.argv) > 6 else 0)
hist, bad, t0, found = {}, 0, time.time(), []
k = 0
while k < n:
    A, B = rng.uniform(lo, hi, 2)
    C = math.pi - A - B
    if not lo <= C <= hi:
        continue
    k += 1
    P = triangle_from_angles(A, B)
    q = np.exp(rng.uniform(-L, L, 3))
    m, s, pts = classify(P, q, alpha)
    hist[(m, s)] = hist.get((m, s), 0) + 1
    if s - m != 2:
        bad += 1
    if m >= 2:
        found.append((A, B, q.tolist(), m, s))
        print('*** minima %d saddles %d  A=%.12f B=%.12f q=%s  pts=%s' % (
            m, s, A, B, q.tolist(), [(np.round(x, 6).tolist(), kk) for x, kk in pts]), flush=True)
print('alpha %.4g angles [%s, %s] deg, L=%.3f: %d configs in %.0fs; histogram %s; parity failures %d; m0>=2: %d'
      % (alpha, sys.argv[2], sys.argv[3], L, n, time.time() - t0, dict(sorted(hist.items())), bad, len(found)),
      flush=True)
