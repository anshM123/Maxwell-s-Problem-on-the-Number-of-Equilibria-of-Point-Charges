"""End-to-end test of the theorem's exact predictions against a direct critical-point search.

Pick a triangle and a point x in it, and set q = (lambda_i r_i^beta) (so x is critical). With E = e2(W) at theta(x) and
k = (b + 1)^2 / b:
  E < k  ->  x is a minimum, it is the ONLY minimum, and there are exactly 4 critical points;
  E > k  ->  x is a saddle, at most one minimum, and #critical points = 2 m0 + 2 (so 2 or 4).
Points are sampled uniformly and also with E close to k (near the fold), where the search is hardest.
usage: python predict_test.py [n] [alpha] [seed]"""
import math
import sys

import numpy as np

from key_identities import W_of, e2
from scan import classify
from separable import tri, bary, angles_at

n = int(sys.argv[1]) if len(sys.argv) > 1 else 400
alpha = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
beta = 2 * alpha + 2
b = beta - 1
k = (b + 1) ** 2 / b
rng = np.random.default_rng(seed)
stats = {'in': 0, 'out': 0, 'fail': 0}
for it in range(n):
    while True:
        A = rng.dirichlet(np.ones(3)) * math.pi
        if A.min() > math.radians(4):
            break
    P = tri(A)
    near = it % 2 == 1
    for _ in range(2000):
        x = rng.dirichlet(np.ones(3)) @ P
        E = e2(W_of(angles_at(P, x), A))
        if not near or abs(E / k - 1) < 0.03:
            break
    lam, r = bary(P, x), np.linalg.norm(x - P, axis=1)
    q = lam * r ** beta
    q /= q.sum()
    m, s, pts = classify(P, q, alpha)
    here = [kd for (y, kd) in pts if np.hypot(*(y - x)) < 1e-6]
    if E < k:
        stats['in'] += 1
        ok = m == 1 and s == 3 and here == ['min']
    else:
        stats['out'] += 1
        ok = m <= 1 and s == m + 2 and here == ['saddle']
    if not ok:
        stats['fail'] += 1
        print('MISMATCH angles %s x %s E/k %.4f: m=%d s=%d at x: %s'
              % (np.round(np.degrees(A), 2).tolist(), np.round(x, 5).tolist(), E / k, m, s, here), flush=True)
print('alpha %.2f: %d points with E < k (predicted 1 min + 3 saddles), %d with E > k; mismatches: %d'
      % (alpha, stats['in'], stats['out'], stats['fail']))
