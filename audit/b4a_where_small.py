"""AUDIT_2 item 4 (harness): replay the branch-and-bound of code/mixed_cert_all.py on the main chart, keeping the
certified boxes, and report where the smallest certified lower bounds occur.  The project's modules ia / mixed_cert are
imported ONLY as the objects under test (read from code, not modified).
usage: python b4a_where_small.py [e0] [r0]"""
import math
import os
import sys

import numpy as np

CODE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'code')
sys.path.insert(0, os.path.abspath(CODE))
import mixed_cert  # noqa: E402  (object under test)

e0 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.2
r0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
domain, grid = [(e0, math.pi), (0.0, 1.0), (-1.0, 1.0)], (24, 16, 16)
edges = [np.linspace(lo, hi, n + 1) for (lo, hi), n in zip(domain, grid)]
lo = np.array(np.meshgrid(*[e[:-1] for e in edges], indexing='ij')).reshape(3, -1)
hi = np.array(np.meshgrid(*[e[1:] for e in edges], indexing='ij')).reshape(3, -1)
for k, (a, b) in enumerate(domain):
    lo[k][lo[k] <= a] = a
    hi[k][hi[k] >= b] = b
print('top edge of the A12 range used by the driver: %r (math.pi); true pi exceeds it by %.3e'
      % (hi[0].max(), math.pi - hi[0].max() + 1.2246467991473532e-16))
scale = np.array([b - a for a, b in domain])[:, None]
cert_lo, cert_hi, cert_val = [], [], []
for rnd in range(60):
    A1min, A2min = lo[1] * lo[0], (1.0 - hi[1]) * lo[0]
    d = (A1min >= math.pi - r0 + 1e-15) | (A2min >= math.pi - r0 + 1e-15)
    lo, hi = lo[:, ~d], hi[:, ~d]
    if lo.shape[1] == 0:
        break
    with np.errstate(all='ignore'):
        res = mixed_cert.I_iv((lo[0], hi[0]), (lo[1], hi[1]), (lo[2], hi[2]))
    ok = (res[0] > 0) & np.isfinite(res[0])
    cert_lo.append(lo[:, ok]); cert_hi.append(hi[:, ok]); cert_val.append(res[0][ok])
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
CL, CH, CV = np.hstack(cert_lo), np.hstack(cert_hi), np.concatenate(cert_val)
print('certified boxes %d, failed %d' % (CL.shape[1], lo.shape[1]))
order = np.argsort(CV)
print('20 smallest certified lower bounds and their boxes (A12, lambda, rho):')
for j in order[:20]:
    print('  %.3e  A12 [%.6f, %.6f]  lam [%.5f, %.5f]  rho [%.4f, %.4f]'
          % (CV[j], CL[0, j], CH[0, j], CL[1, j], CH[1, j], CL[2, j], CH[2, j]))
# distribution of small bounds
for thr in (1e-8, 1e-6, 1e-4, 1e-2):
    sel = CV < thr
    if sel.any():
        print('lower bound < %.0e: %d boxes; A12 range [%.4f, %.4f], lambda range [%.4f, %.4f], rho range [%.3f, %.3f]'
              % (thr, sel.sum(), CL[0, sel].min(), CH[0, sel].max(), CL[1, sel].min(), CH[1, sel].max(),
                 CL[2, sel].min(), CH[2, sel].max()))
np.savez_compressed(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'b4a_main_certified_boxes.npz'),
                    lo=CL, hi=CH, val=CV)
