"""Certificate for Claim M at the Coulomb level (K = 9/2): the boundary of M' has positive curvature for every
triangle, so M' is convex (MIXED.md sections 6-7b).

Shape space {A1, A2, A3 >= 0, sum pi} is covered by three charts (e0, r0 below):
  band   : A1 + A2 <= e0                        -> mixed_charts.I4_band  (I/e^4, upper arc; lower arc = mirror)
  corner : A1 >= pi - r0                        -> mixed_charts.J1n_corner (vertex-1 form, both arcs)
           A2 >= pi - r0 is the mirror image of the corner chart (A1 <-> A2), so it needs no separate run
  main   : A1 + A2 >= e0, A1, A2 <= pi - r0     -> mixed_cert.I_iv (upper arc; lower arc = mirror)
On every box the evaluator returns an interval enclosure; a box is certified when its lower end is > 0. Failed boxes
are bisected along their widest normalised side, up to a round limit. The run must end with no failed box.
usage: python mixed_cert_all.py [e0] [r0] [max_rounds]
"""
import math
import os
import sys
import time

import numpy as np

import ia
import mixed_cert
import mixed_charts


def branch_and_bound(name, f, domain, grid, drop=None, max_rounds=60, cap=6_000_000):
    edges = [np.linspace(lo, hi, n + 1) for (lo, hi), n in zip(domain, grid)]
    lo = np.array(np.meshgrid(*[e[:-1] for e in edges], indexing='ij')).reshape(len(domain), -1)
    hi = np.array(np.meshgrid(*[e[1:] for e in edges], indexing='ij')).reshape(len(domain), -1)
    for k, (a, b) in enumerate(domain):                  # make the outer edges exact
        lo[k][lo[k] <= a] = a
        hi[k][hi[k] >= b] = b
    scale = np.array([b - a for a, b in domain])[:, None]
    n_cert, n_drop, minval, t0 = 0, 0, np.inf, time.time()
    for rnd in range(max_rounds):
        if drop is not None:
            d = drop(lo, hi)
            n_drop += d.sum()
            lo, hi = lo[:, ~d], hi[:, ~d]
        if lo.shape[1] == 0:
            break
        with np.errstate(all='ignore'):
            res = f(*[(lo[k], hi[k]) for k in range(lo.shape[0])])
        ok = res[0] > 0
        ok &= np.isfinite(res[0])
        if ok.any():
            minval = min(minval, res[0][ok].min())
        n_cert += ok.sum()
        lo, hi = lo[:, ~ok], hi[:, ~ok]
        if lo.shape[1] == 0 or 2 * lo.shape[1] > cap:
            break
        w = (hi - lo) / scale
        k = np.argmax(w, axis=0)
        m = 0.5 * (lo + hi)
        idx = np.arange(lo.shape[1])
        lo2, hi2 = lo.copy(), hi.copy()
        hi[k, idx] = m[k, idx]
        lo2[k, idx] = m[k, idx]
        lo, hi = np.hstack([lo, lo2]), np.hstack([hi, hi2])
    print('%-14s certified %9d boxes, dropped (other chart) %8d, FAILED %d, min lower bound %.3e  (%d rounds, %.0fs)'
          % (name, n_cert, n_drop, lo.shape[1], minval, rnd + 1, time.time() - t0), flush=True)
    return lo, hi, n_cert


if __name__ == '__main__':
    e0 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.2
    r0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
    rounds = int(sys.argv[3]) if len(sys.argv) > 3 else 60
    pi_lo = ia.PI[0]
    # main chart: drop boxes lying entirely in a corner region (A1 >= pi - r0 or A2 >= pi - r0 on the whole box)
    def drop_main(lo, hi):
        A1min = lo[1] * lo[0]
        A2min = (1.0 - hi[1]) * lo[0]
        return (A1min >= math.pi - r0 + 1e-15) | (A2min >= math.pi - r0 + 1e-15)
    failures = {}
    # upper edge: the float just above pi (math.pi < pi), so that triangles with tiny A3 > 0 are covered (AUDIT_2)
    res = branch_and_bound('main', mixed_cert.I_iv, [(e0, float(ia.up(np.pi))), (0.0, 1.0), (-1.0, 1.0)], (24, 16, 16),
                           drop=drop_main, max_rounds=rounds)
    failures['main'] = res[:2]
    res = branch_and_bound('band', mixed_charts.I4_band, [(0.0, e0), (0.0, 1.0), (-1.0, 1.0)], (8, 16, 16),
                           max_rounds=rounds)
    failures['band'] = res[:2]
    for sig in (+1, -1):
        res = branch_and_bound('corner %+d' % sig, lambda r, th, rh: mixed_charts.J1n_corner(r, th, rh, sig),
                               [(0.0, r0), (0.0, 1.0), (-1.0, 1.0)], (8, 8, 16), max_rounds=rounds)
        failures['corner%+d' % sig] = res[:2]
    total_fail = sum(v[0].shape[1] for v in failures.values())
    os.makedirs('../cert', exist_ok=True)
    np.savez_compressed('../cert/claimM_e0_%g_r0_%g.npz' % (e0, r0),
                        **{k + '_lo': v[0] for k, v in failures.items()}, **{k + '_hi': v[1] for k, v in failures.items()})
    print('TOTAL failed boxes: %d  ->  %s' % (total_fail, 'CLAIM M CERTIFIED (Coulomb)' if total_fail == 0 else 'not complete'))
    for k, (lo, hi) in failures.items():
        if lo.shape[1]:
            print('  %s failures span: %s' % (k, [(round(float(lo[i].min()), 6), round(float(hi[i].max()), 6))
                                                 for i in range(lo.shape[0])]))
