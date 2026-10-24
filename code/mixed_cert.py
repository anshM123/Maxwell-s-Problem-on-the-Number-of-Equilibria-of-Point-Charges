"""Certified check of the boundary inequality I > 0 (MIXED.md section 7b) on the main chart, by interval arithmetic.

Shape: A12 = A1 + A2 in [e0, pi], lambda = A1/A12 in [0, 1] (A3 = pi - A12).
Boundary point: rho in [-1, 1], psi = rho * psi_max, s = A12/2 + psi/2, on the upper arc Delta = +B(s).
The lower arc is the upper arc of the mirrored shape (A1 <-> A2, lambda -> 1 - lambda); I is invariant under this.
Coulomb: K = 9/2, F = K/2 - 1 = 5/4.

    beta - 1 = (F+1) sin(psi_max (1+rho)/2) sin(psi_max (1-rho)/2),  S = sqrt((beta-1)(beta+1)),  C = beta,
    beta' = -(F+1) sin psi,  P = 2 (F+1) cos psi,  Delta = arcsinh S,  e^{-Delta} = 1/(beta + S),
    zeta = 2 c0 / (b + sqrt(b^2 + 4 a c0)),  a = e^{-Delta} sin((A12 - psi)/2),
    b = sin((A1-A2-psi)/2) + e^{-Delta} sin((A2-A1-psi)/2),  c0 = sin((A12 + psi)/2),
    Lambda1 = log zeta, Lambda2 = Lambda1 - Delta,  h_i = cosh Lambda_i + cos A_i,
    c_i = sin A_i / (2 h_i),  a_i = sinh Lambda_i / h_i,
    I = (c1+c2)^2 (P S^2 + C beta'^2) + S [a1 (beta' c2 + S)^2 - a2 (S - beta' c1)^2].

usage: python mixed_cert.py e0 [max_rounds] [cap]
Writes the certified/failed box lists to ../cert/main_<e0>.npz and prints a summary.
"""
import math
import os
import sys
import time

import numpy as np

import ia

FM, FP = 0.25, 2.25          # F - 1, F + 1 for F = 5/4 (exact binary fractions)


def I_iv(A12, lam, rho):
    one = ia.const(1.0, A12[0].shape)
    two = ia.const(2.0, A12[0].shape)
    A1 = ia.mul(lam, A12)
    A2 = ia.mul(ia.sub(one, lam), A12)
    A3 = ia.sub(ia.PI, A12)
    g3 = ia.cos(A3)
    kap = ia.div(ia.sub(two, ia.scal(FM, g3)), ia.const(FP, A12[0].shape))
    psim = ia.arccos(kap)
    psi = ia.mul(rho, psim)
    sp = ia.sin(ia.scal(0.5, ia.mul(ia.add(one, rho), psim)))
    sm = ia.sin(ia.scal(0.5, ia.mul(ia.sub(one, rho), psim)))
    bm1 = ia.scal(FP, ia.mul(sp, sm))
    bm1 = (np.maximum(bm1[0], 0.0), bm1[1])
    beta = ia.add(one, bm1)
    S = ia.sqrt(ia.mul(bm1, ia.add(bm1, two)))
    emD = ia.recip(ia.add(beta, S))
    D = ia.arcsinh(S)
    bp = ia.scal(-FP, ia.sin(psi))
    P = ia.scal(2 * FP, ia.cos(psi))
    hpsi = ia.scal(0.5, psi)
    half12 = ia.scal(0.5, A12)
    dA = ia.scal(0.5, ia.mul(ia.sub(ia.scal(2.0, lam), one), A12))     # (A1 - A2)/2
    qa = ia.mul(emD, ia.sin(ia.sub(half12, hpsi)))
    qb = ia.add(ia.sin(ia.sub(dA, hpsi)), ia.mul(emD, ia.sin(ia.sub(ia.neg(dA), hpsi))))
    qc = ia.sin(ia.add(half12, hpsi))
    disc = ia.add(ia.sq(qb), ia.scal(4.0, ia.mul(qa, qc)))
    den = ia.add(qb, ia.sqrt(disc))
    zeta = ia.div(ia.scal(2.0, qc), den)
    L1 = ia.log(zeta)
    L2 = ia.sub(L1, D)
    h1 = ia.add(ia.cosh(L1), ia.cos(A1))
    h2 = ia.add(ia.cosh(L2), ia.cos(A2))
    c1 = ia.div(ia.sin(A1), ia.scal(2.0, h1))
    c2 = ia.div(ia.sin(A2), ia.scal(2.0, h2))
    a1 = ia.div(ia.sinh(L1), h1)
    a2 = ia.div(ia.sinh(L2), h2)
    main = ia.mul(ia.sq(ia.add(c1, c2)), ia.add(ia.mul(P, ia.sq(S)), ia.mul(beta, ia.sq(bp))))
    t1 = ia.mul(a1, ia.sq(ia.add(ia.mul(bp, c2), S)))
    t2 = ia.mul(a2, ia.sq(ia.sub(S, ia.mul(bp, c1))))
    return ia.add(main, ia.mul(S, ia.sub(t1, t2)))


def run(e0, max_rounds=40, cap=4_000_000, grid=(24, 16, 16), domain=None):
    if domain is None:
        domain = ((e0, math.pi), (0.0, 1.0), (-1.0, 1.0))
    edges = [np.linspace(lo, hi, n + 1) for (lo, hi), n in zip(domain, grid)]
    lo = np.array(np.meshgrid(edges[0][:-1], edges[1][:-1], edges[2][:-1], indexing='ij')).reshape(3, -1)
    hi = np.array(np.meshgrid(edges[0][1:], edges[1][1:], edges[2][1:], indexing='ij')).reshape(3, -1)
    hi[0][hi[0] >= math.pi] = ia.up(np.pi)        # cover the exact pi end
    certified, failed = [], []
    scale = np.array([domain[0][1] - domain[0][0], 1.0, 2.0])[:, None]
    t0 = time.time()
    for rnd in range(max_rounds):
        if lo.shape[1] == 0:
            break
        res = I_iv((lo[0], hi[0]), (lo[1], hi[1]), (lo[2], hi[2]))
        ok = res[0] > 0
        certified.append(np.vstack([lo[:, ok], hi[:, ok]]))
        lo, hi = lo[:, ~ok], hi[:, ~ok]
        print('round %2d: %8d boxes tested, %8d certified, %8d to split  (%.0fs)'
              % (rnd, ok.size, ok.sum(), (~ok).sum(), time.time() - t0), flush=True)
        if lo.shape[1] == 0:
            break
        if 2 * lo.shape[1] > cap:
            break
        w = (hi - lo) / scale
        d = np.argmax(w, axis=0)
        m = 0.5 * (lo + hi)
        idx = np.arange(lo.shape[1])
        lo2, hi2 = lo.copy(), hi.copy()
        hi[d, idx] = m[d, idx]
        lo2[d, idx] = m[d, idx]
        lo, hi = np.hstack([lo, lo2]), np.hstack([hi, hi2])
    failed = np.vstack([lo, hi])
    cert = np.hstack(certified) if certified else np.zeros((6, 0))
    return cert, failed


if __name__ == '__main__':
    e0 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.2
    rounds = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    cap = int(sys.argv[3]) if len(sys.argv) > 3 else 4_000_000
    cert, failed = run(e0, rounds, cap)
    os.makedirs('../cert', exist_ok=True)
    np.savez_compressed('../cert/main_%g.npz' % e0, cert=cert, failed=failed)
    print('certified boxes: %d; uncertified: %d' % (cert.shape[1], failed.shape[1]))
    if failed.shape[1]:
        lo, hi = failed[:3], failed[3:]
        print('uncertified region: A12 in [%.6f, %.6f], lambda in [%.6f, %.6f], rho in [%.4f, %.4f]'
              % (lo[0].min(), hi[0].max(), lo[1].min(), hi[1].max(), lo[2].min(), hi[2].max()))
        # where are they? histogram by (A12 near pi, lambda near 0/1, A12 near e0)
        nearpi = lo[0] > math.pi - 0.3
        print('  of which A12 > pi - 0.3: %d; lambda < 0.05 or > 0.95: %d; A12 < e0 + 0.05: %d'
              % (nearpi.sum(), ((lo[1] < 0.05) | (hi[1] > 0.95)).sum(), (lo[0] < e0 + 0.05).sum()))
