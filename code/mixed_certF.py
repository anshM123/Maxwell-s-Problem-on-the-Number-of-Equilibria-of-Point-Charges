"""Claim M for a RANGE of exponents: the same three charts as mixed_cert_all.py (main / band / corner), with
F = K/2 - 1 an interval parameter. Here K = (b+1)^2/b, b = 2 alpha + 1, so F = (b+1)^2/(2b) - 1 is increasing in
alpha. Boxes are 4-dimensional, (shape, shape, rho, F).

The formulas are those of mixed_cert.I_iv and mixed_charts.I4_band / J1n_corner, with (F - 1), (F + 1) and
q = sqrt((F - 1)/(F + 1)) as intervals.
usage: python mixed_certF.py alpha_lo alpha_hi [e0] [r0] [rounds]
"""
import math
import sys
import time

import numpy as np

import ia


def _F_of_alpha(a):
    b = 2 * a + 1
    return (b + 1) ** 2 / (2 * b) - 1


def I_main(A12, lam, rho, F):
    sh = A12[0].shape
    one, two = ia.const(1.0, sh), ia.const(2.0, sh)
    Fm, Fp = ia.sub(F, one), ia.add(F, one)
    A1 = ia.mul(lam, A12)
    A2 = ia.mul(ia.sub(one, lam), A12)
    g3 = ia.cos(ia.sub(ia.PI, A12))
    kap = ia.div(ia.sub(two, ia.mul(Fm, g3)), Fp)
    psim = ia.arccos(kap)
    psi = ia.mul(rho, psim)
    sp = ia.sin(ia.scal(0.5, ia.mul(ia.add(one, rho), psim)))
    sm = ia.sin(ia.scal(0.5, ia.mul(ia.sub(one, rho), psim)))
    bm1 = ia.mul(Fp, ia.mul(sp, sm))
    bm1 = (np.maximum(bm1[0], 0.0), bm1[1])
    beta = ia.add(one, bm1)
    S = ia.sqrt(ia.mul(bm1, ia.add(bm1, two)))
    emD = ia.recip(ia.add(beta, S))
    D = ia.arcsinh(S)
    bp = ia.neg(ia.mul(Fp, ia.sin(psi)))
    P = ia.scal(2.0, ia.mul(Fp, ia.cos(psi)))
    hpsi = ia.scal(0.5, psi)
    half12 = ia.scal(0.5, A12)
    dA = ia.scal(0.5, ia.mul(ia.sub(ia.scal(2.0, lam), one), A12))
    qa = ia.mul(emD, ia.sin(ia.sub(half12, hpsi)))
    qb = ia.add(ia.sin(ia.sub(dA, hpsi)), ia.mul(emD, ia.sin(ia.sub(ia.neg(dA), hpsi))))
    qc = ia.sin(ia.add(half12, hpsi))
    zeta = ia.div(ia.scal(2.0, qc), ia.add(qb, ia.sqrt(ia.add(ia.sq(qb), ia.scal(4.0, ia.mul(qa, qc))))))
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


def I4_band(e, lam, rho, F):
    sh = e[0].shape
    one = ia.const(1.0, sh)
    Fm, Fp = ia.sub(F, one), ia.add(F, one)
    q = ia.sqrt(ia.div(Fm, Fp))
    he = ia.scal(0.5, e)
    y = ia.mul(q, ia.mul(he, ia.sinc(he)))
    pmh = ia.mul(ia.mul(q, ia.sinc(he)), ia.asinc(y))
    psih = ia.mul(rho, pmh)
    xp = ia.scal(0.5, ia.mul(ia.add(one, rho), pmh))
    xm = ia.scal(0.5, ia.mul(ia.sub(one, rho), pmh))
    bm1h = ia.mul(Fp, ia.mul(ia.mul(xp, ia.sinc(ia.mul(e, xp))), ia.mul(xm, ia.sinc(ia.mul(e, xm)))))
    bm1h = (np.maximum(bm1h[0], 0.0), bm1h[1])
    beta = ia.add(one, ia.mul(ia.sq(e), bm1h))
    Sh = ia.sqrt(ia.mul(bm1h, ia.add(beta, one)))
    bph = ia.neg(ia.mul(Fp, ia.mul(psih, ia.sinc(ia.mul(e, psih)))))
    P = ia.scal(2.0, ia.mul(Fp, ia.cos(ia.mul(e, psih))))
    eS = ia.mul(e, Sh)
    Dh = ia.mul(Sh, ia.asinhc(eS))
    emD = ia.recip(ia.add(beta, eS))
    u1 = ia.scal(0.5, ia.sub(one, psih))
    qa = ia.mul(emD, ia.mul(u1, ia.sinc(ia.mul(e, u1))))
    l2 = ia.sub(ia.scal(2.0, lam), one)
    v1 = ia.scal(0.5, ia.sub(l2, psih))
    v2 = ia.scal(0.5, ia.sub(ia.neg(l2), psih))
    qb = ia.add(ia.mul(v1, ia.sinc(ia.mul(e, v1))), ia.mul(emD, ia.mul(v2, ia.sinc(ia.mul(e, v2)))))
    w1 = ia.scal(0.5, ia.add(one, psih))
    qc = ia.mul(w1, ia.sinc(ia.mul(e, w1)))
    zeta = ia.div(ia.scal(2.0, qc), ia.add(qb, ia.sqrt(ia.add(ia.sq(qb), ia.scal(4.0, ia.mul(qa, qc))))))
    L1 = ia.log(zeta)
    L2 = ia.sub(L1, ia.mul(e, Dh))
    lam2 = ia.sub(one, lam)
    A1, A2 = ia.mul(lam, e), ia.mul(lam2, e)
    cA1, cA2 = ia.cos(A1), ia.cos(A2)
    ch1, ch2 = ia.cosh(L1), ia.cosh(L2)
    h1, h2 = ia.add(ch1, cA1), ia.add(ch2, cA2)
    c1h = ia.div(ia.mul(lam, ia.sinc(A1)), ia.scal(2.0, h1))
    c2h = ia.div(ia.mul(lam2, ia.sinc(A2)), ia.scal(2.0, h2))
    a1 = ia.div(ia.sinh(L1), h1)
    a2 = ia.div(ia.sinh(L2), h2)
    chx = ia.cosh(ia.hull(L1, L2))
    dadL = ia.div(ia.add(one, ia.mul(cA1, chx)), ia.sq(ia.add(chx, cA1)))
    dcos = ia.scal(0.5, ia.mul(ia.mul(ia.sinc(he), l2), ia.sinc(ia.scal(0.5, ia.mul(l2, e)))))
    part2 = ia.mul(e, ia.div(ia.mul(ia.sinh(L2), dcos), ia.mul(ia.add(ch2, cA1), ia.add(ch2, cA2))))
    da = ia.add(ia.mul(Dh, dadL), part2)
    S2 = ia.sq(Sh)
    term1 = ia.mul(ia.sq(ia.add(c1h, c2h)), ia.add(ia.mul(P, S2), ia.mul(beta, ia.sq(bph))))
    term2 = ia.mul(ia.mul(S2, Sh), da)
    term3 = ia.scal(2.0, ia.mul(ia.mul(bph, S2), ia.add(ia.mul(a1, c2h), ia.mul(a2, c1h))))
    term4 = ia.mul(ia.mul(e, ia.mul(Sh, ia.sq(bph))), ia.sub(ia.mul(a1, ia.sq(c2h)), ia.mul(a2, ia.sq(c1h))))
    return ia.add(ia.add(term1, term2), ia.add(term3, term4))


def J1n_corner(r, th, rho1, F, sig):
    sh = r[0].shape
    one, two = ia.const(1.0, sh), ia.const(2.0, sh)
    Fm, Fp = ia.sub(F, one), ia.add(F, one)
    A2 = ia.mul(th, r)
    A3 = ia.mul(ia.sub(one, th), r)
    A1 = ia.sub(ia.PI, r)
    cA1 = ia.neg(ia.cos(r))
    k1 = ia.div(ia.add(two, ia.mul(Fm, cA1)), Fp)
    pm = ia.arccos(k1)
    psi = ia.mul(rho1, pm)
    t1 = ia.scal(0.5, ia.add(A1, psi))
    sp = ia.sin(ia.scal(0.5, ia.mul(ia.add(one, rho1), pm)))
    sm = ia.sin(ia.scal(0.5, ia.mul(ia.sub(one, rho1), pm)))
    bm1 = ia.mul(Fp, ia.mul(sp, sm))
    bm1 = (np.maximum(bm1[0], 0.0), bm1[1])
    beta = ia.add(one, bm1)
    Sabs = ia.sqrt(ia.mul(bm1, ia.add(bm1, two)))
    S = Sabs if sig > 0 else ia.neg(Sabs)
    bp = ia.neg(ia.mul(Fp, ia.sin(psi)))
    P = ia.scal(2.0, ia.mul(Fp, ia.cos(psi)))
    eD = ia.add(beta, S)
    A1mt1 = ia.scal(0.5, ia.sub(A1, psi))
    qa = ia.mul(eD, ia.sin(A1mt1))
    qb = ia.sub(ia.mul(eD, ia.sin(ia.add(A2, A1mt1))), ia.sin(ia.add(A2, t1)))
    qc = ia.sin(t1)
    zeta = ia.div(ia.scal(2.0, qc), ia.add(qb, ia.sqrt(ia.add(ia.sq(qb), ia.scal(4.0, ia.mul(qa, qc))))))
    L2 = ia.log(zeta)
    sA2, cA2 = ia.sin(A2), ia.cos(A2)
    t2 = ia.arctan(ia.div(ia.mul(zeta, sA2), ia.add(one, ia.mul(zeta, cA2))))
    s = ia.add(t1, t2)
    chL2 = ia.cosh(L2)
    h2 = ia.add(chL2, cA2)
    c2 = ia.div(sA2, ia.scal(2.0, h2))
    a2 = ia.div(ia.sinh(L2), h2)
    sA3 = ia.sin(A3)
    ss, sa = ia.sin(s), ia.sin(ia.add(s, A3))
    U3 = ia.div(sA3, ia.mul(ss, sa))
    cot_sum = ia.add(ia.div(ia.cos(s), ss), ia.div(ia.cos(ia.add(s, A3)), sa))
    U3p = ia.neg(ia.mul(U3, cot_sum))
    f1 = ia.sq(ia.sub(one, ia.mul(U3, c2)))
    main = ia.mul(f1, ia.add(ia.mul(P, ia.sq(S)), ia.mul(beta, ia.sq(bp))))
    tA = ia.mul(U3p, ia.sq(ia.sub(ia.mul(bp, c2), S)))
    tB = ia.mul(a2, ia.sq(ia.sub(ia.mul(S, U3), bp)))
    return ia.add(main, ia.mul(S, ia.sub(tA, tB)))


def bnb(name, f, domain, grid, drop=None, max_rounds=60, cap=8_000_000):
    edges = [np.linspace(lo, hi, n + 1) for (lo, hi), n in zip(domain, grid)]
    lo = np.array(np.meshgrid(*[e[:-1] for e in edges], indexing='ij')).reshape(len(domain), -1)
    hi = np.array(np.meshgrid(*[e[1:] for e in edges], indexing='ij')).reshape(len(domain), -1)
    for k, (a, b) in enumerate(domain):
        lo[k][lo[k] <= a] = a
        hi[k][hi[k] >= b] = b
    scale = np.array([max(b - a, 0.05) for a, b in domain])[:, None]   # floor: never split a near-degenerate F
    n_cert, t0 = 0, time.time()
    for rnd in range(max_rounds):
        if drop is not None:
            d = drop(lo, hi)
            lo, hi = lo[:, ~d], hi[:, ~d]
        if lo.shape[1] == 0:
            break
        with np.errstate(all='ignore'):
            res = f(*[(lo[k], hi[k]) for k in range(lo.shape[0])])
        ok = (res[0] > 0) & np.isfinite(res[0])
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
    print('%-12s certified %9d, FAILED %d  (%d rounds, %.0fs)' % (name, n_cert, lo.shape[1], rnd + 1, time.time() - t0),
          flush=True)
    return lo, hi


if __name__ == '__main__':
    a_lo, a_hi = float(sys.argv[1]), float(sys.argv[2])
    e0 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.2
    r0 = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05
    rounds = int(sys.argv[5]) if len(sys.argv) > 5 else 60
    F_lo, F_hi = _F_of_alpha(a_lo), _F_of_alpha(a_hi)
    F_dom = (ia.dn(F_lo), ia.up(F_hi))      # binary64 F(alpha) widened by 2 ulps; sweep_coverage.py checks the union
    print('alpha in [%g, %g]  ->  F in [%.6f, %.6f]' % (a_lo, a_hi, F_dom[0], F_dom[1]))

    def drop_main(lo, hi):
        return ((lo[1] * lo[0] >= math.pi - r0 + 1e-15) | ((1.0 - hi[1]) * lo[0] >= math.pi - r0 + 1e-15))
    fails = 0
    lo, hi = bnb('main', I_main, [(e0, math.pi), (0.0, 1.0), (-1.0, 1.0), F_dom], (24, 16, 16, 2), drop=drop_main,
                 max_rounds=rounds)
    fails += lo.shape[1]
    lo, hi = bnb('band', I4_band, [(0.0, e0), (0.0, 1.0), (-1.0, 1.0), F_dom], (8, 16, 16, 2), max_rounds=rounds)
    fails += lo.shape[1]
    for sig in (+1, -1):
        lo, hi = bnb('corner %+d' % sig, lambda r, th, rh, F: J1n_corner(r, th, rh, F, sig),
                     [(0.0, r0), (0.0, 1.0), (-1.0, 1.0), F_dom], (8, 8, 16, 2), max_rounds=rounds)
        fails += lo.shape[1]
    print('TOTAL failed: %d  ->  %s' % (fails, 'CLAIM M CERTIFIED for alpha in [%g, %g]' % (a_lo, a_hi) if fails == 0
                                         else 'not complete'))
