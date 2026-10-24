"""Interval evaluations of the boundary-convexity inequality in the two singular charts (MIXED.md section 7b).

band:   A12 = e in [0, e0], lambda = A1/e, rho in [-1, 1], upper arc.  Returns I / e^4 through scaled quantities
        (psi = e*psih, S = e*Sh, beta' = e*bph, c_i = e*ch_i, Delta = e*Dh), all regular at e = 0, and
        (a1 - a2)/e enclosed by the mean value theorem:
          a1 - a2 = (L1 - L2) * da/dL(xi, A1) + sinh L2 (cos A2 - cos A1) / ((cosh L2 + cos A1)(cosh L2 + cos A2)),
          a(L, A) = sinh L/(cosh L + cos A),  da/dL = (1 + cos A cosh L)/(cosh L + cos A)^2,  xi in hull(L1, L2).
        I/e^4 = (ch1+ch2)^2 (P Sh^2 + C bph^2) + Sh^3 (a1-a2)/e + 2 bph Sh^2 (a1 ch2 + a2 ch1)
                + e Sh bph^2 (a1 ch2^2 - a2 ch1^2).

corner: A1 = pi - r, A2 = theta r, A3 = (1 - theta) r, r in [0, r0], theta in [0, 1], rho1 in [-1, 1], arc sig = +-1.
        Vertex-1 representation: M' = {cosh(Lambda3(s) - Lambda2(t2)) < beta1(t1)},
        beta1 = [(F+1) cos psi1 - (F-1) cos A1]/2, psi1 = 2 t1 - A1.  Boundary inequality (normalised by U2^2):
          J1n = (1 - U3 c2)^2 (P S^2 + C beta1'^2) + S [U3' (beta1' c2 - S)^2 - a2 (S U3 - beta1')^2]
        with S = sinh(Lambda3 - Lambda2) = sig sqrt(beta1^2 - 1), C = beta1, beta1' = -(F+1) sin psi1,
        P = 2(F+1) cos psi1, zeta = e^{Lambda2} the positive root of
          e^D sin(A1 - t1) z^2 + [e^D sin(A12 - t1) - sin(A2 + t1)] z - sin t1 = 0   (e^D = beta1 + S),
        t2 = arctan(zeta sin A2 / (1 + zeta cos A2)), s = t1 + t2, U3 = sin A3/(sin s sin(s + A3)),
        U3' = -U3 (cot s + cot(s + A3)), c2 = sin A2 / (2(cosh L2 + cos A2)), a2 = sinh L2/(cosh L2 + cos A2).
"""
import numpy as np

import ia

FM, FP = 0.25, 2.25


def I4_band(e, lam, rho):
    sh = e[0].shape
    one, two = ia.const(1.0, sh), ia.const(2.0, sh)
    q = ia.div(one, ia.const(3.0, sh))
    he = ia.scal(0.5, e)
    y = ia.mul(q, ia.mul(he, ia.sinc(he)))                      # q sin(e/2)
    pmh = ia.mul(ia.mul(q, ia.sinc(he)), ia.asinc(y))            # psi_max / e
    psih = ia.mul(rho, pmh)                                      # psi / e
    xp = ia.scal(0.5, ia.mul(ia.add(one, rho), pmh))
    xm = ia.scal(0.5, ia.mul(ia.sub(one, rho), pmh))
    bm1h = ia.scal(FP, ia.mul(ia.mul(xp, ia.sinc(ia.mul(e, xp))), ia.mul(xm, ia.sinc(ia.mul(e, xm)))))
    bm1h = (np.maximum(bm1h[0], 0.0), bm1h[1])
    beta = ia.add(one, ia.mul(ia.sq(e), bm1h))
    Sh = ia.sqrt(ia.mul(bm1h, ia.add(beta, one)))
    bph = ia.scal(-FP, ia.mul(psih, ia.sinc(ia.mul(e, psih))))
    P = ia.scal(2 * FP, ia.cos(ia.mul(e, psih)))
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


def J1n_corner(r, th, rho1, sig):
    sh = r[0].shape
    one, two = ia.const(1.0, sh), ia.const(2.0, sh)
    A2 = ia.mul(th, r)
    A3 = ia.mul(ia.sub(one, th), r)
    A1 = ia.sub(ia.PI, r)
    cA1 = ia.neg(ia.cos(r))
    k1 = ia.div(ia.add(two, ia.scal(FM, cA1)), ia.const(FP, sh))
    pm = ia.arccos(k1)
    psi = ia.mul(rho1, pm)
    t1 = ia.scal(0.5, ia.add(A1, psi))
    sp = ia.sin(ia.scal(0.5, ia.mul(ia.add(one, rho1), pm)))
    sm = ia.sin(ia.scal(0.5, ia.mul(ia.sub(one, rho1), pm)))
    bm1 = ia.scal(FP, ia.mul(sp, sm))
    bm1 = (np.maximum(bm1[0], 0.0), bm1[1])
    beta = ia.add(one, bm1)
    Sabs = ia.sqrt(ia.mul(bm1, ia.add(bm1, two)))
    S = Sabs if sig > 0 else ia.neg(Sabs)
    bp = ia.scal(-FP, ia.sin(psi))
    P = ia.scal(2 * FP, ia.cos(psi))
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
