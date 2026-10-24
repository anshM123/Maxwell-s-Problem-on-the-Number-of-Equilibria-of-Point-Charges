"""AUDIT_2 item 3: MIXED.md section 6.2, Lemma C (bordered Hessian), independent re-derivation and tests.

 (S) sympy, with ARBITRARY smooth functions beta(.), L1(.), L2(.): G(X) = beta(X1+X2) - cosh(L1(X1) - L2(X2)),
     tau = (beta' + S L2', S L1' - beta') with S = sinh(L1 - L2); check tau . grad G = 0 and
     -Hess G(tau, tau) - J = 0 identically, J = (U1+U2)^2 (P S^2 + C beta'^2) + S[U1'(beta'+S U2)^2 - U2'(S U1 - beta')^2]
     with C = cosh(L1 - L2) (= beta on the boundary), P = -beta''.
 (N) numerical, 40 digits, at boundary points of M' = {E < K} found by the auditor's own bisection on E itself (rays
     from t0 = (A1/2, A2/2)), for alpha = 1/2 (K = 9/2) and alpha = 10 (K = 484/21 = 23.048):
     - the curvature kappa of {E = K} computed DIRECTLY from derivatives of E (sublevel orientation) is compared with
       the vertex-3 formula: kappa = J/|grad G|^3 = I (U1 U2)^2 / |grad G|^3  (quantitative, not only the sign);
     - vertex 1: the same comparison in the coordinates X = (s, -t2): kappa_X computed directly from E(X1+X2, -X2) vs
       J1/|grad G1|^3 with J1 = U2^2 * (J1/U2^2) from the MIXED.md formula;
     - sign agreement of I, J1/U2^2 and kappa; number of negative values at alpha = 10;
     - mirror: I at the lower-arc point of shape (A1, A2) equals I at the upper-arc point of shape (A2, A1).
usage: python b3_lemmaC.py [n_shapes]"""
import sys

import mpmath as mp
import numpy as np
import sympy as sp

import b0_common as C

# ------------------------------------------------------------------ (S) symbolic
X1, X2 = sp.symbols('X1 X2', real=True)
bet = sp.Function('beta')
L1f, L2f = sp.Function('L1'), sp.Function('L2')
G = bet(X1 + X2) - sp.cosh(L1f(X1) - L2f(X2))
grad = [sp.diff(G, X1), sp.diff(G, X2)]
Hs = sp.Matrix([[sp.diff(G, X1, 2), sp.diff(G, X1, X2)], [sp.diff(G, X2, X1), sp.diff(G, X2, 2)]])
S_ = sp.sinh(L1f(X1) - L2f(X2))
C_ = sp.cosh(L1f(X1) - L2f(X2))
bp = sp.diff(bet(X1 + X2), X1)
bpp = sp.diff(bet(X1 + X2), X1, 2)
U1, U2 = sp.diff(L1f(X1), X1), sp.diff(L2f(X2), X2)
U1p, U2p = sp.diff(U1, X1), sp.diff(U2, X2)
tau = sp.Matrix([bp + S_ * U2, S_ * U1 - bp])
orth = sp.simplify(grad[0] * tau[0] + grad[1] * tau[1])
lhs = -(tau.T * Hs * tau)[0, 0]
J = (U1 + U2) ** 2 * (-bpp * S_ ** 2 + C_ * bp ** 2) + S_ * (U1p * (bp + S_ * U2) ** 2 - U2p * (S_ * U1 - bp) ** 2)
diff = sp.simplify(sp.expand(lhs - J))
print('(S) tau . grad G = %s;  -Hess G(tau,tau) - J = %s   (arbitrary beta, L1, L2)' % (orth, diff))
# vertex-3 normalisation: J/(U1 U2)^2 = I with c_i = 1/U_i, a_i = U_i'/U_i^2
c1, c2, a1, a2 = 1 / U1, 1 / U2, U1p / U1 ** 2, U2p / U2 ** 2
I_expr = (c1 + c2) ** 2 * (-bpp * S_ ** 2 + C_ * bp ** 2) + S_ * (a1 * (bp * c2 + S_) ** 2 - a2 * (S_ - bp * c1) ** 2)
print('(S) J/(U1 U2)^2 - I = %s' % sp.simplify(sp.expand(J / (U1 * U2) ** 2 - I_expr)))

# ------------------------------------------------------------------ (N) numerical
mp.mp.dps = 40
rng = np.random.default_rng(3)
nsh = int(sys.argv[1]) if len(sys.argv) > 1 else 40


def rand_shape(min_deg=1.0):
    while True:
        w = rng.dirichlet(np.ones(3)) * 180
        if w.min() > min_deg:
            A1, A2 = [mp.mpf(float(v)) * mp.pi / 180 for v in w[:2]]
            return (A1, A2, mp.pi - A1 - A2)


def kappa_X_vertex1(A, K, t1, t2):
    """curvature of {E = K} in X = (s, -t2) coordinates, sublevel orientation"""
    X1, X2 = t1 + t2, -t2
    f = lambda u, v: C.E_of_t(A, u + v, -v)
    E1 = mp.diff(f, (X1, X2), (1, 0)); E2 = mp.diff(f, (X1, X2), (0, 1))
    E11 = mp.diff(f, (X1, X2), (2, 0)); E22 = mp.diff(f, (X1, X2), (0, 2)); E12 = mp.diff(f, (X1, X2), (1, 1))
    return (E11 * E2 ** 2 - 2 * E12 * E1 * E2 + E22 * E1 ** 2) / (E1 ** 2 + E2 ** 2) ** mp.mpf(1.5)


for alpha in ('0.5', '10'):
    K = C.K_of_alpha(alpha)
    F = K / 2 - 1
    stats = dict(n=0, sign3=0, sign1=0, negI=0, negJ1=0, negk=0, relk3=mp.mpf(0), relk1=mp.mpf(0), Cbeta=mp.mpf(0),
                 mirror=mp.mpf(0), nonmono_rays=0)
    for k in range(nsh):
        A = rand_shape(1.0 if alpha == '0.5' else 2.0)
        A1, A2, A3 = A
        for j in range(10):
            ang = mp.mpf(float(rng.uniform(0, 2 * np.pi)))
            t, mono = C.boundary_point(A, K, ang)
            if t is None:
                continue
            stats['nonmono_rays'] += not mono
            t1, t2 = t
            s = t1 + t2
            kap = C.curvature_sign_direct(A, K, t1, t2)
            I, beta3, Cc = C.I_vertex3(A, K, t1, t2)
            stats['Cbeta'] = max(stats['Cbeta'], abs(Cc - beta3))
            U1v, U2v = C.U(t1, A1), C.U(t2, A2)
            psi = 2 * s - (A1 + A2)
            bp3 = -(F + 1) * mp.sin(psi)
            Sv = mp.sinh(C.Lam(t1, A1) - C.Lam(t2, A2))
            gradG = mp.sqrt((bp3 - Sv * U1v) ** 2 + (bp3 + Sv * U2v) ** 2)
            k3 = I * (U1v * U2v) ** 2 / gradG ** 3
            stats['relk3'] = max(stats['relk3'], abs(k3 - kap) / (abs(kap) + mp.mpf('1e-30')))
            J1n, beta1, C1 = C.J1n_vertex1(A, K, t1, t2)
            U3v = C.U(s, A1 + A2)
            bp1 = -(F + 1) * mp.sin(2 * t1 - A1)
            S1 = mp.sinh(C.Lam(s, A1 + A2) - C.Lam(t2, A2))
            gradG1 = mp.sqrt((bp1 - S1 * U3v) ** 2 + (bp1 + S1 * (-U2v)) ** 2)
            k1 = J1n * U2v ** 2 / gradG1 ** 3
            kX = kappa_X_vertex1(A, K, t1, t2)
            stats['relk1'] = max(stats['relk1'], abs(k1 - kX) / (abs(kX) + mp.mpf('1e-30')))
            stats['n'] += 1
            stats['sign3'] += mp.sign(I) != mp.sign(kap)
            stats['sign1'] += mp.sign(J1n) != mp.sign(kap)
            stats['negI'] += I < 0
            stats['negJ1'] += J1n < 0
            stats['negk'] += kap < 0
            # mirror: the same geometric point seen from the mirrored shape (A2, A1) is (t2, t1); it lies on the
            # other arc there; I must agree
            Am = (A2, A1, A3)
            Im, _, _ = C.I_vertex3(Am, K, t2, t1)
            stats['mirror'] = max(stats['mirror'], abs(Im - I) / (abs(I) + mp.mpf('1e-30')))
    print('alpha %s (K = %s): %d boundary points on %d shapes; rays where E was not monotone before the crossing: %d'
          % (alpha, mp.nstr(K, 6), stats['n'], nsh, stats['nonmono_rays']))
    print('   max |cosh Delta - beta(s)| at the points: %s' % mp.nstr(stats['Cbeta'], 3))
    print('   vertex 3: max rel |I (U1U2)^2/|grad G|^3 - kappa_direct| = %s; sign(I) != sign(kappa): %d'
          % (mp.nstr(stats['relk3'], 3), stats['sign3']))
    print('   vertex 1: max rel |J1/|grad G1|^3 - kappa_X direct| = %s; sign(J1) != sign(kappa): %d'
          % (mp.nstr(stats['relk1'], 3), stats['sign1']))
    print('   negative values: I %d, J1 %d, kappa %d;  mirror |I(shape, t) - I(mirror shape, swapped t)| rel max %s'
          % (stats['negI'], stats['negJ1'], stats['negk'], mp.nstr(stats['mirror'], 3)))
