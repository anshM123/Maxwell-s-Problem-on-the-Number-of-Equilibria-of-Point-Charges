"""AUDIT_3, item 1/5: self-test of the reference c0_ref.py and the F-generality of the formulas used by the certificate.

For F in {13/12 (alpha = 1/4), 5/4 (Coulomb), 1.8, 13/5 (alpha = 2)} and random shapes (including near-degenerate ones):
 (a) product representation  E = 2 (cosh Delta + cos A3)/(sin s sin(s + A3))  at random t in Q (no boundary needed);
 (b) the boundary point found by root finding on E satisfies cosh(Lambda1 - Lambda2) = beta(s) with
     beta(s) = [(F-1) cos A3 - (F+1) cos(2s + A3)]/2   (region formula (eq:region) for general F) and Delta > 0;
 (c) I_ref = g3^3 N/(U1 U2)^2 (from E's derivatives) equals the closed formula (eq:I) of the paper, transcribed here
     by the auditor, at that point: this checks Lemma C's reduction for general F;
 (d) the same for the vertex-1 quantity J1/U2^2 (eq:J1), with the boundary point on t1 = (A1 + psi1)/2;
 (e) psi_max < A1 + A2 and psi1_max < A1 (Lemma disc (2) and its vertex-1 analogue);
 (f) the vertex-1 region formula cosh(Lambda2 - Lambda3) = beta1(t1) at the vertex-1 boundary point.
Output: audit/c1_ref_selftest.out"""
import random
import sys
from fractions import Fraction

import mpmath as mp

sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import c0_ref as R  # noqa: E402

mp.mp.dps = 60
rnd = random.Random(3)


def Lam(t, a):
    return mp.log(mp.sin(t) / mp.sin(a - t))


def I_formula(A1, A2, A3, F, t1, t2):
    s = t1 + t2
    psi = 2 * s + A3 - mp.pi
    bp = -(F + 1) * mp.sin(psi)
    P = 2 * (F + 1) * mp.cos(psi)
    L1, L2 = Lam(t1, A1), Lam(t2, A2)
    S, C = mp.sinh(L1 - L2), mp.cosh(L1 - L2)
    h1, h2 = mp.cosh(L1) + mp.cos(A1), mp.cosh(L2) + mp.cos(A2)
    c1, c2 = mp.sin(A1) / (2 * h1), mp.sin(A2) / (2 * h2)
    a1, a2 = mp.sinh(L1) / h1, mp.sinh(L2) / h2
    return (c1 + c2) ** 2 * (P * S ** 2 + C * bp ** 2) + S * (a1 * (bp * c2 + S) ** 2 - a2 * (S - bp * c1) ** 2)


def J1n_formula(A1, A2, A3, F, t1, t2):
    s = t1 + t2
    psi1 = 2 * t1 - A1
    bp = -(F + 1) * mp.sin(psi1)
    P = 2 * (F + 1) * mp.cos(psi1)
    L2, L3 = Lam(t2, A2), Lam(s, A1 + A2)
    S, C = mp.sinh(L3 - L2), mp.cosh(L3 - L2)
    U3 = R.cot(s) + R.cot(A1 + A2 - s)
    U3p = -U3 * (R.cot(s) + R.cot(s + A3))
    h2 = mp.cosh(L2) + mp.cos(A2)
    c2, a2 = mp.sin(A2) / (2 * h2), mp.sinh(L2) / h2
    return (1 - U3 * c2) ** 2 * (P * S ** 2 + C * bp ** 2) + S * (U3p * (bp * c2 - S) ** 2 - a2 * (S * U3 - bp) ** 2)


def rand_shape(kind):
    if kind == 'generic':
        A12 = mp.mpf(rnd.uniform(0.05, 3.1))
        lam = mp.mpf(rnd.uniform(0.02, 0.98))
    elif kind == 'thin':          # A1 + A2 small
        A12 = mp.mpf(10) ** rnd.uniform(-5, -1)
        lam = mp.mpf(rnd.uniform(0.0, 1.0))
    elif kind == 'flatA3':        # A3 small
        A12 = mp.pi - mp.mpf(10) ** rnd.uniform(-9, -2)
        lam = mp.mpf(rnd.uniform(0.05, 0.95))
    else:                          # 'lopsided': one of A1, A2 tiny relative to the other
        A12 = mp.mpf(rnd.uniform(0.05, 3.0))
        u = mp.mpf(10) ** rnd.uniform(-7, -1)
        lam = u if rnd.random() < 0.5 else 1 - u
    return lam * A12, (1 - lam) * A12, mp.pi - A12


def main():
    Fs = [Fraction(13, 12), Fraction(5, 4), Fraction(9, 5), Fraction(13, 5)]
    worst = dict(prod=0, region=0, I=0, J=0, region1=0)
    fails = dict(psimax=0, psimax1=0, sign=0, none=0)
    n_pts = 0
    for Fq in Fs:
        F = mp.mpf(Fq.numerator) / Fq.denominator
        K = 2 * (F + 1)
        for kind in ('generic', 'thin', 'flatA3', 'lopsided'):
            for _ in range(25):
                A1, A2, A3 = rand_shape(kind)
                # (a) product representation at a random t in Q
                t1, t2 = A1 * mp.mpf(rnd.uniform(0.01, 0.99)), A2 * mp.mpf(rnd.uniform(0.01, 0.99))
                s = t1 + t2
                Dl = Lam(t1, A1) - Lam(t2, A2)
                E1 = R.E_val(A1, A2, t1, t2)
                E2 = 2 * (mp.cosh(Dl) + mp.cos(A3)) / (mp.sin(s) * mp.sin(s + A3))
                worst['prod'] = max(worst['prod'], abs(E1 / E2 - 1))
                # (e) psi_max bounds
                pm = R.psi_max3(A3, F)
                if not pm < A1 + A2:
                    fails['psimax'] += 1
                # (b), (c): a boundary point on the upper arc
                rho = mp.mpf(rnd.uniform(-1, 1))
                s = (A1 + A2 + rho * pm) / 2
                pt = R.upper_arc_point(A1, A2, s, K)
                if pt is None:
                    fails['none'] += 1
                    continue
                t1, t2 = pt
                Dl = Lam(t1, A1) - Lam(t2, A2)
                beta = ((F - 1) * mp.cos(A3) - (F + 1) * mp.cos(2 * s + A3)) / 2
                worst['region'] = max(worst['region'], abs(mp.cosh(Dl) / beta - 1))
                if not Dl > 0:
                    fails['sign'] += 1
                Iref = R.I_ref_main(A1 + A2, A1 / (A1 + A2), rho, F)
                If = I_formula(A1, A2, A3, F, t1, t2)
                worst['I'] = max(worst['I'], abs(Iref - If) / (abs(If) + mp.mpf(10) ** -40))
                n_pts += 1
            # (d), (f): vertex-1 representation on corner-type shapes and on generic shapes
            for k in range(25):
                if k % 2:
                    r = mp.mpf(10) ** rnd.uniform(-5, -1.3)
                    A1 = mp.pi - r
                    th = mp.mpf(rnd.uniform(0.001, 0.999))
                    A2, A3 = th * r, (1 - th) * r
                else:
                    A1, A2, A3 = rand_shape(kind)
                pm1 = R.psi_max1(A1, F)
                if not pm1 < A1:
                    fails['psimax1'] += 1
                rho1 = mp.mpf(rnd.uniform(-1, 1))
                t1 = (A1 + rho1 * pm1) / 2
                for sig in (1, -1):
                    pt = R.corner_point(A1, A2, t1, K, sig)
                    if pt is None:
                        fails['none'] += 1
                        continue
                    _, t2 = pt
                    s = t1 + t2
                    D1 = Lam(s, A1 + A2) - Lam(t2, A2)
                    if not (D1 > 0) == (sig > 0):
                        fails['sign'] += 1
                    beta1 = ((F + 1) * mp.cos(2 * t1 - A1) - (F - 1) * mp.cos(A1)) / 2
                    worst['region1'] = max(worst['region1'], abs(mp.cosh(D1) / beta1 - 1))
                    d = R.E_derivs(A1, A2, t1, t2)
                    g1 = mp.sin(t1) * mp.sin(A1 - t1) / 2
                    Jref = g1 ** 3 * R.curv_num(d) / d['U2'] ** 2
                    Jf = J1n_formula(A1, A2, A3, F, t1, t2)
                    worst['J'] = max(worst['J'], abs(Jref - Jf) / (abs(Jf) + mp.mpf(10) ** -40))
                    n_pts += 1
    out = ['c1_ref_selftest (60 digits): F in {13/12, 5/4, 9/5, 13/5}; generic, thin (A1+A2 down to 1e-5), '
           'flat (A3 down to 1e-9), lopsided (A2/A1 down to 1e-7) shapes; corner shapes r down to 1e-5',
           'boundary points tested: %d' % n_pts,
           '(a) product representation, max relative error: %s' % mp.nstr(worst['prod'], 3),
           '(b) region formula cosh Delta = beta(s) at the vertex-3 boundary point, max rel. error: %s'
           % mp.nstr(worst['region'], 3),
           '(f) vertex-1 region formula cosh(L3 - L2) = beta1(t1), max rel. error: %s' % mp.nstr(worst['region1'], 3),
           '(c) I from E-derivatives vs closed formula (eq:I), max rel. difference: %s' % mp.nstr(worst['I'], 3),
           '(d) J1/U2^2 from E-derivatives vs closed formula (eq:J1), max rel. difference: %s' % mp.nstr(worst['J'], 3),
           '(e) failures of psi_max < A1+A2: %d; of psi1_max < A1: %d' % (fails['psimax'], fails['psimax1']),
           'arc-sign failures: %d; boundary point not found: %d' % (fails['sign'], fails['none'])]
    print('\n'.join(out))


if __name__ == '__main__':
    main()
