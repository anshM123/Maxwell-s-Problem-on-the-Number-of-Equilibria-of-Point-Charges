"""AUDIT_3 reference implementation (the auditor's own; nothing from code and nothing from audit/b0_common.py).

Everything is computed from the definition of E in angle coordinates,
    U_i(t) = cot t + cot(a_i - t)   (a_1 = A1, a_2 = A2, a_3 = A1 + A2, U_3 evaluated at s = t1 + t2),
    E = U1 U2 - U3 (U1 + U2),        M' = {E < K},  K = 2(F + 1),
with exact analytic first and second derivatives of E. The boundary point is found by root finding on E itself along
the relevant line (no closed-form quadratic, no Lambda-bisection):
  vertex 3 (main / band charts): the line s = t1 + t2 = (A12 + psi)/2, psi = rho psi_max; on this line E is unimodal
      (minimum where Lambda1 = Lambda2), and the upper arc (Delta > 0) is the root with the larger t1;
  vertex 1 (corner chart): the line t1 = (A1 + psi1)/2, psi1 = rho1 psi1_max; E is unimodal in t2, and sig = +1
      (S = sinh(Lambda3 - Lambda2) > 0) is the root with the smaller t2.
The curvature numerator of the sublevel set {E < K} is
    N = E11 E2^2 - 2 E12 E1 E2 + E22 E1^2      (> 0  <=>  {E < K} is locally strictly convex there).
With G3 = beta(s) - cosh(Lambda1 - Lambda2) = g3 (K - E), g3 = sin s sin(s + A3)/2 (product representation), the
bordered-Hessian quantity of Lemma C is J = g3^3 N, so the certificate's I = J/(U1 U2)^2 = g3^3 N/(U1 U2)^2.
In the vertex-1 coordinates X = (s, -t2), G1 = g1 (K - E) with g1 = sin t1 sin(A1 - t1)/2, and the same computation
(the linear change of coordinates maps (E_X2, -E_X1) to (-E2, E1)) gives J1 = g1^3 N, so J1/U2^2 = g1^3 N/U2^2.
These normalisations are derived independently of MIXED.md's Lemma C algebra; agreement with the certificate's
enclosures therefore tests both the closed forms and the algebraic reduction."""
import mpmath as mp


def cot(x):
    return mp.cos(x) / mp.sin(x)


def csc2(x):
    return 1 / mp.sin(x) ** 2


def U_parts(t, a):
    """U, U', U'' of U(t) = cot t + cot(a - t)"""
    u = cot(t) + cot(a - t)
    up = -csc2(t) + csc2(a - t)
    upp = 2 * csc2(t) * cot(t) + 2 * csc2(a - t) * cot(a - t)
    return u, up, upp


def E_val(A1, A2, t1, t2):
    U1 = cot(t1) + cot(A1 - t1)
    U2 = cot(t2) + cot(A2 - t2)
    s = t1 + t2
    U3 = cot(s) + cot(A1 + A2 - s)
    return U1 * U2 - U3 * (U1 + U2)


def E_derivs(A1, A2, t1, t2):
    U1, U1p, U1pp = U_parts(t1, A1)
    U2, U2p, U2pp = U_parts(t2, A2)
    U3, U3p, U3pp = U_parts(t1 + t2, A1 + A2)
    E = U1 * U2 - U3 * (U1 + U2)
    E1 = U1p * U2 - U3p * (U1 + U2) - U3 * U1p
    E2 = U1 * U2p - U3p * (U1 + U2) - U3 * U2p
    E11 = U1pp * U2 - U3pp * (U1 + U2) - 2 * U3p * U1p - U3 * U1pp
    E22 = U1 * U2pp - U3pp * (U1 + U2) - 2 * U3p * U2p - U3 * U2pp
    E12 = U1p * U2p - U3pp * (U1 + U2) - U3p * U2p - U3p * U1p
    return dict(E=E, E1=E1, E2=E2, E11=E11, E12=E12, E22=E22, U1=U1, U2=U2, U3=U3)


def curv_num(d):
    return d['E11'] * d['E2'] ** 2 - 2 * d['E12'] * d['E1'] * d['E2'] + d['E22'] * d['E1'] ** 2


def _golden_min(f, lo, hi, iters):
    g = (mp.sqrt(5) - 1) / 2
    a, b = lo, hi
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(iters):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = f(d)
    return (a + b) / 2


def _bisect_cross(f, inside, outside, K, iters):
    """f(inside) < K <= f(outside); returns the crossing point"""
    a, b = inside, outside
    for _ in range(iters):
        m = (a + b) / 2
        if f(m) < K:
            a = m
        else:
            b = m
    return (a + b) / 2


def iters_for_prec():
    return int(mp.mp.prec * 1.6) + 40


def upper_arc_point(A1, A2, s, K):
    """the point of {E = K} on the line t1 + t2 = s with Delta > 0 (larger t1)"""
    lo = max(mp.mpf(0), s - A2)
    hi = min(A1, s)
    w = hi - lo
    f = lambda t1: E_val(A1, A2, t1, s - t1)
    # keep away from the singular ends by a relative margin far below the working precision of the answer
    eps = w * mp.mpf(2) ** (-mp.mp.prec // 2 - 10)
    tmin = _golden_min(f, lo + eps, hi - eps, iters_for_prec())
    if f(tmin) >= K:
        return None
    t1 = _bisect_cross(f, tmin, hi - eps * mp.mpf(2) ** (-mp.mp.prec // 3), K, iters_for_prec())
    return t1, s - t1


def corner_point(A1, A2, t1, K, sig):
    """the point of {E = K} on the line t1 = const; sig = +1 -> smaller t2, sig = -1 -> larger t2"""
    f = lambda t2: E_val(A1, A2, t1, t2)
    eps = A2 * mp.mpf(2) ** (-mp.mp.prec // 2 - 10)
    tmin = _golden_min(f, eps, A2 - eps, iters_for_prec())
    if f(tmin) >= K:
        return None
    tiny = A2 * mp.mpf(2) ** (-mp.mp.prec + 20)
    out = tiny if sig > 0 else A2 - tiny
    t2 = _bisect_cross(f, tmin, out, K, iters_for_prec())
    return t1, t2


def psi_max3(A3, F):
    return mp.acos((2 - (F - 1) * mp.cos(A3)) / (F + 1))


def psi_max1(A1, F):
    return mp.acos((2 + (F - 1) * mp.cos(A1)) / (F + 1))


def I_ref_main(A12, lam, rho, F):
    """certificate quantity I (vertex 3, upper arc) at shape (A12, lam), boundary parameter rho, level F"""
    A12, lam, rho, F = mp.mpf(A12), mp.mpf(lam), mp.mpf(rho), mp.mpf(F)
    A1, A2 = lam * A12, (1 - lam) * A12
    A3 = mp.pi - A12
    K = 2 * (F + 1)
    psi = rho * psi_max3(A3, F)
    s = (A12 + psi) / 2
    pt = upper_arc_point(A1, A2, s, K)
    if pt is None:
        return None
    t1, t2 = pt
    d = E_derivs(A1, A2, t1, t2)
    g3 = mp.sin(s) * mp.sin(s + A3) / 2
    return g3 ** 3 * curv_num(d) / (d['U1'] * d['U2']) ** 2


def I4_ref_band(e, lam, rho, F):
    v = I_ref_main(e, lam, rho, F)
    return None if v is None else v / mp.mpf(e) ** 4


def J1n_ref_corner(r, th, rho1, F, sig):
    r, th, rho1, F = mp.mpf(r), mp.mpf(th), mp.mpf(rho1), mp.mpf(F)
    A1, A2, A3 = mp.pi - r, th * r, (1 - th) * r
    K = 2 * (F + 1)
    psi1 = rho1 * psi_max1(A1, F)
    t1 = (A1 + psi1) / 2
    pt = corner_point(A1, A2, t1, K, sig)
    if pt is None:
        return None
    t1, t2 = pt
    d = E_derivs(A1, A2, t1, t2)
    g1 = mp.sin(t1) * mp.sin(A1 - t1) / 2
    return g1 ** 3 * curv_num(d) / d['U2'] ** 2


def F_of_alpha_exact(alpha):
    """F = (b+1)^2/(2b) - 1 with b = 2 alpha + 1, as an exact Fraction for rational alpha"""
    from fractions import Fraction
    b = 2 * Fraction(alpha) + 1
    return (b + 1) ** 2 / (2 * b) - 1
