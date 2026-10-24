"""AUDIT_2 shared helpers (the auditor's own implementations; nothing from code).

Mixed-sign setting of MIXED.md: q1, q2 > 0 > q3, triangle p1 p2 p3 counter-clockwise with angles A1, A2, A3.
Everything is computed from first principles in mpmath:
  - signed angles sigma_i = arg((p_k - x)/(p_j - x)) for cyclic (i, j, k), signed barycentrics, circumcircle power,
    D = det S from the vectors;
  - the inverse of the angle map t -> x on R+ is computed GEOMETRICALLY (inversion at p3 turns the two angle loci into
    lines), so Lemma P is not assumed;
  - E, U_i, Lambda_i, phi_i from their definitions;
  - boundary points of M' = {E < K} are found by bisection on E itself along rays from t0 = (A1/2, A2/2).
"""
import mpmath as mp

mp.mp.dps = 40


def shape(A1, A2):
    A1, A2 = mp.mpf(A1), mp.mpf(A2)
    A3 = mp.pi - A1 - A2
    L13 = mp.sin(A2) / mp.sin(A3)
    P = [mp.mpc(0, 0), mp.mpc(1, 0), mp.mpc(L13 * mp.cos(A1), L13 * mp.sin(A1))]
    return (A1, A2, A3), P


def cross(u, v):
    return (u.conjugate() * v).imag


def circum(P):
    a, b, c = P
    d = 2 * (a.real * (b.imag - c.imag) + b.real * (c.imag - a.imag) + c.real * (a.imag - b.imag))
    ux = (abs(a) ** 2 * (b.imag - c.imag) + abs(b) ** 2 * (c.imag - a.imag) + abs(c) ** 2 * (a.imag - b.imag)) / d
    uy = (abs(a) ** 2 * (c.real - b.real) + abs(b) ** 2 * (a.real - c.real) + abs(c) ** 2 * (b.real - a.real)) / d
    O = mp.mpc(ux, uy)
    return O, abs(a - O)


def geom(P, x):
    """signed angles, signed barycentrics, distances, P_x, D (= det S), 2*Delta"""
    D2 = cross(P[1] - P[0], P[2] - P[0])
    sig, lam, r = [], [], []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        sig.append(mp.arg((P[k] - x) / (P[j] - x)))
        lam.append(cross(P[j] - x, P[k] - x) / D2)
        r.append(abs(x - P[i]))
    O, R = circum(P)
    Px = R ** 2 - abs(x - O) ** 2
    S = mp.matrix(2, 2)
    for i in range(3):
        u = (x - P[i]) / r[i]
        S[0, 0] += lam[i] * u.real ** 2
        S[0, 1] += lam[i] * u.real * u.imag
        S[1, 1] += lam[i] * u.imag ** 2
    S[1, 0] = S[0, 1]
    return dict(sig=sig, lam=lam, r=r, Px=Px, D=mp.det(S), S=S, D2=D2, O=O, R=R)


def x_of_t(P, t1, t2):
    """the point x in R+ with sigma1 = t1, sigma2 = t2 (inversion at p3: both loci become lines)"""
    d1, d2 = P[0] - P[2], P[1] - P[2]
    a = d1 * mp.expj(-t1)
    bb = -d2 * mp.expj(t2)
    rhs = d1 - d2
    M = mp.matrix([[a.real, bb.real], [a.imag, bb.imag]])
    sol = mp.lu_solve(M, mp.matrix([rhs.real, rhs.imag]))
    r1, r2 = sol[0], sol[1]
    Y = (1 - r1 * mp.expj(-t1)) / d2
    return P[2] + 1 / Y, (r1, r2)


def U(t, a):
    return mp.cot(t) + mp.cot(a - t)


def Lam(t, a):
    return mp.log(mp.sin(t) / mp.sin(a - t))


def E_of_t(A, t1, t2):
    A1, A2, A3 = A
    s = t1 + t2
    U1, U2, U3 = U(t1, A1), U(t2, A2), U(s, A1 + A2)
    return U1 * U2 - U3 * (U1 + U2)


def b_of_alpha(alpha):
    return 2 * mp.mpf(alpha) + 1


def K_of_alpha(alpha):
    b = b_of_alpha(alpha)
    return (b + 1) ** 2 / b


def gradH(A, t1, t2, b):
    A1, A2, A3 = A
    s = t1 + t2
    phi1 = mp.log(mp.sin(t1)) + b * mp.log(mp.sin(A1 - t1))
    phi2 = mp.log(mp.sin(t2)) + b * mp.log(mp.sin(A2 - t2))
    phi3 = mp.log(mp.sin(s)) + b * mp.log(mp.sin(s + A3))
    return phi1 - phi3, phi2 - phi3


def boundary_point(A, K, ang, iters=150):
    """first crossing of E = K along the ray t0 + r (cos ang, sin ang), t0 = (A1/2, A2/2); returns t (or None if the
    ray leaves Q first) and a flag telling whether E was monotone along the sampled part of the ray before the crossing"""
    A1, A2, A3 = A
    t0 = (A1 / 2, A2 / 2)
    e = (mp.cos(ang), mp.sin(ang))

    def inQ(r):
        t1, t2 = t0[0] + r * e[0], t0[1] + r * e[1]
        return 0 < t1 < A1 and 0 < t2 < A2

    # largest r staying in Q
    rmax = mp.inf
    for k in range(2):
        if e[k] > 0:
            rmax = min(rmax, ((A1, A2)[k] - t0[k]) / e[k])
        elif e[k] < 0:
            rmax = min(rmax, -t0[k] / e[k])
    lo, hi = mp.mpf(0), rmax * (1 - mp.mpf('1e-30'))
    # coarse scan for the first crossing
    n = 400
    prev = E_of_t(A, *t0)
    mono = True
    found = None
    for k in range(1, n + 1):
        r = hi * mp.mpf(k) / n
        v = E_of_t(A, t0[0] + r * e[0], t0[1] + r * e[1])
        if v < prev:
            mono = False
        if v >= K:
            found = (hi * mp.mpf(k - 1) / n, r)
            break
        prev = v
    if found is None:
        return None, mono
    lo, hi = found
    for _ in range(iters):
        m = (lo + hi) / 2
        if E_of_t(A, t0[0] + m * e[0], t0[1] + m * e[1]) < K:
            lo = m
        else:
            hi = m
    r = (lo + hi) / 2
    return (t0[0] + r * e[0], t0[1] + r * e[1]), mono


def curvature_sign_direct(A, K, t1, t2):
    """curvature of the level curve {E = K} at (t1, t2), oriented so that > 0 means {E < K} is locally strictly convex:
    kappa = (E11 E2^2 - 2 E12 E1 E2 + E22 E1^2)/|grad E|^3 (for a sublevel set)"""
    f = lambda u, v: E_of_t(A, u, v)
    E1 = mp.diff(f, (t1, t2), (1, 0))
    E2 = mp.diff(f, (t1, t2), (0, 1))
    E11 = mp.diff(f, (t1, t2), (2, 0))
    E22 = mp.diff(f, (t1, t2), (0, 2))
    E12 = mp.diff(f, (t1, t2), (1, 1))
    return (E11 * E2 ** 2 - 2 * E12 * E1 * E2 + E22 * E1 ** 2) / (E1 ** 2 + E2 ** 2) ** mp.mpf(1.5)


def I_vertex3(A, K, t1, t2):
    """the MIXED.md Lemma C quantity I (vertex 3), evaluated from t by the auditor (F = K/2 - 1)"""
    A1, A2, A3 = A
    F = K / 2 - 1
    s = t1 + t2
    psi = 2 * s - (A1 + A2)
    beta = ((F - 1) * mp.cos(A3) + (F + 1) * mp.cos(psi)) / 2
    bp = -(F + 1) * mp.sin(psi)
    P = 2 * (F + 1) * mp.cos(psi)
    L1, L2 = Lam(t1, A1), Lam(t2, A2)
    Dl = L1 - L2
    S, C = mp.sinh(Dl), mp.cosh(Dl)
    h1, h2 = mp.cosh(L1) + mp.cos(A1), mp.cosh(L2) + mp.cos(A2)
    c1, c2 = mp.sin(A1) / (2 * h1), mp.sin(A2) / (2 * h2)
    a1, a2 = mp.sinh(L1) / h1, mp.sinh(L2) / h2
    I = (c1 + c2) ** 2 * (P * S ** 2 + C * bp ** 2) + S * (a1 * (bp * c2 + S) ** 2 - a2 * (S - bp * c1) ** 2)
    return I, beta, C


def J1n_vertex1(A, K, t1, t2):
    """the MIXED.md vertex-1 quantity J1/U2^2 evaluated from t by the auditor"""
    A1, A2, A3 = A
    F = K / 2 - 1
    s = t1 + t2
    psi1 = 2 * t1 - A1
    beta1 = ((F + 1) * mp.cos(psi1) - (F - 1) * mp.cos(A1)) / 2
    bp = -(F + 1) * mp.sin(psi1)
    P = 2 * (F + 1) * mp.cos(psi1)
    L2, L3 = Lam(t2, A2), Lam(s, A1 + A2)
    Dl = L3 - L2
    S, C = mp.sinh(Dl), mp.cosh(Dl)
    U3 = U(s, A1 + A2)
    U3p = -U3 * (mp.cot(s) + mp.cot(s + A3))
    h2 = mp.cosh(L2) + mp.cos(A2)
    c2 = mp.sin(A2) / (2 * h2)
    a2 = mp.sinh(L2) / h2
    J = (1 - U3 * c2) ** 2 * (P * S ** 2 + C * bp ** 2) + S * (U3p * (bp * c2 - S) ** 2 - a2 * (S * U3 - bp) ** 2)
    return J, beta1, C
