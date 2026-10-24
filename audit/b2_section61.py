"""AUDIT_2 item 2: MIXED.md section 6.1 (Coulomb unless stated), independent checks with 40-digit mpmath.

 (1) Lambda_i' = U_i, U_i = 2(cosh Lambda_i + cos a_i)/sin a_i, U_i' = 2 sinh Lambda_i U_i / sin a_i,
     t = arg(1 + e^{Lambda + i a});
 (2) product representation E = (2 U3(s)/sin A3)(cosh Delta + cos A3), 1/U3 = sin s sin(s+A3)/sin A3, and
     K - E = (2U3/sin A3)(beta(s) - cosh Delta) with beta = [(F-1) cos A3 - (F+1) cos(2s + A3)]/2;
 (3) vertex-1 and vertex-2 representations E = (2U1/sin A1)(cosh(Lambda2 - Lambda3) - cos A1), and the 1<->2 analogue;
 (4) (s, Delta) map: Jacobian det = -(U1 + U2); Delta strictly increasing along t1 + t2 = s, tending to -inf/+inf;
 (5) beta in terms of psi, cos psi_max = (2 - (F-1) cos A3)/(F+1), psi_max < A1 + A2, beta >= 1 <=> |psi| <= psi_max,
     beta - 1 = (F+1) sin(psi_max(1+rho)/2) sin(psi_max(1-rho)/2);
 (6) on M': min P = 2(F+1) cos psi, max |Lambda_i|, max |Delta| (claims: P >= 3.5, |Lambda_i|, |Delta| <= log 2);
 (7) vertex-3 quadratic for zeta = e^{Lambda1} at (s, Delta = +-B(s)): the root gives t with t1 + t2 = s, E(t) = K,
     Lambda1 - Lambda2 = Delta; compared with a boundary point found independently by bisection in t1 on the line
     t1 + t2 = s;
 (8) vertex-1 quadratic for zeta = e^{Lambda2} at (t1, D = +-arccosh beta1(t1)): E(t1, t2) = K and
     Lambda3(s) - Lambda2(t2) = D.
usage: python b2_section61.py"""
import mpmath as mp
import numpy as np

import b0_common as C

rng = np.random.default_rng(61)
alpha = mp.mpf('0.5')
b = C.b_of_alpha(alpha)
K = C.K_of_alpha(alpha)
F = K / 2 - 1
worst = {}


def upd(k, v):
    worst[k] = max(worst.get(k, mp.mpf(0)), abs(v))


def cnt(k, flag):
    worst[k] = worst.get(k, 0) + (1 if flag else 0)


def rand_shape(min_deg=0.5):
    while True:
        w = rng.dirichlet(np.ones(3)) * 180
        if w.min() > min_deg:
            A1, A2 = [mp.mpf(float(v)) * mp.pi / 180 for v in w[:2]]
            return (A1, A2, mp.pi - A1 - A2)


def beta_s(A, s):
    return ((F - 1) * mp.cos(A[2]) - (F + 1) * mp.cos(2 * s + A[2])) / 2


# (1)
for k in range(300):
    a = mp.mpf(float(rng.uniform(0.01, 3.1)))
    t = a * mp.mpf(float(rng.uniform(0.01, 0.99)))
    L = C.Lam(t, a)
    Uv = C.U(t, a)
    upd('(1) |Lambda\' - U|/U', (mp.diff(lambda u: C.Lam(u, a), t) - Uv) / Uv)
    upd('(1) |U - 2(cosh L + cos a)/sin a|/U', (Uv - 2 * (mp.cosh(L) + mp.cos(a)) / mp.sin(a)) / Uv)
    upd('(1) |U\' - 2 sinh L U/sin a|/|U\'|+1', (mp.diff(lambda u: C.U(u, a), t) - 2 * mp.sinh(L) * Uv / mp.sin(a)) / (abs(2 * mp.sinh(L) * Uv / mp.sin(a)) + 1))
    upd('(1) |t - arg(1 + e^(L + i a))|', t - mp.arg(1 + mp.exp(L + 1j * a)))

# (2), (3), (4)
for k in range(2000):
    A = rand_shape()
    A1, A2, A3 = A
    t1 = A1 * mp.mpf(float(rng.uniform(0.001, 0.999)))
    t2 = A2 * mp.mpf(float(rng.uniform(0.001, 0.999)))
    s = t1 + t2
    U1, U2, U3 = C.U(t1, A1), C.U(t2, A2), C.U(s, A1 + A2)
    E = U1 * U2 - U3 * (U1 + U2)
    L1, L2, L3 = C.Lam(t1, A1), C.Lam(t2, A2), C.Lam(s, A1 + A2)
    Dl = L1 - L2
    upd('(2) rel |E - (2U3/sin A3)(cosh Delta + cos A3)|', (E - 2 * U3 / mp.sin(A3) * (mp.cosh(Dl) + mp.cos(A3))) / E)
    upd('(2) rel |1/U3 - sin s sin(s+A3)/sin A3|', (1 / U3 - mp.sin(s) * mp.sin(s + A3) / mp.sin(A3)) * U3)
    upd('(2) rel |(K - E) - (2U3/sin A3)(beta(s) - cosh Delta)|',
        ((K - E) - 2 * U3 / mp.sin(A3) * (beta_s(A, s) - mp.cosh(Dl))) / (abs(K) + abs(E)))
    upd('(3) rel |E - (2U1/sin A1)(cosh(L2 - L3) - cos A1)|', (E - 2 * U1 / mp.sin(A1) * (mp.cosh(L2 - L3) - mp.cos(A1))) / E)
    upd('(3) rel |E - (2U2/sin A2)(cosh(L1 - L3) - cos A2)|', (E - 2 * U2 / mp.sin(A2) * (mp.cosh(L1 - L3) - mp.cos(A2))) / E)
    if k % 20 == 0:
        Jm = mp.matrix([[1, 1], [mp.diff(lambda u: C.Lam(u, A1), t1), -mp.diff(lambda u: C.Lam(u, A2), t2)]])
        upd('(4) |det d(s,Delta)/dt + (U1 + U2)| / (U1+U2)', (mp.det(Jm) + (U1 + U2)) / (U1 + U2))
# monotonicity of Delta along t1 + t2 = s, and limits
nonmono = 0
for k in range(200):
    A = rand_shape()
    A1, A2, A3 = A
    s = (A1 + A2) * mp.mpf(float(rng.uniform(0.01, 0.99)))
    lo, hi = max(mp.mpf(0), s - A2), min(A1, s)
    vals = []
    for j in range(1, 400):
        t1 = lo + (hi - lo) * mp.mpf(j) / 400
        vals.append(C.Lam(t1, A1) - C.Lam(s - t1, A2))
    nonmono += any(vals[j + 1] <= vals[j] for j in range(len(vals) - 1))
    ends = (C.Lam(lo + (hi - lo) * mp.mpf('1e-30'), A1) - C.Lam(s - lo - (hi - lo) * mp.mpf('1e-30'), A2),
            C.Lam(hi - (hi - lo) * mp.mpf('1e-30'), A1) - C.Lam(s - hi + (hi - lo) * mp.mpf('1e-30'), A2))
    cnt('(4) Delta at the two ends of the s-line not (< -60, > 60)', not (ends[0] < -60 and ends[1] > 60))
cnt('(4) Delta not strictly increasing along t1 + t2 = s (200 lines x 400 points)', nonmono)

# (5)
for k in range(3000):
    A = rand_shape(0.01)
    A1, A2, A3 = A
    A12 = A1 + A2
    cpm = (2 - (F - 1) * mp.cos(A3)) / (F + 1)
    pm = mp.acos(cpm)
    cnt('(5) psi_max >= A1 + A2', pm >= A12)
    s = A12 * mp.mpf(float(rng.uniform(0, 1)))
    psi = 2 * s - A12
    bet = beta_s(A, s)
    upd('(5) |beta - [(F-1)cos A3 + (F+1)cos psi]/2|', bet - ((F - 1) * mp.cos(A3) + (F + 1) * mp.cos(psi)) / 2)
    cnt('(5) (beta >= 1) != (|psi| <= psi_max)', (bet >= 1) != (abs(psi) <= pm))
    rho = mp.mpf(float(rng.uniform(-1, 1)))
    bet2 = ((F - 1) * mp.cos(A3) + (F + 1) * mp.cos(rho * pm)) / 2
    upd('(5) |beta - 1 - (F+1) sin(pm(1+rho)/2) sin(pm(1-rho)/2)|',
        bet2 - 1 - (F + 1) * mp.sin(pm * (1 + rho) / 2) * mp.sin(pm * (1 - rho) / 2))


# (6) sample M' through (s, Delta): t1 from Delta by bisection on the s-line
def t_of_sDelta(A, s, Dl):
    A1, A2, A3 = A
    lo, hi = max(mp.mpf(0), s - A2), min(A1, s)
    for _ in range(200):
        m = (lo + hi) / 2
        if C.Lam(m, A1) - C.Lam(s - m, A2) < Dl:
            lo = m
        else:
            hi = m
    t1 = (lo + hi) / 2
    return t1, s - t1


minP, maxL, maxD = mp.inf, mp.mpf(0), mp.mpf(0)
for k in range(300):
    A = rand_shape(0.2)
    A1, A2, A3 = A
    pm = mp.acos((2 - (F - 1) * mp.cos(A3)) / (F + 1))
    for j in range(12):
        rho = mp.mpf(float(rng.uniform(-1, 1)))
        s = (A1 + A2 + rho * pm) / 2
        B = mp.acosh(max(beta_s(A, s), mp.mpf(1)))
        Dl = B * mp.mpf(float(rng.uniform(-1, 1)))
        t1, t2 = t_of_sDelta(A, s, Dl)
        minP = min(minP, 2 * (F + 1) * mp.cos(2 * s - A1 - A2))
        maxL = max(maxL, abs(C.Lam(t1, A1)), abs(C.Lam(t2, A2)), abs(C.Lam(s, A1 + A2)))
        maxD = max(maxD, abs(Dl))
        cnt('(6) sampled point of {|Delta| < B(s)} has E >= K', C.E_of_t(A, t1, t2) >= K)
print('(6) on M\' (3600 sampled points, 300 shapes): min P = %s (claim >= 3.5), max |Lambda_1|,|Lambda_2|,|Lambda_3| = %s, '
      'max |Delta| = %s (log 2 = %s)' % (mp.nstr(minP, 6), mp.nstr(maxL, 6), mp.nstr(maxD, 6), mp.nstr(mp.log(2), 6)))

# (7) vertex-3 quadratic, compared with bisection
for k in range(400):
    A = rand_shape(0.3)
    A1, A2, A3 = A
    A12 = A1 + A2
    pm = mp.acos((2 - (F - 1) * mp.cos(A3)) / (F + 1))
    rho = mp.mpf(float(rng.uniform(-0.999, 0.999)))
    s = (A12 + rho * pm) / 2
    sign = 1 if k % 2 == 0 else -1
    Dl = sign * mp.acosh(beta_s(A, s))
    qa = mp.exp(-Dl) * mp.sin(A12 - s)
    qb = mp.sin(A1 - s) + mp.exp(-Dl) * mp.sin(A2 - s)
    qc = mp.sin(s)
    roots = mp.polyroots([qa, qb, -qc])
    pos = [mp.re(z) for z in roots if abs(mp.im(z)) < mp.mpf('1e-30') and mp.re(z) > 0]
    cnt('(7) number of positive roots != 1', len(pos) != 1)
    z = pos[0]
    t1 = mp.arg(1 + z * mp.expj(A1))
    t2 = mp.arg(1 + z * mp.exp(-Dl) * mp.expj(A2))
    upd('(7) |t1 + t2 - s|', t1 + t2 - s)
    upd('(7) |E(t) - K| / K', (C.E_of_t(A, t1, t2) - K) / K)
    upd('(7) |Lambda1 - Lambda2 - Delta|', C.Lam(t1, A1) - C.Lam(t2, A2) - Dl)
    t1b, t2b = t_of_sDelta(A, s, Dl)
    upd('(7) |t(quadratic) - t(bisection)|', max(abs(t1 - t1b), abs(t2 - t2b)))

# (8) vertex-1 quadratic
for k in range(400):
    A = rand_shape(0.3)
    A1, A2, A3 = A
    A12 = A1 + A2
    pm1 = mp.acos((2 + (F - 1) * mp.cos(A1)) / (F + 1))
    cnt('(8) psi1_max >= A1', pm1 >= A1)
    rho = mp.mpf(float(rng.uniform(-0.999, 0.999)))
    t1 = (A1 + rho * pm1) / 2
    beta1 = ((F + 1) * mp.cos(2 * t1 - A1) - (F - 1) * mp.cos(A1)) / 2
    sign = 1 if k % 2 == 0 else -1
    Dv = sign * mp.acosh(beta1)
    eD = mp.exp(Dv)
    roots = mp.polyroots([eD * mp.sin(A1 - t1), eD * mp.sin(A12 - t1) - mp.sin(A2 + t1), -mp.sin(t1)])
    pos = [mp.re(z) for z in roots if abs(mp.im(z)) < mp.mpf('1e-30') and mp.re(z) > 0]
    cnt('(8) number of positive roots != 1', len(pos) != 1)
    z = pos[0]
    t2 = mp.atan(z * mp.sin(A2) / (1 + z * mp.cos(A2)))
    cnt('(8) t2 not in (0, A2)', not (0 < t2 < A2))
    upd('(8) |E(t1, t2) - K| / K', (C.E_of_t(A, t1, t2) - K) / K)
    upd('(8) |Lambda3(s) - Lambda2(t2) - D|', C.Lam(t1 + t2, A12) - C.Lam(t2, A2) - Dv)

print('\nTABLE (max abs/rel errors, and counts of failures):')
for key in sorted(worst):
    v = worst[key]
    print('  %-90s %s' % (key, mp.nstr(v, 3) if not isinstance(v, int) else v))
