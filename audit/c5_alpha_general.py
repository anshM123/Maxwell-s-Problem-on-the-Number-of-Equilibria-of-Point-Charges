"""AUDIT_3, item 1: numerical checks (40 digits) of every exponent-dependent step of the reduction, at
alpha in {1/4, 0.7, 1, 3/2, 2}; geometry from first principles (own code; nothing from code/ or audit/b0_common.py).

 (1) critical points by construction: x in R+ (kappa > 0) or R- (kappa < 0), q_i = kappa lambda_i r_i^beta; then
     grad V(x) = 0, and Hess V (by numerical differentiation of V itself) = 2 alpha kappa (beta S - I);
 (2) type: extremum <=> D > c = b/(b+1)^2 (eigenvalues of the numerical Hessian), minimum iff kappa > 0; no extremum
     in R-; every extremum in R+ has P_x < 0 and t = (sigma1, sigma2) in Q;
 (3) separated charges: for x in R+^out, grad H(t(x)) = psi_q with b = 2 alpha + 1 (H's gradient from phi_i);
 (4) d = (phi1', phi2', -phi3') by numerical differentiation of phi_i, and e2(d) = (b+1)^2 - b E at random t in Q;
     det Hess H (numerical second derivatives of H) = e2(d);
 (5) Hess H negative definite at t0 = (A1/2, A2/2), and E(t0) = 4 < K;
 (6) Poincare-Hopf degrees for V|plane at general alpha: winding of grad V on small circles (1 for either sign) and on
     a large circle (1 if Q_tot != 0, 2 if Q_tot = 0);
 (7) planarity in R^4: Newton iterations for grad V = 0 in R^4 from off-plane seeds end in the plane.
Output: audit/c5_alpha_general.out"""
import random

import mpmath as mp

mp.mp.dps = 40
rnd = random.Random(55)


def tri(A1, A2):
    A3 = mp.pi - A1 - A2
    L = mp.sin(A2) / mp.sin(A3)          # |p1 p3| with |p1 p2| = 1
    return [mp.mpc(0, 0), mp.mpc(1, 0), mp.mpc(L * mp.cos(A1), L * mp.sin(A1))]


def cross(u, v):
    return (u.conjugate() * v).imag


def circ(P):
    a, b, c = P
    # circumcentre by solving |O - a| = |O - b| = |O - c|
    M = mp.matrix([[2 * (b.real - a.real), 2 * (b.imag - a.imag)], [2 * (c.real - a.real), 2 * (c.imag - a.imag)]])
    rhs = mp.matrix([abs(b) ** 2 - abs(a) ** 2, abs(c) ** 2 - abs(a) ** 2])
    o = mp.lu_solve(M, rhs)
    O = mp.mpc(o[0], o[1])
    return O, abs(a - O)


def geo(P, x):
    D2 = cross(P[1] - P[0], P[2] - P[0])
    lam, sig, r = [], [], []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        lam.append(cross(P[j] - x, P[k] - x) / D2)
        sig.append(mp.arg((P[k] - x) / (P[j] - x)))
        r.append(abs(x - P[i]))
    O, R = circ(P)
    S = mp.matrix(2, 2)
    for i in range(3):
        u = (x - P[i]) / r[i]
        S[0, 0] += lam[i] * u.real ** 2
        S[0, 1] += lam[i] * u.real * u.imag
        S[1, 1] += lam[i] * u.imag ** 2
    S[1, 0] = S[0, 1]
    return lam, sig, r, R ** 2 - abs(x - O) ** 2, S


def V(P, q, alpha, x, y):
    return sum(qi * ((x - p.real) ** 2 + (y - p.imag) ** 2) ** (-alpha) for qi, p in zip(q, P))


def grad_hess_num(P, q, alpha, x):
    f = lambda u, v: V(P, q, alpha, u, v)
    g = (mp.diff(f, (x.real, x.imag), (1, 0)), mp.diff(f, (x.real, x.imag), (0, 1)))
    H = mp.matrix([[mp.diff(f, (x.real, x.imag), (2, 0)), mp.diff(f, (x.real, x.imag), (1, 1))],
                   [mp.diff(f, (x.real, x.imag), (1, 1)), mp.diff(f, (x.real, x.imag), (0, 2))]])
    return g, H


def phi(t, a, b):
    return mp.log(mp.sin(t)) + b * mp.log(mp.sin(a - t))


def phi3(s, A3, b):
    return mp.log(mp.sin(s)) + b * mp.log(mp.sin(s + A3))


def gradH(A, t1, t2, b):
    A1, A2, A3 = A
    return phi(t1, A1, b) - phi3(t1 + t2, A3, b), phi(t2, A2, b) - phi3(t1 + t2, A3, b)


def E_of_t(A, t1, t2):
    A1, A2, _ = A
    U = lambda t, a: mp.cot(t) + mp.cot(a - t)
    U1, U2, U3 = U(t1, A1), U(t2, A2), U(t1 + t2, A1 + A2)
    return U1 * U2 - U3 * (U1 + U2)


def winding(P, q, alpha, centre, rad, n=720):
    tot = mp.mpf(0)
    prev = None
    for k in range(n + 1):
        th = 2 * mp.pi * k / n
        x = centre + rad * mp.expj(th)
        g = mp.mpc(0, 0)
        for qi, p in zip(q, P):
            d = x - p
            g += -2 * alpha * qi * abs(d) ** (-2 * alpha - 2) * d
        a = mp.arg(g)
        if prev is not None:
            da = a - prev
            while da > mp.pi:
                da -= 2 * mp.pi
            while da < -mp.pi:
                da += 2 * mp.pi
            tot += da
        prev = a
    return int(mp.nint(tot / (2 * mp.pi)))


def main():
    out = []
    alphas = [mp.mpf(1) / 4, mp.mpf('0.7'), mp.mpf(1), mp.mpf(3) / 2, mp.mpf(2)]
    for alpha in alphas:
        b = 2 * alpha + 1
        beta = b + 1
        c = b / (b + 1) ** 2
        K = (b + 1) ** 2 / b
        st = dict(grad=0, hess=0, type_fail=0, max_in_Rm=0, ext_not_out=0, psi=0, e2d=0, detH=0, t0_fail=0,
                  n_cp=0, n_ext=0, wind_small=0, wind_large_fail=0, n_wind=0)
        for trial in range(60):
            A1 = mp.mpf(rnd.uniform(0.05, 2.9))
            A2 = mp.mpf(rnd.uniform(0.02, float(mp.pi - A1) - 0.02))
            A = (A1, A2, mp.pi - A1 - A2)
            P = tri(A1, A2)
            # a point in R+ (lambda = (+,+,-)) or R- (lambda = (-,-,+)): x = sum mu_i p_i with prescribed sign pattern
            if trial % 4 == 0:
                mu = [mp.mpf(rnd.uniform(0.2, 3)), mp.mpf(rnd.uniform(0.2, 3))]
                mu3 = 1 - mu[0] - mu[1]
                kap = 1
            elif trial % 2 == 0:
                # near the excentre opposite p3 (barycentrics proportional to (sin A1, sin A2, -sin A3)), where D -> 1/4
                sA = [mp.sin(A[0]), mp.sin(A[1]), mp.sin(A[2])]
                den = sA[0] + sA[1] - sA[2]
                mu = [sA[0] / den * (1 + mp.mpf(rnd.uniform(-0.3, 0.3))), sA[1] / den * (1 + mp.mpf(rnd.uniform(-0.3, 0.3)))]
                mu3 = 1 - mu[0] - mu[1]
                kap = 1
            else:
                mu = [-mp.mpf(rnd.uniform(0.05, 2)), -mp.mpf(rnd.uniform(0.05, 2))]
                mu3 = 1 - mu[0] - mu[1]
                kap = -1
            x = mu[0] * P[0] + mu[1] * P[1] + mu3 * P[2]
            lam, sig, r, Px, S = geo(P, x)
            q = [kap * lam[i] * r[i] ** beta for i in range(3)]
            if not (q[0] > 0 and q[1] > 0 and q[2] < 0):
                st['rejected'] = st.get('rejected', 0) + 1      # construction outside R+/R- (not a test failure)
                continue
            g, H = grad_hess_num(P, q, alpha, x)
            scale = max(abs(qi) * ri ** (-2 * alpha - 1) for qi, ri in zip(q, r))
            st['grad'] = max(st['grad'], max(abs(g[0]), abs(g[1])) / scale)
            Hf = 2 * alpha * kap * (beta * S - mp.eye(2))
            st['hess'] = max(st['hess'], mp.mnorm(H - Hf, 1) / mp.mnorm(Hf, 1))
            ev = mp.eigsy(H)[0]
            D = mp.det(S)
            ext = ev[0] * ev[1] > 0
            st['n_cp'] += 1
            if ext != (D > c):
                st['type_fail'] += 1
            if ext:
                st['n_ext'] += 1
                if kap < 0:
                    st['max_in_Rm'] += 1
                if not (ev[0] > 0 and ev[1] > 0):
                    st['type_fail'] += 1
                t1, t2 = sig[0], sig[1]
                if not (Px < 0 and 0 < t1 < A1 and 0 < t2 < A2):
                    st['ext_not_out'] += 1
            if kap > 0 and Px < 0:
                t1, t2 = sig[0], sig[1]
                gH = gradH(A, t1, t2, b)
                psi = (mp.log(q[0] / -q[2]) + b * mp.log(mp.sin(A1) / mp.sin(A[2])),
                       mp.log(q[1] / -q[2]) + b * mp.log(mp.sin(A2) / mp.sin(A[2])))
                st['psi'] = max(st['psi'], abs(gH[0] - psi[0]) + abs(gH[1] - psi[1]))
            # (4) e2(d) and det Hess H at a random t in Q
            t1, t2 = A1 * mp.mpf(rnd.uniform(0.02, 0.98)), A2 * mp.mpf(rnd.uniform(0.02, 0.98))
            d1 = mp.diff(lambda u: phi(u, A1, b), t1)
            d2 = mp.diff(lambda u: phi(u, A2, b), t2)
            d3 = -mp.diff(lambda u: phi3(u, A[2], b), t1 + t2)
            e2 = d1 * d2 + d2 * d3 + d3 * d1
            E = E_of_t(A, t1, t2)
            st['e2d'] = max(st['e2d'], abs(e2 - ((b + 1) ** 2 - b * E)) / abs(e2))
            Hh = mp.matrix(2, 2)
            for i in range(2):
                for j in range(2):
                    Hh[i, j] = mp.diff(lambda u, v: gradH(A, u, v, b)[i], (t1, t2), (1 if j == 0 else 0, 1 if j == 1 else 0))
            st['detH'] = max(st['detH'], abs(mp.det(Hh) - e2) / abs(e2))
            # (5) negative definiteness at t0, E(t0) = 4 < K
            u1, u2 = A1 / 2, A2 / 2
            H0 = mp.matrix(2, 2)
            for i in range(2):
                for j in range(2):
                    H0[i, j] = mp.diff(lambda u, v: gradH(A, u, v, b)[i], (u1, u2), (1 if j == 0 else 0, 1 if j == 1 else 0))
            H0 = (H0 + H0.T) / 2
            e0 = mp.eigsy(H0)[0]
            if not (e0[0] < 0 and e0[1] < 0 and abs(E_of_t(A, u1, u2) - 4) < mp.mpf(10) ** -30 and K > 4):
                st['t0_fail'] += 1
            # (6) winding numbers (every 6th trial; one configuration with Q_tot = 0 by rescaling q3)
            if trial % 6 == 0:
                for qq in (q, [q[0], q[1], -(q[0] + q[1])]):
                    Qt = sum(qq)
                    ws = [winding(P, qq, alpha, P[i], mp.mpf('1e-3') * min(abs(P[i] - P[j]) for j in range(3) if j != i),
                                  240) for i in range(3)]
                    st['wind_small'] += sum(1 for w in ws if w != 1)
                    big = winding(P, qq, alpha, mp.mpc(0, 0), mp.mpf(10) ** 6, 360)
                    want = 2 if abs(Qt) < mp.mpf(10) ** -30 else 1
                    st['n_wind'] += 1
                    if big != want:
                        st['wind_large_fail'] += 1
                        out.append('   large-circle winding %d (want %d) at Q_tot = %s' % (big, want, mp.nstr(Qt, 5)))
        out.append('alpha %-5s (b %-4s K %s): %d critical points built (%d extrema); max |grad V|/scale %s; '
                   'Hess V vs 2 alpha kappa (beta S - I) rel %s; type/sign failures %d; extrema in R- %d; '
                   'extrema not in R+out with t in Q %d; |grad H - psi_q| %s; e2(d) rel %s; det Hess H vs e2(d) rel %s; '
                   't0 failures %d; small-circle winding != 1: %d; large-circle failures %d of %d; rejected constructions %d'
                   % (mp.nstr(alpha, 3), mp.nstr(b, 3), mp.nstr(K, 5), st['n_cp'], st['n_ext'], mp.nstr(st['grad'], 2),
                      mp.nstr(st['hess'], 2), st['type_fail'], st['max_in_Rm'], st['ext_not_out'],
                      mp.nstr(st['psi'], 2), mp.nstr(st['e2d'], 2), mp.nstr(st['detH'], 2), st['t0_fail'],
                      st['wind_small'], st['wind_large_fail'], st['n_wind'], st.get('rejected', 0)))
        print(out[-1], flush=True)
    # (7) planarity in R^4 (float Newton is enough to see where the iterations end)
    import numpy as np
    rng = np.random.default_rng(7)
    worst_h, n_conv = 0.0, 0
    for alpha in (0.25, 1.0, 2.0):
        for _ in range(40):
            p = np.zeros((3, 4))
            p[:, :2] = rng.normal(size=(3, 2))
            q = np.array([rng.uniform(0.2, 2), rng.uniform(0.2, 2), -rng.uniform(0.2, 2)])
            rng.shuffle(q)
            beta = 2 * alpha + 2
            for _ in range(30):
                x = np.concatenate([rng.normal(size=2) * 2, rng.normal(size=2) * 0.5])
                for it in range(200):
                    d = x[None, :] - p
                    r = np.linalg.norm(d, axis=1)
                    w = q * r ** (-beta)
                    g = -2 * alpha * (w[:, None] * d).sum(0)
                    H = -2 * alpha * sum(w[i] * (np.eye(4) - beta * np.outer(d[i], d[i]) / r[i] ** 2) for i in range(3))
                    try:
                        step = np.linalg.solve(H, g)
                    except np.linalg.LinAlgError:
                        break
                    nrm = np.linalg.norm(step)
                    if nrm > 0.5:
                        step *= 0.5 / nrm
                    x = x - step
                    if nrm < 1e-14 * (1 + np.linalg.norm(x)):
                        break
                d = x[None, :] - p
                r = np.linalg.norm(d, axis=1)
                g = -2 * alpha * ((q * r ** (-beta))[:, None] * d).sum(0)
                if np.linalg.norm(g) < 1e-9 * np.abs(q * r ** (-beta - 1)).max() and r.min() > 1e-6:
                    n_conv += 1
                    worst_h = max(worst_h, np.linalg.norm(x[2:]) / (1 + np.linalg.norm(x[:2])))
    out.append('planarity in R^4 (alpha 0.25, 1, 2; 120 configurations, 3600 seeds): %d converged critical points, '
               'max off-plane component %.1e' % (n_conv, worst_h))
    print(out[-1])


if __name__ == '__main__':
    main()
