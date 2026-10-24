"""AUDIT of the REWRITTEN PROOF.md (2026-10-24 10:01): independent checks of Lemmas 3-6 and Section 7 identities.
Nothing is imported from code.  40-digit mpmath unless stated.

 (L3) pedal identity and charge formula, for x in T (acute and obtuse triangles, points near edges/vertices):
      delta_i = theta_i - A_i > 0, sum delta = pi; r_j r_k sin delta_i = P_x sin A_i (P_x = R^2 - |x-O|^2);
      lambda_i r_i^beta = C(x) sin theta_i sin^b delta_i / sin^b A_i, C = (r1 r2 r3)^(b+1)/(2 Delta P_x^b);
      inversion facts used in the proof: the angles of T* = (p_i*) are delta_i, circumradius of T* = R/P_x.
 (L4) W_i = sin A_i/(sin delta_i sin theta_i) = cot delta_i - cot theta_i; Sum lambda_k r_k^2 = P_x (Leibniz);
      D (computed as det S from the vectors) = 1/E, E = W1W2 + W2W3 + W3W1.
 (L5) on all of H_A (theta NOT required to come from a point of T): l_i'' = csc^2(t - A_i) + csc^2 t; the second
      derivative of E along random directions v (sum v = 0) is positive and >= sum_pairs W_iW_j (l_i'' v_i^2 + l_j'' v_j^2).
 (L6) on H_A: e2(cot delta) = 1, e2(cot theta) = 1, e2(d) = (b+1)^2 - b E with d_i = cot theta_i + b cot delta_i;
      Gram determinant in the basis (e1-e3, e2-e3) = e2(d); incentre: theta = ((pi + A_i)/2), W_i = 2 tan(A_i/2), E = 4,
      D = 1/4, d_i = (b-1) tan(A_i/2).
 (S1) log Q_i(x) - f_i(theta_i) + b log sin A_i - log C(x) = 0 for each i (the mod-(1,1,1) bookkeeping).
 (S3) at critical points (q = Q(x)): Hess V_q = 2 alpha kappa L o DPsi, kappa = q_i r_i^-beta / lambda_i,
      L(e) = sum lambda_i e_i (x - p_i); L has kernel exactly R(1,1,1).
usage: python a8_new_identities.py"""
import math

import mpmath as mp
import numpy as np

mp.mp.dps = 40
rng = np.random.default_rng(424242)
worst = {'L3a: min(delta) < 0 count': 0, 'S1/L6: sign(D - c) != sign(e2(d)) count': 0,
         'H_A L5: E'' <= 0 count': 0, 'H_A L5: E'' < proof lower bound count': 0}


def upd(key, val):
    worst[key] = max(worst.get(key, mp.mpf(0)), abs(val))


def rand_angles(rng, min_deg=0.5):
    while True:
        w = rng.dirichlet(np.ones(3)) * 180
        if w.min() > min_deg:
            return [mp.mpf(float(v)) * mp.pi / 180 for v in w[:2]]


def triangle(A1, A2):
    A3 = mp.pi - A1 - A2
    L13 = mp.sin(A2) / mp.sin(A3)
    return [mp.matrix([0, 0]), mp.matrix([1, 0]), mp.matrix([L13 * mp.cos(A1), L13 * mp.sin(A1)])], [A1, A2, A3]


cross = lambda u, v: u[0] * v[1] - u[1] * v[0]
dot = lambda u, v: u[0] * v[0] + u[1] * v[1]
norm = lambda u: mp.sqrt(dot(u, u))


def angle(o, a, b):
    u, v = a - o, b - o
    return mp.acos(dot(u, v) / (norm(u) * norm(v)))


def circum(P):
    (ax, ay), (bx, by), (cx, cy) = [(p[0], p[1]) for p in P]
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    ux = ((ax ** 2 + ay ** 2) * (by - cy) + (bx ** 2 + by ** 2) * (cy - ay) + (cx ** 2 + cy ** 2) * (ay - by)) / d
    uy = ((ax ** 2 + ay ** 2) * (cx - bx) + (bx ** 2 + by ** 2) * (ax - cx) + (cx ** 2 + cy ** 2) * (bx - ax)) / d
    O = mp.matrix([ux, uy])
    return O, norm(P[0] - O)


def quantities(P, A, x, alpha):
    b = 2 * alpha + 1
    beta = b + 1
    Delta2 = cross(P[1] - P[0], P[2] - P[0])                       # 2 Delta (ccw)
    lam = [cross(P[(i + 1) % 3] - x, P[(i + 2) % 3] - x) / Delta2 for i in range(3)]
    r = [norm(x - P[i]) for i in range(3)]
    th = [angle(x, P[(i + 1) % 3], P[(i + 2) % 3]) for i in range(3)]
    de = [th[i] - A[i] for i in range(3)]
    O, R = circum(P)
    Px = R ** 2 - norm(x - O) ** 2
    return b, beta, Delta2, lam, r, th, de, O, R, Px


def e2(v):
    return v[0] * v[1] + v[1] * v[2] + v[2] * v[0]


# ---------------------------------------------------------------- (L3), (L4), (L6 on theta(T)), (S1)
npts = 0
for k in range(400):
    A1, A2 = rand_angles(rng, 0.5 if k % 2 else 10)
    P, A = triangle(A1, A2)
    w = rng.dirichlet(np.ones(3))
    if k % 4 == 2:
        w[rng.integers(3)] = 1e-6; w /= w.sum()
    if k % 4 == 3:
        w = np.full(3, 1e-6); w[rng.integers(3)] = 1; w /= w.sum()
    x = sum((mp.mpf(float(w[i])) * P[i] for i in range(3)), mp.matrix([0, 0]))
    alpha = mp.mpf(float([0.5, 0.05, 0.25, 1.0, 2.5, 10.0][k % 6]))
    b, beta, Delta2, lam, r, th, de, O, R, Px = quantities(P, A, x, alpha)
    npts += 1
    upd('L3a: min(delta) < 0 count', 1 if min(de) <= 0 else 0)
    upd('L3a: |sum delta - pi|', sum(de) - mp.pi)
    for i in range(3):
        j, kk = (i + 1) % 3, (i + 2) % 3
        upd('L3b: rel |r_j r_k sin delta_i - P_x sin A_i|', (r[j] * r[kk] * mp.sin(de[i]) - Px * mp.sin(A[i])) / (Px * mp.sin(A[i])))
    C = (r[0] * r[1] * r[2]) ** (b + 1) / (Delta2 * Px ** b)
    for i in range(3):
        lhs = lam[i] * r[i] ** beta
        rhs = C * mp.sin(th[i]) * mp.sin(de[i]) ** b / mp.sin(A[i]) ** b
        upd('L3c: rel |lambda_i r_i^beta - C sin th sin^b de / sin^b A|', lhs / rhs - 1)
    # inversion facts
    Pst = [x + (P[i] - x) / r[i] ** 2 for i in range(3)]
    for i in range(3):
        upd('L3 proof: |angle of T* at p_i* - delta_i|', angle(Pst[i], Pst[(i + 1) % 3], Pst[(i + 2) % 3]) - de[i])
    Ost, Rst = circum(Pst)
    upd('L3 proof: rel |R* - R/P_x|', Rst / (R / Px) - 1)
    # Lemma 4
    W = [mp.sin(A[i]) / (mp.sin(de[i]) * mp.sin(th[i])) for i in range(3)]
    for i in range(3):
        upd('L4: |W_i - (cot delta_i - cot theta_i)| / W_i', (W[i] - (mp.cot(de[i]) - mp.cot(th[i]))) / W[i])
    upd('L4: rel |sum lambda_k r_k^2 - P_x|', (sum(lam[i] * r[i] ** 2 for i in range(3)) - Px) / Px)
    S = mp.matrix(2, 2)
    for i in range(3):
        u = (x - P[i]) / r[i]
        S += lam[i] * (u * u.T)
    Dv = mp.det(S)
    E = e2(W)
    upd('L4: |D E - 1|', Dv * E - 1)
    # Lemma 6 identities on theta(T)
    X = [mp.cot(v) for v in de]
    Cc = [mp.cot(v) for v in th]
    d = [Cc[i] + b * X[i] for i in range(3)]
    upd('L6: |e2(cot delta) - 1|', e2(X) - 1)
    upd('L6: |e2(cot theta) - 1|', e2(Cc) - 1)
    upd('L6: rel |e2(d) - ((b+1)^2 - b E)|', (e2(d) - ((b + 1) ** 2 - b * E)) / ((b + 1) ** 2 + b * E))
    # sign(D - c) == sign(e2(d))
    c = b / (b + 1) ** 2
    if (Dv - c) * e2(d) < 0:
        upd('S1/L6: sign(D - c) != sign(e2(d)) count', 1)
    # Step 1 bookkeeping
    for i in range(3):
        fi = mp.log(mp.sin(th[i])) + b * mp.log(mp.sin(th[i] - A[i]))
        upd('S1: |log Q_i - f_i(theta_i) + b log sin A_i - log C|',
            mp.log(lam[i] * r[i] ** beta) - fi + b * mp.log(mp.sin(A[i])) - mp.log(C))

# ---------------------------------------------------------------- (L5), (L6) on all of H_A
nH = 0
for k in range(400):
    A1, A2 = rand_angles(rng, 0.5 if k % 2 else 10)
    A = [A1, A2, mp.pi - A1 - A2]
    b = 2 * mp.mpf(float([0.5, 0.05, 0.25, 1.0, 2.5, 10.0][k % 6])) + 1
    while True:                                    # uniform-ish random point of the hexagon H_A
        t12 = [A[i] + (mp.pi - A[i]) * mp.mpf(float(rng.uniform())) for i in range(2)]
        t3 = 2 * mp.pi - t12[0] - t12[1]
        if A[2] < t3 < mp.pi:
            break
    th = [t12[0], t12[1], t3]
    nH += 1
    de = [th[i] - A[i] for i in range(3)]
    X = [mp.cot(v) for v in de]
    Cc = [mp.cot(v) for v in th]
    W = [X[i] - Cc[i] for i in range(3)]
    E = e2(W)
    d = [Cc[i] + b * X[i] for i in range(3)]
    upd('H_A L6: |e2(cot delta) - 1|', e2(X) - 1)
    upd('H_A L6: |e2(cot theta) - 1|', e2(Cc) - 1)
    upd('H_A L6: rel |e2(d) - ((b+1)^2 - bE)|', (e2(d) - ((b + 1) ** 2 - b * E)) / ((b + 1) ** 2 + b * E))
    gram = mp.matrix([[d[0] + d[2], d[2]], [d[2], d[1] + d[2]]])
    upd('H_A L6: |det Gram - e2(d)| (rel)', (mp.det(gram) - e2(d)) / (abs(e2(d)) + 1))
    # l_i''
    for i in range(3):
        li = lambda t, i=i: mp.log(mp.sin(A[i])) - mp.log(mp.sin(t - A[i])) - mp.log(mp.sin(t))
        upd('H_A L5: rel |l_i\'\' - (csc^2(t-A_i) + csc^2 t)|',
            mp.diff(li, th[i], 2) / (mp.csc(th[i] - A[i]) ** 2 + mp.csc(th[i]) ** 2) - 1)
    # E'' along a random direction v (sum v = 0), and the lower bound of the proof
    v = rng.normal(size=3); v -= v.mean(); v = [mp.mpf(float(t)) for t in v]

    def Eline(s):
        tt = [th[i] + s * v[i] for i in range(3)]
        return e2([mp.sin(A[i]) / (mp.sin(tt[i] - A[i]) * mp.sin(tt[i])) for i in range(3)])

    E2 = mp.diff(Eline, 0, 2)
    lpp = [mp.csc(de[i]) ** 2 + mp.csc(th[i]) ** 2 for i in range(3)]
    lower = sum(W[i] * W[j] * (lpp[i] * v[i] ** 2 + lpp[j] * v[j] ** 2) for i, j in ((0, 1), (1, 2), (2, 0)))
    if E2 <= 0:
        upd('H_A L5: E\'\' <= 0 count', 1)
    if E2 < lower * (1 - mp.mpf('1e-25')):
        upd('H_A L5: E\'\' < proof lower bound count', 1)
    upd('H_A L5: min E\'\'/lower bound (reported as 1/x)', 1 / (E2 / lower))

# ---------------------------------------------------------------- incentre facts
for k in range(100):
    A1, A2 = rand_angles(rng, 0.5)
    P, A = triangle(A1, A2)
    a = [norm(P[(i + 1) % 3] - P[(i + 2) % 3]) for i in range(3)]
    I = sum((a[i] * P[i] for i in range(3)), mp.matrix([0, 0])) / sum(a)
    alpha = mp.mpf('0.5')
    b, beta, Delta2, lam, r, th, de, O, R, Px = quantities(P, A, I, alpha)
    W = [mp.sin(A[i]) / (mp.sin(de[i]) * mp.sin(th[i])) for i in range(3)]
    for i in range(3):
        upd('incentre: |theta_i - (pi + A_i)/2|', th[i] - (mp.pi + A[i]) / 2)
        upd('incentre: |W_i - 2 tan(A_i/2)|', W[i] - 2 * mp.tan(A[i] / 2))
        upd('incentre: |d_i - (b-1) tan(A_i/2)|', mp.cot(th[i]) + b * mp.cot(de[i]) - (b - 1) * mp.tan(A[i] / 2))
    upd('incentre: |E - 4|', e2(W) - 4)

# ---------------------------------------------------------------- (S3) Hess V_q = 2 alpha kappa L o DPsi
mp.mp.dps = 30
for k in range(120):
    A1, A2 = rand_angles(rng, 2)
    P, A = triangle(A1, A2)
    w = rng.dirichlet(np.ones(3))
    x0 = sum((mp.mpf(float(w[i])) * P[i] for i in range(3)), mp.matrix([0, 0]))
    alpha = mp.mpf(float([0.5, 0.25, 1.0, 2.5][k % 4]))
    beta = 2 * alpha + 2
    Delta2 = cross(P[1] - P[0], P[2] - P[0])

    def logQ(xv, yv, i):
        xx = mp.matrix([xv, yv])
        lam_i = cross(P[(i + 1) % 3] - xx, P[(i + 2) % 3] - xx) / Delta2
        return mp.log(lam_i) + beta * mp.log(norm(xx - P[i]))

    lam = [cross(P[(i + 1) % 3] - x0, P[(i + 2) % 3] - x0) / Delta2 for i in range(3)]
    r = [norm(x0 - P[i]) for i in range(3)]
    q = [lam[i] * r[i] ** beta for i in range(3)]
    kappa = q[0] * r[0] ** (-beta) / lam[0]
    H = mp.matrix(2, 2)
    for i in range(3):
        u = (x0 - P[i]) / r[i]
        H += -2 * alpha * q[i] * r[i] ** (-beta) * (mp.eye(2) - beta * (u * u.T))
    LD = mp.matrix(2, 2)
    for i in range(3):
        g = [mp.diff(lambda s: logQ(s, x0[1], i), x0[0]), mp.diff(lambda s: logQ(x0[0], s, i), x0[1])]
        for a_ in range(2):
            for c_ in range(2):
                LD[a_, c_] += lam[i] * (x0[a_] - P[i][a_]) * g[c_]
    M = 2 * alpha * kappa * LD
    upd('S3: max |Hess V_q - 2 alpha kappa L DPsi| / |Hess|', max(abs(H[i, j] - M[i, j]) for i in range(2) for j in range(2))
        / max(abs(H[i, j]) for i in range(2) for j in range(2)))
    Lmat = mp.matrix(2, 3)
    for i in range(3):
        for a_ in range(2):
            Lmat[a_, i] = lam[i] * (x0[a_] - P[i][a_])
    sv = mp.svd_r(Lmat, compute_uv=False)
    upd('S3: L(1,1,1) (abs)', max(abs(sum(Lmat[a_, i] for i in range(3))) for a_ in range(2)))
    upd('S3: 1/(second singular value of L) (L has rank 2)', 1 / min(sv[0], sv[1]))

print('%d points of T (half with an angle down to 0.5 deg; a quarter near an edge, a quarter near a vertex; '
      'alpha in {0.05,0.25,0.5,1,2.5,10}); %d random points of H_A; 100 incentres; 120 critical points for S3' % (npts, nH))
for key in sorted(worst):
    print('  %-75s %s' % (key, mp.nstr(worst[key], 3)))
