"""Exact symbolic proofs, valid for EVERY point x off the side lines (inside or outside the triangle), of the signed
identities used for charges of mixed signs (MIXED.md sections 2-4).

Triangle p1 = (0,0), p2 = (1,0), p3 = (u,v) with v > 0 (counter-clockwise); x = (X,Y).  Signed angles at x:
  sin(sigma_i) = cross(p_j - x, p_k - x)/(r_j r_k),  cos(sigma_i) = dot(p_j - x, p_k - x)/(r_j r_k)   ((i,j,k) cyclic),
  cos A_i = dot(p_j - p_i, p_k - p_i)/(a_k a_j),  sin A_i = cross(p_j - p_i, p_k - p_i)/(a_k a_j).
Identities:
  (P)  r_j r_k sin(sigma_i - A_i) = P_x sin A_i,  P_x = R^2 - |x - O|^2      (signed pedal identity)
  (D)  D(x) = det sum_i lambda_i u_i u_i^T = P_x sin(sigma_1) sin(sigma_2) sin(sigma_3) / (2 Delta)
  (E)  D(x) * e2(W) = 1,  W_i = cot(sigma_i - A_i) - cot(sigma_i)
Square roots (r_i and the side lengths a_i) are carried as independent symbols, together with their squares.
"""
import sympy as sp

X, Y, u, v = sp.symbols('X Y u v', real=True)
P = [sp.Matrix([0, 0]), sp.Matrix([1, 0]), sp.Matrix([u, v])]
x = sp.Matrix([X, Y])
r = sp.symbols('r1 r2 r3', positive=True)            # r_i = |x - p_i|
al = sp.symbols('a1 a2 a3', positive=True)            # a_i = |p_j - p_k| (side opposite p_i)


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


r2 = [dot(x - p, x - p) for p in P]
a2 = [dot(P[(i + 1) % 3] - P[(i + 2) % 3], P[(i + 1) % 3] - P[(i + 2) % 3]) for i in range(3)]
subs_sq = {r[i] ** 2: r2[i] for i in range(3)}
subs_sq.update({al[i] ** 2: a2[i] for i in range(3)})


def reduce(expr):
    """replace even powers of the root symbols by polynomials, repeatedly"""
    expr = sp.expand(expr)
    for _ in range(6):
        new = expr
        for s, val in subs_sq.items():
            base = s.base
            new = new.replace(lambda e: e.is_Pow and e.base == base and e.exp.is_Integer and e.exp >= 2,
                              lambda e: val ** (e.exp // 2) * base ** (e.exp % 2))
        new = sp.expand(new)
        if new == expr:
            break
        expr = new
    return expr


Delta2 = cross(P[1] - P[0], P[2] - P[0])          # 2 Delta (> 0 for counter-clockwise labelling)
# circumcentre and P_x
ox, oy = sp.symbols('ox oy')
sol = sp.solve([dot(sp.Matrix([ox, oy]) - P[0], sp.Matrix([ox, oy]) - P[0]) - dot(sp.Matrix([ox, oy]) - P[1], sp.Matrix([ox, oy]) - P[1]),
                dot(sp.Matrix([ox, oy]) - P[0], sp.Matrix([ox, oy]) - P[0]) - dot(sp.Matrix([ox, oy]) - P[2], sp.Matrix([ox, oy]) - P[2])], [ox, oy])
O = sp.Matrix([sol[ox], sol[oy]])
R2 = dot(O - P[0], O - P[0])
Px = sp.simplify(R2 - dot(x - O, x - O))

ok = True
sins, coss, cosA, sinA = [], [], [], []
for i in range(3):
    j, k = (i + 1) % 3, (i + 2) % 3
    sins.append(cross(P[j] - x, P[k] - x) / (r[j] * r[k]))
    coss.append(dot(P[j] - x, P[k] - x) / (r[j] * r[k]))
    cosA.append(dot(P[j] - P[i], P[k] - P[i]) / (al[k] * al[j]))
    sinA.append(cross(P[j] - P[i], P[k] - P[i]) / (al[k] * al[j]))
# (P): r_j r_k (sin s cos A - cos s sin A) - P_x sin A, times a_j a_k
for i in range(3):
    j, k = (i + 1) % 3, (i + 2) % 3
    expr = (r[j] * r[k] * (sins[i] * cosA[i] - coss[i] * sinA[i]) - Px * sinA[i]) * al[j] * al[k]
    d = sp.simplify(reduce(sp.together(expr)))
    print('(P) i=%d: difference = %s' % (i + 1, d))
    ok &= d == 0
# (D): det S with signed barycentrics lambda_i = cross(p_j - x, p_k - x)/(2 Delta)
lam = [cross(P[(i + 1) % 3] - x, P[(i + 2) % 3] - x) / Delta2 for i in range(3)]
S = sp.zeros(2, 2)
for i in range(3):
    ui = (x - P[i]) / r[i]
    S += lam[i] * ui * ui.T
detS = reduce(sp.together(S.det()) * 1)
rhsD = Px * sins[0] * sins[1] * sins[2] / Delta2
d = sp.simplify(reduce(sp.together(detS - rhsD) * (r[0] * r[1] * r[2]) ** 2))
print('(D) difference =', d)
ok &= d == 0
# (E): D * e2(W) = 1 with W_i = cot(sigma_i - A_i) - cot sigma_i = sin A_i / (sin(sigma_i - A_i) sin sigma_i)
W = []
for i in range(3):
    s_minus = sins[i] * cosA[i] - coss[i] * sinA[i]
    W.append(sinA[i] / (s_minus * sins[i]))
e2W = W[0] * W[1] + W[1] * W[2] + W[2] * W[0]
num, den = sp.fraction(sp.together(rhsD * e2W - 1))
d = sp.simplify(reduce(sp.expand(num)))
print('(E) numerator of D e2(W) - 1 =', d)
ok &= d == 0
print('ALL SIGNED IDENTITIES PROVED' if ok else 'CHECK FAILED')
