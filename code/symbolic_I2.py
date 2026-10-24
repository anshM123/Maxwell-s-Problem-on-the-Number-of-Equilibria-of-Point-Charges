"""Exact proof of identity (I2) by polynomial expansion (no simplify):

  det DPsi = beta^2 (D - c) / (2 Delta lambda1 lambda2 lambda3),   c = (beta - 1)/beta^2,

for the triangle p1 = (0,0), p2 = (1,0), p3 = (a,b) and x = (X,Y), with
  DPsi rows = grad L_i - grad L_3 (i = 1, 2),   grad L_i = g_i / lambda_i + beta d_i / r_i^2,
  g_i = grad lambda_i (constant), d_i = x - p_i, r_i^2 = |d_i|^2,
  D = sum_{i<j} lambda_i lambda_j (d_i x d_j)^2 / (r_i^2 r_j^2),   Delta = b/2.
Both sides are multiplied by the nonzero factor W = lambda1 lambda2 lambda3 r1^2 r2^2 r3^2 (times b for the
barycentric denominators), and the difference of the two polynomials is expanded and must be identically zero."""
import sympy as sp

a, b, X, Y, beta = sp.symbols('a b X Y beta')


def cross(u, v):
    return u[0] * v[1] - u[1] * v[0]


P = [(0, 0), (1, 0), (a, b)]
x = (X, Y)
d = [(X - p[0], Y - p[1]) for p in P]
r2 = [di[0] ** 2 + di[1] ** 2 for di in d]
# b * lambda_i are polynomials (signed areas times 2); write lambda_i = l_i / b
l = [cross((P[1][0] - X, P[1][1] - Y), (P[2][0] - X, P[2][1] - Y)),
     cross((P[2][0] - X, P[2][1] - Y), (P[0][0] - X, P[0][1] - Y)),
     cross((P[0][0] - X, P[0][1] - Y), (P[1][0] - X, P[1][1] - Y))]
assert sp.expand(sum(l) - b) == 0
g = [(sp.diff(li, X), sp.diff(li, Y)) for li in l]            # gradients of b*lambda_i, so g_i/lambda_i = g~_i / l_i
# det DPsi = sum_cyc cross(G_i, G_j), G_i = g~_i / l_i + beta d_i / r_i^2
# multiply by W = l1 l2 l3 r1^2 r2^2 r3^2: G_i * (l_i r_i^2) = g~_i r_i^2 + beta l_i d_i  =: N_i  (polynomial vectors)
N = [(g[i][0] * r2[i] + beta * l[i] * d[i][0], g[i][1] * r2[i] + beta * l[i] * d[i][1]) for i in range(3)]
# cross(G_i, G_j) * (l_i r_i^2)(l_j r_j^2) = cross(N_i, N_j); the cyclic sum times W needs factors l_k r_k^2
lhs = sum(cross(N[i], N[j]) * l[k] * r2[k] for (i, j, k) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)))
# right side times W: beta^2 (D - c)/(2 Delta lambda1 lambda2 lambda3) * W, with lambda = l/b, Delta = b/2:
#   = beta^2 (D - c) * b^3/(b l1 l2 l3) * l1 l2 l3 r1^2 r2^2 r3^2 = beta^2 b^2 (D - c) r1^2 r2^2 r3^2
#   D r1^2 r2^2 r3^2 = sum_{i<j} (l_i l_j / b^2) cross(d_i, d_j)^2 r_k^2
Dterm = sum(l[i] * l[j] * cross(d[i], d[j]) ** 2 * r2[k] for (i, j, k) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)))
rhs = beta ** 2 * Dterm - (beta - 1) * b ** 2 * r2[0] * r2[1] * r2[2]
diff = sp.expand(lhs - rhs)
print('I2 identity, expanded difference:', diff)
print('number of terms in lhs:', len(sp.Add.make_args(sp.expand(lhs))))
