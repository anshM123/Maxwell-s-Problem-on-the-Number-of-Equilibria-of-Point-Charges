"""Exact symbolic proof of the vertex-3 product representation used for mixed signs (MIXED.md section 6.1):

    E * sin(s) * sin(s + A3) = 2 (cosh Delta + cos A3),
    E = U1 U2 - U3 (U1 + U2),   U1 = cot t1 + cot(A1 - t1),   U2 = cot t2 + cot(A2 - t2),
    U3 = cot s + cot(A1 + A2 - s),   s = t1 + t2,   A3 = pi - A1 - A2,
    e^Delta = sin t1 sin(A2 - t2) / (sin(A1 - t1) sin t2).

Every sine and cosine is written through unit complex numbers x = e^{i t1}, y = e^{i t2}, a = e^{i A1}, b = e^{i A2}
(so e^{i A3} = -1/(a b)). Both sides then become rational functions, and their difference cancels to 0.
The vertex-1 representation E * sin t1 * sin(A1 - t1) = 2 (cosh(Lambda2 - Lambda3) - cos A1) is checked the same way.
"""
import sympy as sp

x, y, a, b = sp.symbols('x y a b', nonzero=True)


def S(z):          # sin of the angle whose exponential is z
    return (z - 1 / z) / (2 * sp.I)


def C(z):
    return (z + 1 / z) / 2


def cot(z):
    return C(z) / S(z)


ez1, ez2 = x, y
ezs = x * y                      # e^{i s}
eA1, eA2 = a, b
eA12 = a * b
eA3 = -1 / (a * b)               # e^{i (pi - A1 - A2)}
U1 = cot(ez1) + cot(eA1 / ez1)
U2 = cot(ez2) + cot(eA2 / ez2)
U3 = cot(ezs) + cot(eA12 / ezs)
E = U1 * U2 - U3 * (U1 + U2)

# vertex 3
eD = S(ez1) * S(eA2 / ez2) / (S(eA1 / ez1) * S(ez2))
coshD = (eD + 1 / eD) / 2
lhs = E * S(ezs) * S(ezs * eA3)
rhs = 2 * (coshD + C(eA3))
print('vertex 3: difference =', sp.simplify(sp.together(lhs - rhs)))

# vertex 1: Lambda2 - Lambda3 with Lambda2 = log(sin t2/sin(A2 - t2)), Lambda3 = log(sin s/sin(A12 - s))
eD1 = (S(ez2) / S(eA2 / ez2)) / (S(ezs) / S(eA12 / ezs))
coshD1 = (eD1 + 1 / eD1) / 2
lhs1 = E * S(ez1) * S(eA1 / ez1)
rhs1 = 2 * (coshD1 - C(eA1))
print('vertex 1: difference =', sp.simplify(sp.together(lhs1 - rhs1)))
