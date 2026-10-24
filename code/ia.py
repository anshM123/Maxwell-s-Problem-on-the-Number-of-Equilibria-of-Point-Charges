"""Vectorised interval arithmetic in IEEE binary64 with outward rounding.

An interval array is a pair (lo, hi) of float64 arrays.
- +, -, *, / and sqrt are correctly rounded by IEEE; each result is widened by two ulps (`dn`/`up`).
- Library transcendentals (sin, cos, log, exp, sinh, cosh, arcsin, arccos, arcsinh) are assumed accurate to within
  2 ulps, the same assumption as astro/orbit_distance_P8/FLOAT_RIGOUR.md. Each one is widened by four ulps (`dn4`/`up4`).
- Monotone pieces are evaluated at the endpoints. sin and cos include interior extrema explicitly.
"""
import numpy as np

INF = np.inf


def dn(x):
    return np.nextafter(np.nextafter(x, -INF), -INF)


def up(x):
    return np.nextafter(np.nextafter(x, INF), INF)


def dn4(x):
    return dn(dn(x))


def up4(x):
    return up(up(x))


def const(c, shape=None):
    lo, hi = dn(np.float64(c)), up(np.float64(c))
    if shape is not None:
        lo, hi = np.full(shape, lo), np.full(shape, hi)
    return lo, hi


PI = const(np.pi)


def iv(lo, hi):
    return np.asarray(lo, float), np.asarray(hi, float)


def add(x, y):
    return dn(x[0] + y[0]), up(x[1] + y[1])


def sub(x, y):
    return dn(x[0] - y[1]), up(x[1] - y[0])


def neg(x):
    return -x[1], -x[0]


def mul(x, y):
    p = np.stack([x[0] * y[0], x[0] * y[1], x[1] * y[0], x[1] * y[1]])
    return dn(p.min(axis=0)), up(p.max(axis=0))


def scal(c, x):
    """exact float constant c times interval x"""
    a, b = c * x[0], c * x[1]
    return dn(np.minimum(a, b)), up(np.maximum(a, b))


def sq(x):
    a, b = x[0] * x[0], x[1] * x[1]
    lo = np.where((x[0] <= 0) & (x[1] >= 0), 0.0, np.minimum(a, b))
    return np.maximum(dn(lo), 0.0) if np.ndim(lo) else max(dn(lo), 0.0), up(np.maximum(a, b))


def recip(y):
    if np.any((y[0] <= 0) & (y[1] >= 0)):
        raise ZeroDivisionError('interval contains 0')
    return dn(1.0 / y[1]), up(1.0 / y[0])


def div(x, y):
    bad = (y[0] <= 0) & (y[1] >= 0)
    lo1 = np.where(bad, 1.0, y[0])
    hi1 = np.where(bad, 1.0, y[1])
    r = mul(x, (dn(1.0 / hi1), up(1.0 / lo1)))
    return np.where(bad, -INF, r[0]), np.where(bad, INF, r[1])


def sqrt(x):
    lo = np.maximum(x[0], 0.0)
    return np.maximum(dn(np.sqrt(lo)), 0.0), up(np.sqrt(np.maximum(x[1], 0.0)))


def log(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        lo = np.where(x[0] > 0, dn4(np.log(np.maximum(x[0], 1e-300))), -INF)
        return lo, up4(np.log(x[1]))


def exp(x):
    return np.maximum(dn4(np.exp(x[0])), 0.0), up4(np.exp(x[1]))


def sinh(x):
    return dn4(np.sinh(x[0])), up4(np.sinh(x[1]))


def cosh(x):
    a, b = np.cosh(x[0]), np.cosh(x[1])
    lo = np.where((x[0] <= 0) & (x[1] >= 0), 1.0, np.minimum(a, b))
    return np.maximum(dn4(lo), 1.0), up4(np.maximum(a, b))


def arcsinh(x):
    return dn4(np.arcsinh(x[0])), up4(np.arcsinh(x[1]))


def arccos(x):
    """decreasing on [-1, 1]"""
    lo = np.clip(x[1], -1.0, 1.0)
    hi = np.clip(x[0], -1.0, 1.0)
    return np.maximum(dn4(np.arccos(lo)), 0.0), up4(np.arccos(hi))


def arcsin(x):
    lo = np.clip(x[0], -1.0, 1.0)
    hi = np.clip(x[1], -1.0, 1.0)
    return dn4(np.arcsin(lo)), up4(np.arcsin(hi))


def _contains(x, c, period):
    """does [x0, x1] contain c + k*period for some integer k?"""
    k = np.ceil((x[0] - c) / period)
    return c + k * period <= x[1]


def sin(x):
    a, b = np.sin(x[0]), np.sin(x[1])
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    wide = (x[1] - x[0]) >= 2 * np.pi
    hi = np.where(_contains(x, np.pi / 2, 2 * np.pi) | wide, 1.0, hi)
    lo = np.where(_contains(x, -np.pi / 2, 2 * np.pi) | wide, -1.0, lo)
    # the tests use fl(pi/2) and fl(2 pi); widen the flags by a tiny guard band
    g = 1e-12
    hi = np.where(_contains((x[0] - g, x[1] + g), np.pi / 2, 2 * np.pi), 1.0, hi)
    lo = np.where(_contains((x[0] - g, x[1] + g), -np.pi / 2, 2 * np.pi), -1.0, lo)
    return np.maximum(dn4(lo), -1.0), np.minimum(up4(hi), 1.0)


def cos(x):
    return sin(add(x, (dn(np.pi / 2), up(np.pi / 2))))


def arctan(x):
    return dn4(np.arctan(x[0])), up4(np.arctan(x[1]))


def _absrange(x):
    """range of |t| for t in [x0, x1]"""
    a0, a1 = np.abs(x[0]), np.abs(x[1])
    lo = np.where((x[0] <= 0) & (x[1] >= 0), 0.0, np.minimum(a0, a1))
    return lo, np.maximum(a0, a1)


def _ratio(f, y):
    """f(y)/y at points y >= 0, with value 1 at y = 0 (f(0) = 0, f'(0) = 1); NaN stays NaN (AUDIT_3, N6)"""
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(y > 0, f(y) / np.where(y > 0, y, 1.0), np.where(np.isnan(y), np.nan, 1.0))


def sinc(x):
    """sin(t)/t: even, decreasing in |t| on [0, pi]; requires |t| < pi"""
    lo_a, hi_a = _absrange(x)
    if np.any(hi_a >= np.pi):
        raise ValueError('sinc argument too large')
    return dn4(dn(_ratio(np.sin, hi_a))), np.minimum(up4(up(_ratio(np.sin, lo_a))), 1.0)


def asinc(x):
    """arcsin(t)/t: even, increasing in |t| on [0, 1)"""
    lo_a, hi_a = _absrange(x)
    if np.any(hi_a >= 1.0):
        raise ValueError('asinc argument too large')
    return np.maximum(dn4(dn(_ratio(np.arcsin, lo_a))), 1.0), up4(up(_ratio(np.arcsin, hi_a)))


def asinhc(x):
    """arcsinh(t)/t: even, decreasing in |t|"""
    lo_a, hi_a = _absrange(x)
    return dn4(dn(_ratio(np.arcsinh, hi_a))), np.minimum(up4(up(_ratio(np.arcsinh, lo_a))), 1.0)


def hull(*xs):
    return np.minimum.reduce([x[0] for x in xs]), np.maximum.reduce([x[1] for x in xs])


def mid(x):
    return 0.5 * (x[0] + x[1])
