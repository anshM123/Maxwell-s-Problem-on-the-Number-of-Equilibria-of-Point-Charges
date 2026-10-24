"""sanity test of the rigorous mp.iv helpers sinhI, coshI, acosI used in b4e_mpiv_recert.py (enclosure of dense
50-digit samples, including intervals containing 0 and touching +-1)"""
import mpmath as mp
import numpy as np
src = open('b4e_mpiv_recert.py', encoding='utf-8').read()
ns = {'__file__': 'b4e_mpiv_recert.py', '__name__': 'helpers'}
exec(src.split('dat = np.load')[0].replace('prec = int(sys.argv[1]) if len(sys.argv) > 1 else 64', 'prec = 64')
     .replace('stride = int(sys.argv[2]) if len(sys.argv) > 2 else 1', 'stride = 1'), ns)
iv = ns['iv']
rng = np.random.default_rng(1)
bad = 0
for k in range(3000):
    a = rng.uniform(-3, 3); w = 10 ** rng.uniform(-8, 0)
    X = iv.mpf([a, a + w])
    for f, g in ((ns['sinhI'], mp.sinh), (ns['coshI'], mp.cosh)):
        r = f(X)
        for t in np.linspace(a, a + w, 9):
            bad += not (r.a <= g(mp.mpf(float(t))) <= r.b)
    a = rng.uniform(-1, 1); w = min(10 ** rng.uniform(-8, 0), 1 - a)
    r = ns['acosI'](iv.mpf([a, a + w]))
    for t in np.linspace(a, a + w, 9):
        bad += not (r.a <= mp.acos(mp.mpf(float(t))) <= r.b)
r = ns['coshI'](iv.mpf([-0.5, 0.25])); bad += not (r.a <= 1 <= r.b)
r = ns['acosI'](iv.mpf([0.9, 1.0])); bad += not (r.a <= 0)
print('helper enclosure failures: %d (of 3000 x 27 samples + edge cases)' % bad)
