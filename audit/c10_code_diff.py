"""AUDIT_3, item 2: line-by-line difference between the interval-F chart functions (code/mixed_certF.py) and the audited
Coulomb chart functions (code/mixed_cert.I_iv, code/mixed_charts.I4_band, code/mixed_charts.J1n_corner).
Expected: only the replacement of the Coulomb constants F - 1 = 1/4, F + 1 = 9/4, q = 1/3 by the intervals Fm, Fp and
sqrt(Fm/Fp), plus cosmetic changes.
Output: audit/c10_code_diff.out"""
import difflib
import inspect
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'code'))
import mixed_cert  # noqa: E402
import mixed_certF  # noqa: E402
import mixed_charts  # noqa: E402

pairs = [(mixed_cert.I_iv, mixed_certF.I_main), (mixed_charts.I4_band, mixed_certF.I4_band),
         (mixed_charts.J1n_corner, mixed_certF.J1n_corner)]
for a, b in pairs:
    la = inspect.getsource(a).splitlines()[1:]
    lb = inspect.getsource(b).splitlines()[1:]
    d = [l for l in difflib.unified_diff(la, lb, lineterm='', n=0) if not l.startswith(('---', '+++', '@@'))]
    print('== %s.%s  ->  %s.%s : %d differing lines' % (a.__module__, a.__name__, b.__module__, b.__name__, len(d)))
    for l in d:
        print('   ' + l.strip())
