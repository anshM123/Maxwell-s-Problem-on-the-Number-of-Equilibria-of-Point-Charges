"""AUDIT item 7: trend of the convexity margin of theta(M) toward degenerate triangles (Coulomb unless alpha given).
Families: thin isosceles (A, A, 180-2A); two small unequal angles (A, 3A, ...); one small angle (A, 60, 120-A) and
(A, (180-A)/2, (180-A)/2).  Uses analyse() of a67_shape_scan.py.
usage: python a67_thin_trend.py [alpha]"""
import math
import sys

alpha = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
sys.argv = [sys.argv[0], str(alpha)]
import a67_shape_scan as m  # noqa: E402

for fam, shapes in (('isosceles (A,A,180-2A)', [(a, a) for a in (0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40, 60)]),
                    ('two small (A,3A)', [(a, 3 * a) for a in (0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40)]),
                    ('one small (A,60)', [(a, 60) for a in (0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40, 60)]),
                    ('one small (A,(180-A)/2)', [(a, (180 - a) / 2) for a in (0.05, 0.1, 0.5, 2, 10, 30, 60)])):
    print('family', fam, flush=True)
    for A1, A2 in shapes:
        r = m.analyse(math.radians(A1), math.radians(A2))
        print('  A=(%7.3f,%7.3f,%8.3f) kw_min %.5f  kmin %.2e  turn %s  wind [%s,%s]  nzeros %s  nonmono %s  reentry %s  '
              'aspect %.1f' % (A1, A2, 180 - A1 - A2, r['kw_min'], r['kmin'], r['turn'], r['wind_min'], r['wind_max'],
                               r['nzeros'], r['nonmono'], r['reentry'], r['kw_max'] / max(r['kw_min'], 1e-300)), flush=True)
