"""Re-runs the band chart of every certified alpha-slice listed in a runs file (cert/runs.txt) with the current ia.py,
after the fix that makes ia._ratio propagate NaN (AUDIT_3, note N6), and compares the certified box count with the count
in the run's log. The band chart is the only caller of ia.sinc / asinc / asinhc, so the main and corner charts of the
logged runs do not depend on the fix. Each band chart is run exactly as in mixed_certF.py (e0 = 0.4, the same grid and
F-box).
usage: python band_recheck.py ../cert/runs.txt"""
import contextlib
import io
import os
import re
import sys

import ia
from mixed_certF import I4_band, _F_of_alpha, bnb
from sweep_coverage import certified_slices

E0 = 0.4


def logged_band_counts(log_path):
    """band counts of the CERTIFIED slices of a log, in file order"""
    counts, status = [], None
    for line in open(log_path, encoding='utf-8'):
        if line.startswith('alpha ['):
            status = 'CERTIFIED' if ': CERTIFIED' in line else 'other'
        m = re.match(r'\s+band\s+certified\s+(\d+), FAILED (\d+)', line)
        if m and status == 'CERTIFIED':
            counts.append((int(m.group(1)), int(m.group(2))))
    return counts


if __name__ == '__main__':
    base = os.path.dirname(os.path.abspath(sys.argv[1]))
    specs = [l.strip() for l in open(sys.argv[1], encoding='utf-8') if l.strip() and not l.lstrip().startswith('#')]
    total, mismatches, failures = 0, 0, 0
    for spec in specs:
        out, a_lo, a_hi, w = spec.rsplit(':', 3)
        sl = certified_slices(os.path.join(base, out), float(a_lo), float(a_hi), float(w))
        logged = logged_band_counts(os.path.join(base, out[:-4] + '.log'))
        if len(logged) != len(sl):
            raise SystemExit('%s: %d certified slices but %d logged band counts' % (out, len(sl), len(logged)))
        for (lo, hi), (n_log, f_log) in zip(sl, logged):
            F_dom = (ia.dn(_F_of_alpha(lo)), ia.up(_F_of_alpha(hi)))
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rest = bnb('band', I4_band, [(0.0, E0), (0.0, 1.0), (-1.0, 1.0), F_dom], (8, 16, 16, 2), max_rounds=60)
            n_now = int(re.search(r'certified\s+(\d+)', buf.getvalue()).group(1))
            total += 1
            failures += rest[0].shape[1]
            if n_now != n_log or rest[0].shape[1] != 0 or f_log != 0:
                mismatches += 1
                print('MISMATCH %s [%r, %r]: logged %d (failed %d), now %d (failed %d)'
                      % (out, lo, hi, n_log, f_log, n_now, rest[0].shape[1]), flush=True)
        print('%-22s %3d slices re-run' % (out, len(sl)), flush=True)
    print('band charts re-run: %d; count mismatches or failures: %d; failed boxes: %d  ->  %s'
          % (total, mismatches, failures, 'ALL IDENTICAL TO THE LOGS' if mismatches == 0 else 'CHECK'))
