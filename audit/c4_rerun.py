"""AUDIT_3, item 3: independent re-runs of certified alpha-slices and comparison with the logged box counts.

For one or two slices of every output in cert/runs.txt (including the slice ending at alpha = 2), code/mixed_certF.py is
run again (as a subprocess, one at a time, from code/, with e0 = 0.4, r0 = 0.05, 60 rounds: the arguments that
mixed_alpha_sweep.py passes), with the slice end points given as exact floats (old-format outputs: the replayed values
of c2_coverage.py). The run must print 'CLAIM M CERTIFIED', and its per-chart box counts must equal the counts in the
matching .log entry. The computation is deterministic, so equal counts also confirm that the logged run used
e0 = 0.4 and r0 = 0.05 (the main-chart counts depend on e0, the corner counts on r0).
The script mixed_certF.py writes no files.
usage: python c4_rerun.py [index ...]   (default: all slices below)"""
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, '..', 'code')
CERT = os.path.join(HERE, '..', 'cert')
sys.path.insert(0, HERE)
from c2_coverage import replay  # noqa: E402

PY = sys.executable


def old(a_lo, a_hi, w, x, y):
    m = [(lo, hi) for lo, hi, _ in replay(a_lo, a_hi, w) if '%.6f' % lo == x and '%.6f' % hi == y]
    assert len(m) == 1
    return m[0]


SLICES = [
    ('alpha_sweep.log', '[0.250000, 0.300000]', old(0.25, 2.0, 0.05, '0.250000', '0.300000')),
    ('alpha_sweep.log', '[1.000000, 1.050000]', old(0.25, 2.0, 0.05, '1.000000', '1.050000')),
    ('alpha_sweep_mid.log', '[1.55, 1.5562500000000001]', (1.55, 1.5562500000000001)),
    ('alpha_sweep_hi.log', '[1.700000, 1.706250]', old(1.7, 2.0, 0.025, '1.700000', '1.706250')),
    ('alpha_sweep_hi2.log', '[1.8343750000000016, 1.8375000000000017]', (1.8343750000000016, 1.8375000000000017)),
    ('alpha_sweep_t1.log', '[1.84375, 1.845703125]', (1.84375, 1.845703125)),
    ('alpha_sweep_t2.log', '[1.947265625, 1.94921875]', (1.947265625, 1.94921875)),
    ('alpha_sweep_t3.log', '[1.998046875, 2.0]', (1.998046875, 2.0)),
]


def logged_counts(logname, header_part):
    lines = open(os.path.join(CERT, logname), encoding='utf-8').read().splitlines()
    idx = [i for i, l in enumerate(lines) if l.startswith('alpha ' + header_part + ': CERTIFIED')]
    if not idx:
        return None
    i = idx[-1]
    out = {}
    for l in lines[i + 1:i + 5]:
        m = re.match(r'\s*(main|band|corner [+-]1)\s+certified\s+(\d+), FAILED (\d+)', l)
        out[m.group(1)] = (int(m.group(2)), int(m.group(3)))
    return out


def main():
    sel = [int(a) for a in sys.argv[1:]] or range(len(SLICES))
    allok = True
    for k in sel:
        logname, hdr, (lo, hi) = SLICES[k]
        want = logged_counts(logname, hdr)
        t0 = time.time()
        p = subprocess.run([PY, 'mixed_certF.py', repr(lo), repr(hi), '0.4', '0.05', '60'], cwd=CODE,
                           capture_output=True, text=True)
        got = {}
        for l in p.stdout.splitlines():
            m = re.match(r'\s*(main|band|corner [+-]1)\s+certified\s+(\d+), FAILED (\d+)', l)
            if m:
                got[m.group(1)] = (int(m.group(2)), int(m.group(3)))
        cert = 'CLAIM M CERTIFIED' in p.stdout
        same = got == want
        allok &= cert and same
        print('slice %-44s (%s): rerun %s; counts %s; logged %s; identical: %s  (%.0fs)'
              % (hdr, logname, 'CERTIFIED' if cert else 'NOT CERTIFIED', got, want, same, time.time() - t0), flush=True)
        if p.returncode:
            print(p.stderr[-2000:])
    print('ALL RE-RUNS CERTIFIED WITH LOGGED COUNTS: %s' % allok)


if __name__ == '__main__':
    main()
