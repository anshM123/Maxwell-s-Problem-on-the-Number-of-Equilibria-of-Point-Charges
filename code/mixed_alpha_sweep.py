"""Certify Claim M for every alpha in [a_lo, a_hi] by slices of width w (runs mixed_certF.py on each slice).
Each slice is retried with half the width (down to w/8) if it does not close. Results go to ../cert/alpha_sweep.log.
usage: python mixed_alpha_sweep.py a_lo a_hi w [e0] [r0] [log]"""
import subprocess
import sys
import time

a_lo, a_hi, w = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
e0 = sys.argv[4] if len(sys.argv) > 4 else '0.4'
r0 = sys.argv[5] if len(sys.argv) > 5 else '0.05'
py = sys.executable
log = open(sys.argv[6] if len(sys.argv) > 6 else '../cert/alpha_sweep.log', 'a', encoding='utf-8')


def run(lo, hi):
    out = subprocess.run([py, 'mixed_certF.py', repr(lo), repr(hi), e0, r0, '60'], capture_output=True, text=True)
    ok = 'CLAIM M CERTIFIED' in out.stdout
    lines = [l for l in out.stdout.splitlines() if 'certified' in l or 'TOTAL' in l]
    return ok, lines


stack = []
a = a_lo
while a < a_hi - 1e-12:
    stack.append((a, min(a + w, a_hi)))
    a += w
stack.reverse()
t0 = time.time()
failed = []
while stack:
    lo, hi = stack.pop()
    ok, lines = run(lo, hi)
    msg = 'alpha [%r, %r]: %s (%.0fs elapsed)' % (lo, hi, 'CERTIFIED' if ok else 'not closed', time.time() - t0)
    print(msg, flush=True)
    log.write(msg + '\n' + '\n'.join('    ' + l for l in lines) + '\n')
    log.flush()
    if not ok:
        if hi - lo > w / 8 + 1e-12:
            mid = 0.5 * (lo + hi)
            stack.append((mid, hi))
            stack.append((lo, mid))
        else:
            failed.append((lo, hi))
summary = 'SWEEP alpha in [%g, %g]: %s' % (a_lo, a_hi, 'ALL CERTIFIED' if not failed else 'uncertified slices %s' % failed)
print(summary, flush=True)
log.write(summary + '\n')
