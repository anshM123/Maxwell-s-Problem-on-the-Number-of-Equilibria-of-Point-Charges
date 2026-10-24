"""AUDIT_3, item 8: totals behind the paper's sentence 'about 1.3 x 10^9 boxes ... about four hours on one core'.
Sums, over the 163 counted slices (the CERTIFIED lines of the outputs in cert/runs.txt), the per-chart box counts and
per-chart times of the matching .log entries; also over all logged slices (including not-closed attempts).
Also reports the smallest and largest slice widths.
Output: audit/c9_log_totals.out"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, '..', 'cert')
sys.path.insert(0, HERE)
from c2_coverage import parse_out  # noqa: E402

tot_boxes = tot_time = 0
all_boxes = all_time = 0
widths = []
for spec in [l.strip() for l in open(os.path.join(CERT, 'runs.txt'), encoding='utf-8') if l.strip() and l[0] != '#']:
    fn, a_lo, a_hi, w = spec.rsplit(':', 3)
    cert, nc, _, _ = parse_out(os.path.join(CERT, fn), float(a_lo), float(a_hi), float(w))
    lines = open(os.path.join(CERT, fn[:-4] + '.log'), encoding='utf-8').read().splitlines()
    heads = {l.strip() for _, _, l in cert}
    cur = None
    for l in lines:
        if l.startswith('alpha ['):
            cur = l.strip()
            continue
        m = re.match(r'\s*(main|band|corner [+-]1)\s+certified\s+(\d+), FAILED (\d+)\s+\((\d+) rounds, (\d+)s\)', l)
        if m and cur is not None:
            n = int(m.group(2)) + int(m.group(3))
            t = int(m.group(5))
            all_boxes += n
            all_time += t
            if cur in heads:
                tot_boxes += n
                tot_time += t
    widths += [hi - lo for lo, hi, _ in cert]
print('counted slices: %d; boxes certified in them: %.3e; chart time (as logged, rounded per chart): %.2f h'
      % (len(widths), tot_boxes, tot_time / 3600))
print('all logged attempts (incl. not closed): boxes evaluated at the end %.3e; chart time %.2f h' % (all_boxes, all_time / 3600))
print('slice widths: min %r, max %r' % (min(widths), max(widths)))
