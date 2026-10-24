"""Coverage check for the alpha sweep: do the certified runs of mixed_certF.py cover every F in [F(1/4), F(2)] =
[13/12, 13/5]?

Each certified slice [lo, hi] of a sweep (mixed_alpha_sweep.py) certified Claim M for the F-box
(ia.dn(F(lo)), ia.up(F(hi))) computed in floating point exactly as mixed_certF.py does; that box is recomputed here. The
union of the boxes must contain [13/12, 13/5]; its two ends are compared with 13/12 and 13/5 as exact rationals. This
needs no error analysis of F(alpha) in floating point. F(alpha) = (b+1)^2/(2b) - 1, b = 2 alpha + 1, is continuous and
increasing, so F in [13/12, 13/5] is the same as alpha in [1/4, 2].

Sweep outputs (stdout of mixed_alpha_sweep.py) are given as  file:a_lo:a_hi:w . Older outputs print the endpoints
with 6 decimals; their exact floating-point values are recovered by replaying the slicing of the sweep (initial
slices a, a + w, ... and halvings down to w/8). Newer outputs print them exactly.
usage: python sweep_coverage.py out1:a_lo:a_hi:w [out2:...] ...
       python sweep_coverage.py runs.txt       (one out:a_lo:a_hi:w per line, paths relative to runs.txt)"""
import os
import re
import sys
from fractions import Fraction

import ia
from mixed_certF import _F_of_alpha


def replay_nodes(a_lo, a_hi, w):
    nodes, a = [], a_lo
    while a < a_hi - 1e-12:
        nodes.append((a, min(a + w, a_hi)))
        a += w
    out, layer = list(nodes), list(nodes)
    for _ in range(3):
        layer = [h for lo, hi in layer for h in ((lo, 0.5 * (lo + hi)), (0.5 * (lo + hi), hi))]
        out += layer
    return out


def certified_slices(path, a_lo, a_hi, w):
    nodes = replay_nodes(a_lo, a_hi, w)
    found = []
    for line in open(path, encoding='utf-8'):
        m = re.match(r'alpha \[([^,]+), ([^\]]+)\]: CERTIFIED', line.strip())
        if not m:
            continue
        x, y = m.group(1), m.group(2)
        if len(x.split('.')[-1]) == 6 and len(y.split('.')[-1]) == 6:      # 6-decimal format: replay
            match = [(lo, hi) for lo, hi in nodes if '%.6f' % lo == x and '%.6f' % hi == y]
            if len(match) != 1:
                raise SystemExit('cannot identify slice %s in %s (%d matches)' % (line.strip(), path, len(match)))
            found.append(match[0])
        else:
            found.append((float(x), float(y)))
    return found


if __name__ == '__main__':
    specs = sys.argv[1:]
    if len(specs) == 1 and specs[0].endswith('.txt'):
        base = os.path.dirname(os.path.abspath(specs[0]))
        specs = [os.path.join(base, l.strip()) for l in open(specs[0], encoding='utf-8')
                 if l.strip() and not l.lstrip().startswith('#')]
    boxes = []
    for spec in specs:
        path, a_lo, a_hi, w = spec.rsplit(':', 3)
        sl = certified_slices(path, float(a_lo), float(a_hi), float(w))
        print('%-40s %3d certified slices' % (os.path.basename(path), len(sl)))
        boxes += [(ia.dn(_F_of_alpha(lo)), ia.up(_F_of_alpha(hi)), lo, hi) for lo, hi in sl]
    boxes.sort()
    target_lo, target_hi = Fraction(13, 12), Fraction(13, 5)
    if Fraction(boxes[0][0]) > target_lo:
        raise SystemExit('lowest F-box starts at %r > 13/12' % boxes[0][0])
    reach, gaps = boxes[0][1], []
    for flo, fhi, lo, hi in boxes[1:]:
        if flo > reach:
            gaps.append((reach, flo))
        reach = max(reach, fhi)
    covered = not gaps and Fraction(reach) >= target_hi
    print('union of certified F-boxes: [%r, %r]; gaps: %s'
          % (float(boxes[0][0]), float(reach), [(float(a), float(b)) for a, b in gaps] if gaps else 'none'))
    print('alpha-slices span [%r, %r]' % (min(b[2] for b in boxes), max(b[3] for b in boxes)))
    print('COVERS F in [13/12, 13/5] (alpha in [1/4, 2]): %s' % ('YES' if covered else 'NO'))
