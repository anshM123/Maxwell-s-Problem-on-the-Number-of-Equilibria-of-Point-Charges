"""AUDIT_3, item 3: independent check of the alpha-sweep bookkeeping (cert/runs.txt, cert/alpha_sweep*.out/.log).

Own implementation (does not import code/sweep_coverage.py):
 1. every output named in cert/runs.txt is parsed; its CERTIFIED slices are identified
    - new format (repr): float(text) is exact;
    - old format (6 decimals): by an own replay of mixed_alpha_sweep.py's slicing (a += w, halving 3 times), matched on
      the 6-decimal text; as a robustness check the coverage is ALSO computed with the decimal values read literally;
 2. every counted slice must have, in the matching .log, its header line followed by four chart lines 'FAILED 0' and
    a 'TOTAL failed: 0 -> CLAIM M CERTIFIED for alpha in [%g, %g]' line whose %g values match the slice;
 3. F-box of a slice = (dn(F(lo)), up(F(hi))), F(a) = (b+1)^2/(2b) - 1 with b = 2a + 1 in binary64, dn/up = two
    nextafter steps (re-implemented here, and compared bit-for-bit with mixed_certF._F_of_alpha / ia.dn / ia.up);
 4. the union of the F-boxes is compared with [13/12, 13/5] in exact rational arithmetic; the smallest overlap between
    consecutive boxes is reported; the sliver run's single F-box is checked the same way;
 5. completeness: does every output end with a 'SWEEP ... ALL CERTIFIED' line, and is every 'not closed' slice
    replaced by certified sub-slices (in alpha) somewhere in the certificate?
 6. (wording) does each F-box contain the exact F-range of its slice, as the paper says?
usage: python c2_coverage.py   (run from audit/)"""
import math
import os
import re
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, '..', 'cert')
sys.path.insert(0, os.path.join(HERE, '..', 'code'))


def nx_dn(x):
    return math.nextafter(math.nextafter(x, -math.inf), -math.inf)


def nx_up(x):
    return math.nextafter(math.nextafter(x, math.inf), math.inf)


def F_float(a):
    b = 2 * a + 1
    return (b + 1) ** 2 / (2 * b) - 1


def F_exact(a):
    b = 2 * Fraction(a) + 1
    return (b + 1) ** 2 / (2 * b) - 1


def replay(a_lo, a_hi, w):
    nodes, a = [], a_lo
    while a < a_hi - 1e-12:
        nodes.append((a, min(a + w, a_hi), 0))
        a += w
    out, layer = list(nodes), list(nodes)
    for lev in range(1, 4):
        nl = []
        for lo, hi, _ in layer:
            m = 0.5 * (lo + hi)
            nl += [(lo, m, lev), (m, hi, lev)]
        layer = nl
        out += layer
    return out


def parse_out(path, a_lo, a_hi, w):
    nodes = replay(a_lo, a_hi, w)
    cert, notclosed, sweep_line, literal = [], [], None, []
    for line in open(path, encoding='utf-8'):
        line = line.strip()
        if line.startswith('SWEEP'):
            sweep_line = line
            continue
        m = re.match(r'alpha \[([^,]+), ([^\]]+)\]: (CERTIFIED|not closed)', line)
        if not m:
            continue
        x, y, st = m.groups()
        old = len(x.split('.')[-1]) == 6 and len(y.split('.')[-1]) == 6
        if old:
            mt = [(lo, hi) for lo, hi, _ in nodes if '%.6f' % lo == x and '%.6f' % hi == y]
            if len(mt) != 1:
                raise SystemExit('ambiguous/unmatched %s in %s: %d' % (line, path, len(mt)))
            lo, hi = mt[0]
        else:
            lo, hi = float(x), float(y)
            if repr(lo) != x or repr(hi) != y:
                raise SystemExit('non-round-trip text %s' % line)
        (cert if st == 'CERTIFIED' else notclosed).append((lo, hi, line))
        if st == 'CERTIFIED':
            literal.append((float(x), float(y)))
    return cert, notclosed, sweep_line, literal


def check_log(logpath, header, lo, hi):
    lines = open(logpath, encoding='utf-8').read().splitlines()
    idx = [i for i, l in enumerate(lines) if l.strip() == header]
    if not idx:
        return 'header missing'
    i = idx[-1]
    body = [l.strip() for l in lines[i + 1:i + 6]]
    names = ['main', 'band', 'corner +1', 'corner -1']
    for nm, l in zip(names, body[:4]):
        m = re.match(r'(main|band|corner [+-]1)\s+certified\s+(\d+), FAILED (\d+)', l)
        if not m or m.group(1) != nm or m.group(3) != '0':
            return 'bad chart line %r' % l
    tot = body[4]
    want = 'TOTAL failed: 0  ->  CLAIM M CERTIFIED for alpha in [%g, %g]' % (lo, hi)
    if tot != want:
        return 'bad TOTAL line %r (want %r)' % (tot, want)
    return None


def union_check(boxes, label):
    boxes = sorted(boxes)
    lo0 = Fraction(boxes[0][0])
    reach = boxes[0][1]
    gaps, min_overlap = [], None
    for flo, fhi, *_ in boxes[1:]:
        if flo > reach:
            gaps.append((reach, flo))
        else:
            ov = Fraction(reach) - Fraction(flo)
            min_overlap = ov if min_overlap is None else min(min_overlap, ov)
        reach = max(reach, fhi)
    ok = lo0 <= Fraction(13, 12) and Fraction(reach) >= Fraction(13, 5) and not gaps
    print('%s: %d boxes, union [%r, %r], gaps %s, smallest overlap %.3e, covers [13/12, 13/5]: %s'
          % (label, len(boxes), boxes[0][0], reach, gaps or 'none', float(min_overlap or 0), 'YES' if ok else 'NO'))
    return ok


def main():
    import ia
    import mixed_certF
    specs = [l.strip() for l in open(os.path.join(CERT, 'runs.txt'), encoding='utf-8')
             if l.strip() and not l.lstrip().startswith('#')]
    all_boxes, lit_boxes, all_alpha, all_nc = [], [], [], []
    problems = 0
    for spec in specs:
        fn, a_lo, a_hi, w = spec.rsplit(':', 3)
        path = os.path.join(CERT, fn)
        logpath = path[:-4] + '.log'
        cert, nc, sweep_line, literal = parse_out(path, float(a_lo), float(a_hi), float(w))
        bad = 0
        for lo, hi, line in cert:
            header = line
            err = check_log(logpath, header, lo, hi)
            if err:
                bad += 1
                print('   LOG PROBLEM %s: %s -> %s' % (fn, line, err))
        problems += bad
        print('%-22s %3d certified, %2d not closed, log problems %d, final line: %s'
              % (fn, len(cert), len(nc), bad, sweep_line or '(none: run stopped)'))
        for lo, hi, _ in cert:
            flo, fhi = nx_dn(F_float(lo)), nx_up(F_float(hi))
            # bit-for-bit comparison with the code under test
            assert flo == ia.dn(mixed_certF._F_of_alpha(lo)) and fhi == ia.up(mixed_certF._F_of_alpha(hi))
            all_boxes.append((flo, fhi, lo, hi, fn))
            all_alpha.append((lo, hi))
        for x, y in literal:
            lit_boxes.append((nx_dn(F_float(x)), nx_up(F_float(y))))
        all_nc += [(lo, hi, fn) for lo, hi, _ in nc]
    print()
    ok1 = union_check(all_boxes, 'certified F-boxes (replayed floats)')
    ok2 = union_check(lit_boxes, 'certified F-boxes (decimal text read literally)')
    # alpha-space view (only a sanity check; the F-union above is what matters)
    all_alpha.sort()
    reach, agaps = all_alpha[0][1], []
    for lo, hi in all_alpha[1:]:
        if lo > reach:
            agaps.append((reach, lo))
        reach = max(reach, hi)
    print('alpha-slices: %d, span [%r, %r], alpha-gaps: %s' % (len(all_alpha), all_alpha[0][0], reach, agaps or 'none'))
    # every not-closed slice must be covered by certified slices
    unc = 0
    for lo, hi, fn in all_nc:
        inside = sorted((a, b) for a, b in all_alpha if b > lo and a < hi)
        r = lo
        for a, b in inside:
            if a <= r:
                r = max(r, b)
        if r < hi:
            unc += 1
            print('   not-closed slice [%r, %r] (%s) is not re-covered (reach %r)' % (lo, hi, fn, r))
    print('not-closed slices re-covered by certified slices: %d of %d' % (len(all_nc) - unc, len(all_nc)))
    # the sliver run's F-box
    sl = (nx_dn(F_float(0.25)), nx_up(F_float(2.0)))
    print('sliver F-box [%r, %r] contains [13/12, 13/5]: %s'
          % (sl[0], sl[1], Fraction(sl[0]) <= Fraction(13, 12) and Fraction(sl[1]) >= Fraction(13, 5)))
    # exactness facts used: F increasing, F(1/4) = 13/12, F(2) = 13/5
    print('F(1/4) = %s, F(2) = %s (exact); F float at 1/4: %r (%s 13/12), at 2: %r (%s 13/5)'
          % (F_exact(Fraction(1, 4)), F_exact(2), F_float(0.25),
             '>' if Fraction(F_float(0.25)) > Fraction(13, 12) else '<=', F_float(2.0),
             '<' if Fraction(F_float(2.0)) < Fraction(13, 5) else '>='))
    # wording check for the paper ('for all F in a box containing the values of F on the slice'): does each float
    # F-box contain the exact F-range F([lo, hi]) of its slice (lo, hi taken as exact binary64 values)?
    worst_lo = worst_hi = None
    bad = 0
    for flo, fhi, lo, hi, fn in all_boxes:
        dl = F_exact(Fraction(lo)) - Fraction(flo)
        dh = Fraction(fhi) - F_exact(Fraction(hi))
        bad += (dl < 0) or (dh < 0)
        worst_lo = dl if worst_lo is None else min(worst_lo, dl)
        worst_hi = dh if worst_hi is None else min(worst_hi, dh)
    print('each F-box contains the exact F-range of its slice: %s (%d exceptions); smallest margins %.2e (low end), '
          '%.2e (high end)' % (bad == 0, bad, float(worst_lo), float(worst_hi)))
    print('VERDICT: log problems %d; F-union covers: %s (replayed), %s (literal)' % (problems, ok1, ok2))


if __name__ == '__main__':
    main()
