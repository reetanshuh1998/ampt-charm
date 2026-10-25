#!/usr/bin/env python3
"""eps09s_errw.py ARMDIR [TOOL] : per-event EPS09s error-set weights (README Sec. 2.3).

For every hq_*.log.gz in ARMDIR ('# hqshw' lines, written by runs with AMPT_HQPDF=1 and AMPT_HQSHW=1, i.e. the
reference configuration) writes ARMDIR/hqerr_<name>.txt.gz: one line per event 'iev w_1 ... w_31'.  Gluon legs
(kf 21) are re-evaluated with TOOL (default bin/eps09s_rg, built by 'make eps-tool'; sets 1..31); quark legs keep the
logged central factor (their error-set variation is not re-evaluated: it needs the free CTEQ6.1M densities; their
share, about 4 %, is printed).  Events HQSHWT left unmodified (logged r1 = r2 = 1 exactly: charm-less) get w = 1 in
every set, as in the in-code weight.  Closure: set 1 must reproduce the logged gluon factors (asserted to 5e-4, the
logging precision).  Then: analysis/analyze_arm.py ARMDIR OUT__epsKK.json --shw tab:K for K = 1..31, and
analysis/eps_hessian.py."""
import glob
import gzip
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = sys.argv[1]
TOOL = os.path.abspath(sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, 'bin', 'eps09s_rg'))
TABLES = os.path.join(ROOT, 'deps', 'eps09s')          # eps09sNLO1..31.dat, read from the working directory
assert os.path.isfile(os.path.join(TABLES, 'eps09sNLO31.dat')), f'{TABLES}/eps09sNLO*.dat missing: run make deps'
assert os.access(TOOL, os.X_OK), f'{TOOL} missing: run make eps-tool'
legs, ev = [], []                    # legs: (x, Q2, s) of gluon legs; ev: (file, iev, [leg index or logged r] x 2)
nq = ng = 0
for hq in sorted(glob.glob(os.path.join(D, 'hq_*.log.gz'))):
    for l in gzip.open(hq, 'rt'):
        if '# hqshw' not in l:
            continue
        f = l.split()
        iev, k1, k2 = int(f[3]), f[5], f[6]
        x, q2, s, r = (float(f[8]), float(f[9])), float(f[11]), (float(f[13]), float(f[14])), (float(f[16]), float(f[17]))
        if r == (1.0, 1.0):          # HQSHWT left the event unmodified (charm-less: x repeated, kf = 0, or IHNT2(16) = 1)
            ev.append((hq, iev, None)); continue
        slot = []
        for b, kf in enumerate((k1, k2)):
            if kf == '21':
                slot.append(('g', len(legs), r[b])); legs.append((x[b], q2, s[b])); ng += 1
            else:
                slot.append(('q', r[b], r[b])); nq += 1
        ev.append((hq, iev, slot))
assert ev, f"no '# hqshw' lines in {D}/hq_*.log.gz (needs AMPT_HQPDF=1, AMPT_HQSHW=1)"
print(f'# {len(ev)} events, gluon legs {ng}, quark legs {nq} ({nq / max(ng + nq, 1):.3f} of legs keep the central factor)')
inp = ''.join(f'{a:.6e} {b:.6f} {c:.6f}\n' for a, b, c in legs)
out = subprocess.run([TOOL, '2.0'], input=inp, capture_output=True, text=True, cwd=TABLES, check=True).stdout
RG = np.array(out.split(), float).reshape(-1, 31)
assert RG.shape == (len(legs), 31), RG.shape
cur, buf, nchk, dmax, wsum = None, [], 0, 0.0, np.zeros(31)


def flush():
    if cur is not None:
        name = os.path.basename(cur).replace('hq_', 'hqerr_').replace('.log.gz', '.txt.gz')
        with gzip.open(os.path.join(D, name), 'wt') as g:
            g.writelines(buf)


for hq, iev, slot in ev:
    if hq != cur:
        flush(); cur, buf = hq, []
    w = np.ones(31)
    for t, a, rl in slot or ():
        if t == 'g':
            d = abs(RG[a, 0] - rl); dmax = max(dmax, d); nchk += 1
            assert d < 5e-4, (hq, iev, RG[a, 0], rl)
            w *= RG[a]
        else:
            w *= a
    wsum += w
    buf.append(f'{iev} ' + ' '.join(f'{v:.6f}' for v in w) + '\n')
flush()
print(f'# closure: set 1 vs logged gluon factors, {nchk} legs, max |diff| {dmax:.2e}')
print('# mean event weight per set:', ' '.join(f'{k + 1}:{v:.4f}' for k, v in enumerate(wsum / len(ev))))
