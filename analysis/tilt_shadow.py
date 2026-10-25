#!/usr/bin/env python3
"""tilt_shadow.py OUT.md LABEL=ARM.json ... : shadowing geometry vs charm tilt vs D0 directed flow (STEP_LOG 95-96).

Per arm: rapidity-odd charm/u-quark initial dipole eps1 (formation points, 0.5 < |y| < 2.5) and their ratio (the
charm-to-bulk tilt), c- and u-quark v1 at freeze-out (same window), D0 (D0 + D0bar) dv1/dy and Delta v1 (D0 - D0bar)
in 10-80 % (STAR window, pT > 1.5 GeV), and the mean shadowing weight.  Charm keys carry the event weight, so
re-weighted variants of one sample (analysis --shw) differ only by the shadowing model.  STAR: dv1/dy
-0.080 +- 0.023, Delta -0.011 +- 0.039."""
import json
import sys

import numpy as np

OUT, SPECS = sys.argv[1], sys.argv[2:]


def odd(r, key, lo=0.5, hi=2.5):
    yq = np.array(r['edges']['Y_Q']); c = 0.5 * (yq[1:] + yq[:-1])
    v, e = np.array(r[key]['val'], float), np.array(r[key]['err'], float)
    s = (np.abs(c) > lo) & (np.abs(c) < hi) & np.isfinite(v) & (e > 0)
    w = 1 / e[s] ** 2
    return np.sum(w * np.sign(c[s]) * v[s]) / np.sum(w), 1 / np.sqrt(np.sum(w))


f = lambda t, d=4: f'{t[0]:+.{d}f} ± {t[1]:.{d}f}'
L = ['# Shadowing geometry, charm tilt and D0 directed flow', '',
     'ε1: rapidity-odd ⟨sign(y) cos φ_x⟩ of formation points, 0.5 < |y| < 2.5 (0–80 %). Tilt ratio = ε1(c)/ε1(u). '
     'STAR D0 dv1/dy (10–80 %, pT > 1.5 GeV) −0.080 ± 0.023, Δv1 −0.011 ± 0.039.', '',
     '| arm | ⟨w_shadow⟩ | ε1(u) | ε1(c) | tilt ratio c/u | v1(u) | v1(c) | D0 dv1/dy | Δv1 D0−D̄0 |',
     '|---|---|---|---|---|---|---|---|---|']
for sp in SPECS:
    lab, fn = sp.split('=', 1)
    r = json.load(open(fn))
    c = r['counts']
    w = c.get('sumw_shw', 0) / c['sumw'] if c.get('sumw_shw') else 1.0
    eu, ec = odd(r, 'eps1_uq_all'), odd(r, 'eps1_cq_all')
    rat = (ec[0] / eu[0], abs(ec[0] / eu[0]) * np.hypot(ec[1] / ec[0], eu[1] / eu[0]))
    vu, vc = odd(r, 'v1_uq_all'), odd(r, 'v1_cq_all')
    L.append(f'| {lab} | {w:.3f} | {f(eu)} | {f(ec)} | {f(rat, 3)} | {f(vu)} | {f(vc)} | '
             f'{f(r["slope_D0avg_10-80"])} | {f(r["slope_dD0_10-80"])} |')
open(OUT, 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
