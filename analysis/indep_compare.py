#!/usr/bin/env python3
"""indep_compare.py INDEP.json PROD_AA.json PROD_PP.json OUT.md : compare the independent re-analysis
(ana/indep_check.py) with analyze_arm.py on the same events.  Feed-down kinematics use different random numbers,
so agreement is statistical: |diff| is compared with the production jackknife error (pull = diff / err)."""
import json
import sys

import numpy as np

I, A, P = (json.load(open(f)) for f in sys.argv[1:4])
OUT = sys.argv[4]


def pv(r, k):
    return np.array(r[k]['val'], float), np.array(r[k]['err'], float)


rows = []


def cmp(name, a, b, e):
    a, b, e = np.atleast_1d(a), np.atleast_1d(b), np.atleast_1d(e)
    ok = np.isfinite(a) & np.isfinite(b) & (e > 0)
    pull = (a - b)[ok] / e[ok]
    rows.append((name, int(ok.sum()), float(np.max(np.abs(a - b)[ok])) if ok.any() else np.nan,
                 float(np.max(np.abs(pull))) if ok.any() else np.nan, float(np.sqrt(np.mean(pull ** 2))) if ok.any() else np.nan))


for k in ('v1_D0avg_10-80', 'v1_D0_10-80', 'v1_D0bar_10-80'):
    cmp(k, pv(I, k)[0], *pv(A, k))
for k in ('slope_D0avg_10-80', 'slope_dD0_10-80'):
    cmp(k, I[k][0], A[k][0], A[k][1])
for cl in ('10-40', '0-80'):
    cmp(f'v2{{RP}} D0 {cl}', pv(I, f'v2_D0_{cl}')[0], *pv(A, f'v2_D0_{cl}'))
    cmp(f'v2{{SP}} D0 {cl}', pv(I, f'v2SPS_D0_{cl}')[0], *pv(A, f'v2SPS_D0_{cl}'))
cmp('p+p D0 per c', pv(I, 'pp_dN_D0_perc')[0], *pv(P, 'dN_D0_mb_perc'))
for cl in ('0-10', '10-40', '40-80'):
    cmp(f'Au+Au D0 per c {cl}', pv(I, f'dN_D0_{cl}_perc')[0], *pv(A, f'dN_D0_{cl}_perc'))
    a, ae = pv(A, f'dN_D0_{cl}_perc'); p, pe = pv(P, 'dN_D0_mb_perc')
    with np.errstate(invalid='ignore', divide='ignore'):
        cmp(f'R_AA {cl}', pv(I, f'RAA_{cl}')[0], a / p, a / p * np.hypot(ae / a, pe / p))
a, ae = pv(A, 'dN_Lc_lcbins_10-80_perc'); b, be = pv(A, 'dN_D0_lcbins_10-80_perc')
cmp('Λc/D0 vs pT 10-80', pv(I, 'LcD0_pt_10-80')[0], a / b, a / b * np.hypot(ae / a, be / b))
for cl in ('50-80', '20-50', '0-20'):
    a, ae = pv(A, f'dN_Lc_3to6_{cl}_perc'); b, be = pv(A, f'dN_D0_3to6_{cl}_perc')
    cmp(f'Λc/D0 3-6 GeV {cl}', I[f'LcD0_3to6_{cl}'][0], a[0] / b[0], a[0] / b[0] * np.hypot(ae[0] / a[0], be[0] / b[0]))
L = ['# Independent re-analysis vs production analysis, same events', '',
     f'Events: {I["events"]}. Production counts: Au+Au {A["counts"]["events"]}, p+p {P["counts"]["events"]}.', '',
     '| observable | points | max |Δ| | max |pull| | rms pull |', '|---|---|---|---|---|']
L += [f'| {n} | {k} | {d:.2e} | {mp:.2f} | {rp:.2f} |' for n, k, d, mp, rp in rows]
open(OUT, 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
