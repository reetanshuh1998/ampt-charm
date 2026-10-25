#!/usr/bin/env python3
"""compare_with_published.py QUICK_AA.json QUICK_PP.json : is a small new production statistically consistent with the
published full-statistics results (observation/results/json/dm25_main.json, 2.4M Au+Au; dm25_ppmain.json, 11.9M p+p)?

Compared: charm-hadron yields per c quark (D0, D+, Ds, Lambda_c; all rapidities) and the D0 mean pT at |y| < 1, in p+p
and in 0-80 % Au+Au, plus the event-level gates (charm conservation, matching).  pull = (new - published) / sqrt(e_new^2
+ e_pub^2).  e_new is the larger of its own jackknife error and the published error scaled to the new sample size,
e_pub sqrt(N_pub / N_new): a jackknife over a handful of job files has too few degrees of freedom to be trusted
(it underestimated a binomial error by a factor 2 for 2 files).  PASS when every |pull| < 3.5 (ten quantities: a
3.5 sigma fluctuation has ~0.5 % probability) and the gates hold."""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
qa, qp = (json.load(open(f)) for f in sys.argv[1:3])
pa = json.load(open(ROOT / 'observation/results/json/dm25_main.json'))
pp = json.load(open(ROOT / 'observation/results/json/dm25_ppmain.json'))
rows, ok = [], True
for sysn, q, p, cl in (('p+p', qp, pp, 'mb'), ('Au+Au 0-80 %', qa, pa, '0-80')):
    scale = math.sqrt(p['counts']['events'] / q['counts']['events'])
    for k, lab in [(f'N_{nm}_{cl}_perc', f'{nm} per c quark') for nm in ('D0', 'Dp', 'Ds', 'Lc')] + \
                  [(f'meanpt_D0_{cl}', 'D0 <pT> |y|<1 (GeV)')]:
        ep = p[k]['err'][0]
        rows.append((sysn, lab, q[k]['val'][0], max(q[k]['err'][0], ep * scale), p[k]['val'][0], ep))
print(f'{"system":14s} {"quantity":22s} {"this run":>18s} {"published":>18s} {"pull":>6s}')
for s, n, a, ea, b, eb in rows:
    pull = (a - b) / math.hypot(ea, eb) if ea > 0 else float('nan')
    ok &= math.isfinite(pull) and abs(pull) < 3.5
    print(f'{s:14s} {n:22s} {a:9.4f} ± {ea:6.4f} {b:9.4f} ± {eb:6.4f} {pull:+6.2f}')
for lab, r in (('Au+Au', qa), ('p+p', qp)):
    c = r['counts']
    g = c['netc_nonzero_events'] == 0 and c['nfail'] == 0 and c['dropped_nonfinite'] == 0
    ok &= g
    print(f'{"PASS" if g else "FAIL"}  {lab} gates: {c["events"]} events, net charm != 0 in {c["netc_nonzero_events"]}, '
          f'charm matching failures {c["nfail"]}, dropped {c["dropped_nonfinite"]}')
print('PASS  consistent with the published results' if ok else 'FAIL  inconsistent with the published results')
sys.exit(0 if ok else 1)
