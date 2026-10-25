#!/usr/bin/env python3
"""Extract the model curves of STAR PRL 127 092301 (arXiv:2101.11793), Figs. 3 and 4, from the
PDF vector paths (pdftocairo -svg of pages 6 and 7) -> star_ds_models.json.

Axes are calibrated on the tick marks and checked on the systematic brackets of the STAR points,
which must reproduce HEPData ins1843268 (checked: 0-10, 10-20, 20-40, 40-80 % first bins and the
0-10 % point of Fig. 4, all to 3 decimals).
usage: extract_star_ds.py p6.svg p7.svg OUT.json
"""
import json
import re
import sys

import numpy as np


def paths(fn):
    s = open(fn).read()
    out = []
    for p in re.findall(r'<path ([^>]*?)/?>', s[s.find('</defs>'):]):
        a = dict(re.findall(r'([\w-]+)="([^"]*)"', p))
        n = list(map(float, re.findall(r'-?\d+\.?\d*(?:e-?\d+)?', a.get('d', ''))))
        out.append(np.array(list(zip(n[0::2], n[1::2]))))
    return out


P6, P7 = paths(sys.argv[1]), paths(sys.argv[2])
# Fig. 3: x major ticks 145.6 .. 522.2 = 1 .. 7 GeV; y major ticks (a) 395.6 .. 574.7 = 0.2 .. 0.8,
# (b) 118.8 .. 294.8 = 0.2 .. 0.8 (local coordinates of the figure's transform, y up)
pt = lambda x: 1 + (x - 145.6) / 62.77
ya = lambda y: 0.2 + (y - 395.6) / 298.5
yb = lambda y: 0.2 + (y - 118.8) / 293.3
# Fig. 4: x ticks 85.5 / 203.3 = 0 / 100; y ticks 64.9 / 385.5 = 0 / 0.8
npart = lambda x: (x - 85.5) * 100 / 117.8
y4 = lambda y: (y - 64.9) / 400.75


def curve(arr, fx, fy):
    a = arr[np.argsort(arr[:, 0])]
    return dict(x=[round(float(v), 3) for v in fx(a[:, 0])], y=[round(float(v), 4) for v in fy(a[:, 1])])


def check(lo, hi, fy, want):
    got = fy(0.5 * (lo + hi))
    assert abs(got - want) < 0.003, (got, want)
    return round(float(got), 4)


chk = dict(fig3a_0_10_first=check(397.0, 428.1, ya, 0.2568), fig3a_10_20_first=check(427.8, 455.3, ya, 0.3539),
           fig3b_40_80_first=check(162.3, 196.4, yb, 0.4063), fig3b_20_40_first=check(160.3, 181.6, yb, 0.3777),
           fig4_0_10=check(153.8, 187.1, y4, 0.2636))
out = dict(
    source='STAR, PRL 127 092301 (2021), arXiv:2101.11793, Figs. 3 and 4; curves from the PDF vector paths '
           '(extract_star_ds.py); x = pT (GeV/c) or <Npart>, y = Ds/D0 = (Ds+ + Ds-)/(D0 + D0bar).',
    calibration_check_vs_hepdata=chk,
    fig3a_0_10={'Tsinghua (seq. coal.)': curve(P6[92], pt, ya), 'Catania (coal.)': curve(P6[42], pt, ya),
                'Catania (coal.+frag.)': curve(P6[43], pt, ya), 'Cao, Ko': curve(P6[88], pt, ya),
                'He, Rapp (0-20 %)': curve(P6[90], pt, ya), 'PYTHIA 8 p+p': curve(P6[41], pt, ya)},
    fig3b={'Tsinghua (seq. coal.), 20-40 %': curve(P6[238], pt, yb), 'Tsinghua (seq. coal.), 40-80 %': curve(P6[240], pt, yb),
           'PYTHIA 8 p+p': curve(P6[236], pt, yb)},
    fig4={'Tsinghua (seq. coal.)': curve(P7[38], npart, y4), 'STAR constant fit': round(float(y4(215.7)), 4)},
    thermus_10_40_integrated=0.36,
    star_integrated_10_40=[0.42, 0.04, 0.11])
json.dump(out, open(sys.argv[3], 'w'), indent=1)
print(json.dumps(chk), 'constant fit', out['fig4']['STAR constant fit'])
for k, v in out['fig3a_0_10'].items():
    print(f'{k:24s} pT {v["x"][0]:.2f}-{v["x"][-1]:.2f}  Ds/D0 at ~2,4,6 GeV:',
          np.round(np.interp([2, 4, 6], v['x'], v['y']), 3))
