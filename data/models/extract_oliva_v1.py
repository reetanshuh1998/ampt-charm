#!/usr/bin/env python3
"""Extract the D0+D0bar v1(y) curves of Oliva, Plumari, Greco, JHEP 05 (2021) 034 (arXiv:2009.11066), Fig. 8
(Au+Au 200 GeV, b = 9 fm, pT > 1.5 GeV/c), from the PDF vector paths -> oliva2021_d0v1.json.

pdftocairo -svg -f 9 -l 9 2009.11066.pdf p9.svg ; usage: extract_oliva_v1.py p9.svg OUT.json

Axes: linear maps fitted to the centres of the four STAR diamonds drawn in each panel (HEPData ins1733225,
Figure 3a: y = -0.6, -0.2, 0.2, 0.6; <v1> = 0.069, 0.018, -0.030, -0.025); max residual asserted < 0.002.
Paths: top panel 309 = eta_m 1.1 (solid blue), 310 = eta_m 1.2 (dashed red); bottom panel 479 = charm NOT tilted
(eta_m 1.1, solid), 480 = charm tilted like the bulk (dashed).
"""
import json
import re
import sys

import numpy as np

s = open(sys.argv[1]).read()
s = s[s.find('</defs>'):]
PS = re.findall(r'<path ([^>]*?)/?>', s)


def pts(i):
    a = dict(re.findall(r'([\w-]+)="([^"]*)"', PS[i]))
    n = list(map(float, re.findall(r'-?\d+\.?\d*(?:e-?\d+)?', a['d'])))
    P = np.array(list(zip(n[0::2], n[1::2])))
    m = re.findall(r'-?\d+\.?\d*(?:e-?\d+)?', a.get('transform', ''))
    if m:
        A, B, C, D, E, F = map(float, m)
        P = np.c_[A * P[:, 0] + C * P[:, 1] + E, B * P[:, 0] + D * P[:, 1] + F]
    return P


STAR = np.array([(-0.6, 0.069), (-0.2, 0.018), (0.2, -0.030), (0.6, -0.025)])


def panel(markers, curves):
    c = np.array([[np.mean([pts(i)[:, 0].min(), pts(i)[:, 0].max()]), np.mean([pts(i)[:, 1].min(), pts(i)[:, 1].max()])]
                  for i in markers])
    kx, ky = np.polyfit(c[:, 0], STAR[:, 0], 1), np.polyfit(c[:, 1], STAR[:, 1], 1)
    r = max(np.abs(np.polyval(kx, c[:, 0]) - STAR[:, 0]).max(), np.abs(np.polyval(ky, c[:, 1]) - STAR[:, 1]).max())
    assert r < 0.002, r
    out = {}
    for nm, i in curves.items():
        P = pts(i)
        P = P[np.argsort(P[:, 0])]
        out[nm] = dict(y=[round(float(v), 3) for v in np.polyval(kx, P[:, 0])],
                       v1=[round(float(v), 5) for v in np.polyval(ky, P[:, 1])])
    return out, float(r)


top, r1 = panel((285, 286, 287, 288), {'eta_m 1.1': 309, 'eta_m 1.2': 310})
bot, r2 = panel((455, 456, 457, 458), {'charm NOT tilted (eta_m 1.1)': 479, 'charm tilted like the bulk': 480})
out = dict(source='L. Oliva, S. Plumari, V. Greco, JHEP 05 (2021) 034, arXiv:2009.11066, Fig. 8: '
                  'Catania Boltzmann transport (quasi-particle bulk with a tilted initial condition fitted to charged-'
                  'hadron v1, eta_m 1.1/1.2 bracket 5-40 % and 40-80 %), b = 9 fm, pT > 1.5 GeV/c, D0+D0bar averaged; '
                  'curves from the PDF vector paths (extract_oliva_v1.py).',
           calibration_max_residual=[round(r1, 5), round(r2, 5)], top=top, bottom=bot)
json.dump(out, open(sys.argv[2], 'w'), indent=1)
yc = np.array([-0.6, -0.2, 0.2, 0.6])
for pan in (top, bot):
    for k, v in pan.items():
        m = np.interp(yc, v['y'], v['v1'])
        print(f'{k:30s} v1 at -0.6..0.6: {np.round(m, 4)}  slope (through 0): {np.sum(yc * m) / np.sum(yc * yc):+.4f}')
print('calibration residuals', out['calibration_max_residual'])
