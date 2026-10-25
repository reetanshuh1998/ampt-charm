#!/usr/bin/env python3
"""Extract the LBT and Duke D0 R_AA curves (0-10 %) of STAR PRC 99 034908 (arXiv:1812.10224), Fig. 32(a), from
the PDF vector paths -> star_raa_models.json.

pdftocairo -svg -f 24 -l 24 1812.10224.pdf p24.svg
usage: extract_star_raa.py p24.svg OUT.json

Axes: least-squares linear maps from the centres of STAR's 11 D0 markers (filled: path 75, open: path 84, compound
paths split per sub-path) to the HEPData values (ins1711377, RAA table, 0-10 %); the max residual is asserted
< 0.005 (pT in GeV/c, R_AA).  Path 52 = LBT (black dotted), path 53 = Duke (blue dashed).
"""
import json
import re
import sys

import numpy as np

s = open(sys.argv[1]).read()
s = s[s.find('</defs>'):]
PS = re.findall(r'<path ([^>]*?)/?>', s)


def subpaths(idx):
    a = dict(re.findall(r'([\w-]+)="([^"]*)"', PS[idx]))
    m = re.findall(r'-?\d+\.?\d*(?:e-?\d+)?', a.get('transform', ''))
    A, B, C, D, E, F = map(float, m) if m else (1, 0, 0, 1, 0, 0)
    out = []
    for sub in re.split(r'M', a['d'])[1:]:
        n = list(map(float, re.findall(r'-?\d+\.?\d*(?:e-?\d+)?', sub)))
        P = np.array(list(zip(n[0::2], n[1::2])))
        out.append(np.c_[A * P[:, 0] + C * P[:, 1] + E, B * P[:, 0] + D * P[:, 1] + F])
    return out


def centres(idx):
    return [(0.5 * (P[:, 0].min() + P[:, 0].max()), 0.5 * (P[:, 1].min() + P[:, 1].max()))
            for P in subpaths(idx) if np.ptp(P[:, 0]) > 1]


# HEPData ins1711377 RAA 0-10 %: (<pT>, R_AA)
HD = [(0.242, 0.312), (0.728, 0.386), (1.22, 0.482), (1.71, 0.622), (2.21, 0.684), (2.71, 0.70), (3.36, 0.539),
      (4.36, 0.473), (5.37, 0.413), (6.59, 0.301), (8.64, 0.184)]
c = np.array(sorted(centres(75) + centres(84)))
assert len(c) == 11, len(c)
w = np.array(HD)
kx, ky = np.polyfit(c[:, 0], w[:, 0], 1), np.polyfit(c[:, 1], w[:, 1], 1)
res = (float(np.max(np.abs(np.polyval(kx, c[:, 0]) - w[:, 0]))), float(np.max(np.abs(np.polyval(ky, c[:, 1]) - w[:, 1]))))
assert max(res) < 0.005, res


def curve(idx):
    P = np.vstack(subpaths(idx))
    P = P[np.argsort(P[:, 0])]
    return dict(pt=[round(float(v), 3) for v in np.polyval(kx, P[:, 0])],
                raa=[round(float(v), 4) for v in np.polyval(ky, P[:, 1])])


out = dict(source='STAR, PRC 99 034908 (2019), arXiv:1812.10224, Fig. 32(a): D0 R_AA 0-10 % with the LBT and Duke '
                  'model calculations (refs. 48-50 there); curves from the PDF vector paths (extract_star_raa.py).',
           calibration_max_residual_pt_raa=[round(v, 4) for v in res],
           RAA_D0_0_10={'LBT': curve(52), 'Duke': curve(53)})
json.dump(out, open(sys.argv[2], 'w'), indent=1)
print('calibration residuals (pT, R_AA):', out['calibration_max_residual_pt_raa'])
for k, v in out['RAA_D0_0_10'].items():
    print(f'{k:5s} pT {v["pt"][0]:.2f}-{v["pt"][-1]:.2f}  R_AA at 1.2/2.2/3.4/5.4 GeV:',
          np.round(np.interp([1.22, 2.21, 3.36, 5.37], v['pt'], v['raa']), 3))
