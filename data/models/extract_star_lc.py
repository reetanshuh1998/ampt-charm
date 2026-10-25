#!/usr/bin/env python3
"""Extract the model curves of STAR PRL 124 172301 (arXiv:1910.14628), Fig. 2(b) (Lambda_c/D0 vs pT, 10-80 %)
and Fig. 3 (vs Npart, 3 < pT < 6 GeV/c), from the PDF vector paths -> star_lc_models.json.

pdftocairo -svg -f 5 -l 5 1910.14628.pdf p5.svg ; pdftocairo -svg -f 6 -l 6 1910.14628.pdf p6.svg
usage: extract_star_lc.py p5.svg p6.svg OUT.json

Axes are calibrated on the centres of STAR's own markers (red filled circles), fitted linearly to the HEPData
values (ins1762441, Figures 2 and 3); the residual of the third point is the calibration check (asserted
< 0.003 in pT / ratio).  Path indices were identified by colour and dash pattern against the legend.
"""
import json
import re
import sys

import numpy as np


def load(fn):
    s = open(fn).read()
    s = s[s.find('</defs>'):]
    out = {}
    for i, p in enumerate(re.findall(r'<path ([^>]*?)/?>', s)):
        a = dict(re.findall(r'([\w-]+)="([^"]*)"', p))
        n = list(map(float, re.findall(r'-?\d+\.?\d*(?:e-?\d+)?', a.get('d', ''))))
        if len(n) < 4:
            continue
        P = np.array(list(zip(n[0::2], n[1::2])))
        m = re.findall(r'-?\d+\.?\d*(?:e-?\d+)?', a.get('transform', ''))
        if m:
            A, B, C, D, E, F = map(float, m)
            P = np.c_[A * P[:, 0] + C * P[:, 1] + E, B * P[:, 0] + D * P[:, 1] + F]
        out[i] = dict(P=P, stroke=a.get('stroke', ''))
    return out


def centre(P):
    return 0.5 * (P[:, 0].min() + P[:, 0].max()), 0.5 * (P[:, 1].min() + P[:, 1].max())


def calib(markers, want):
    """linear page->data maps from marker centres; returns fx, fy and the max residual"""
    c = np.array([centre(m) for m in markers])
    w = np.array(want)
    kx = np.polyfit(c[:, 0], w[:, 0], 1)
    ky = np.polyfit(c[:, 1], w[:, 1], 1)
    res = max(np.max(np.abs(np.polyval(kx, c[:, 0]) - w[:, 0])), np.max(np.abs(np.polyval(ky, c[:, 1]) - w[:, 1])))
    return (lambda x: np.polyval(kx, x)), (lambda y: np.polyval(ky, y)), float(res)


def curve(P, fx, fy, nd=4):
    a = P[np.argsort(P[:, 0])]
    return dict(x=[round(float(v), 3) for v in fx(a[:, 0])], y=[round(float(v), nd) for v in fy(a[:, 1])])


def band_centre(P, fx, fy, nb=40):
    """filled band (closed polygon): mid-line and half-width in nb x bins"""
    x, y = fx(P[:, 0]), fy(P[:, 1])
    e = np.linspace(x.min(), x.max(), nb + 1)
    xs, lo, hi = [], [], []
    for a, b in zip(e[:-1], e[1:]):
        s = (x >= a) & (x <= b)
        if s.sum() >= 2:
            xs.append(0.5 * (a + b)); lo.append(y[s].min()); hi.append(y[s].max())
    return dict(x=[round(v, 3) for v in xs], y=[round(0.5 * (a + b), 4) for a, b in zip(lo, hi)],
                lo=[round(v, 4) for v in lo], hi=[round(v, 4) for v in hi])


P5, P6 = load(sys.argv[1]), load(sys.argv[2])
# ---- Fig. 2(b): markers 309, 308, 307 = pT 2.92, 4.07, 5.94 (HEPData Figure 2)
fx, fy, r2 = calib([P5[309]['P'], P5[308]['P'], P5[307]['P']], [(2.92, 1.1903), (4.07, 0.7995), (5.94, 0.4818)])
assert r2 < 0.003, r2
# PYTHIA (CR) is drawn as many 2-point segments in magenta: collect them
cr = np.vstack([p['P'] for i, p in P5.items() if p['stroke'] == 'rgb(100%, 0%, 100%)' and len(p['P']) == 2
                and p['P'][:, 1].max() < 311 and p['P'][:, 1].min() > 177])
cr = cr[fy(cr[:, 1]) < 1.0]          # drop the legend sample line (drawn at ratio ~3.3)
cr = cr[np.argsort(cr[:, 0])]
_, iu = np.unique(np.round(cr[:, 0], 1), return_index=True)
fig2 = {'Ko et al.: three quark (0-5 %)': curve(P5[272]['P'], fx, fy),
        'Ko et al.: di-quark (0-5 %)': curve(P5[273]['P'], fx, fy),
        'Ko et al.: with flow (0-10 %)': curve(P5[278]['P'], fx, fy),
        'Catania, coal.+frag. (10-80 %)': curve(P5[274]['P'], fx, fy),
        'Catania, coal. (10-80 %)': curve(P5[275]['P'], fx, fy),
        'Tsinghua (10-80 %)': curve(P5[276]['P'], fx, fy),
        'Rapp et al. (0-20 %)': band_centre(P5[277]['P'], fx, fy),
        'Cao et al. (10-80 %)': curve(P5[279]['P'], fx, fy),
        'PYTHIA 8': curve(P5[185]['P'], fx, fy),
        'PYTHIA 8, CR': curve(cr[iu], fx, fy)}
# ---- Fig. 3: markers 56, 55, 54 = Npart 30.1, 119.6, 279.8 (HEPData Figure 3)
gx, gy, r3 = calib([P6[56]['P'], P6[55]['P'], P6[54]['P']], [(30.1, 0.4991), (119.6, 0.7948), (279.8, 1.0772)])
assert r3 < 0.003 * 100, r3      # x residual is in Npart units
cy = lambda P: float(gy(P[:, 1].mean() if len(P) == 3 else centre(P)[1]))
fig3 = {'Catania (coal.)': curve(P6[107]['P'], gx, gy), 'Catania (coal.+frag.)': curve(P6[108]['P'], gx, gy),
        'PYTHIA 8 (p+p, marker)': round(cy(P6[126]['P']), 4), 'PYTHIA 8, CR (p+p, marker)': round(cy(P6[127]['P']), 4)}
out = dict(source='STAR, PRL 124 172301 (2020), arXiv:1910.14628, Fig. 2(b) (Lambda_c/D0 vs pT; data 10-80 %, model '
                  'centralities as labelled) and Fig. 3 (vs Npart, 3 < pT < 6 GeV/c); curves from the PDF vector paths '
                  '(extract_star_lc.py). Rapp et al. is a band: y = mid-line, lo/hi = edges.',
           calibration_max_residual=dict(fig2=round(r2, 5), fig3_npart_or_ratio=round(r3, 5)),
           star_quoted_chi2=dict(note='STAR text, Fig. 2(b), with finite-bin-width treatment: PYTHIA 20.7/3, PYTHIA CR '
                                      '8.2/3, Catania coal. reduced chi2 26.1, Ko with flow reduced chi2 4.8; others '
                                      '"consistent within uncertainties"'),
           fig2b=fig2, fig3=fig3)
json.dump(out, open(sys.argv[3], 'w'), indent=1)
print('calibration residuals', out['calibration_max_residual'])
for k, v in fig2.items():
    print(f'{k:32s} pT {v["x"][0]:.2f}-{v["x"][-1]:.2f}  at 2.92/4.07/5.94:', np.round(np.interp([2.92, 4.07, 5.94], v['x'], v['y']), 3))
print('fig3', {k: (v if not isinstance(v, dict) else v['y']) for k, v in fig3.items()})
