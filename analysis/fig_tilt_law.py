#!/usr/bin/env python3
"""fig_tilt_law.py OUT_STEM : headline figure (PAPER_PLAN Fig. 7, STEP_LOG 98-101).

(a) c-quark v1 vs the charm-to-bulk tilt ratio r = eps1(c)/eps1(u); (b) D0 dv1/dy (10-80 %, pT > 1.5 GeV) vs r with
STAR.  Duke-Owens points: one 600k sample re-weighted (none / HIJING x 0.4 / full HIJING); CTEQ6.1M points: free and
EPS09s, the EPS09s point with its Hessian 90 % CL band (observation/results/tables/EPS09S_BAND_DM25.md) added in quadrature to the statistical
error.  Line: linear fit to the Duke-Owens points.
2026-10-25: EPS09s point = final 2.4M sample (dm25_main), plus public AMPT (dm25_pub) as an independent test of the law
(backup fig_tilt_law_v1.py.bak)."""
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

D = str(__import__('pathlib').Path(__file__).resolve().parents[1] / 'observation' / 'results' / 'json')
OUT = sys.argv[1]


def odd(r, key, lo=0.5, hi=2.5):
    yq = np.array(r['edges']['Y_Q']); c = 0.5 * (yq[1:] + yq[:-1])
    v, e = np.array(r[key]['val'], float), np.array(r[key]['err'], float)
    s = (np.abs(c) > lo) & (np.abs(c) < hi) & np.isfinite(v) & (e > 0)
    w = 1 / e[s] ** 2
    return np.sum(w * np.sign(c[s]) * v[s]) / np.sum(w), 1 / np.sqrt(np.sum(w))


def point(fn):
    r = json.load(open(fn))
    eu, ec = odd(r, 'eps1_uq_all'), odd(r, 'eps1_cq_all')
    rat = ec[0] / eu[0], abs(ec[0] / eu[0]) * np.hypot(ec[1] / ec[0], eu[1] / eu[0])
    return rat, odd(r, 'v1_cq_all'), tuple(r['slope_D0avg_10-80'])


DO = [point(f'{D}/{f}.json') for f in ('dm23_ab__none', 'dm23_ab', 'dm23_ab__hij10')]
CT = [point(f'{D}/dm24_e_ab__cteq.json'), point(f'{D}/dm25_main.json')]   # free CTEQ (600k), x EPS09s (2.4M)
PUB = point(f'{D}/dm25_pub.json')                                                  # public AMPT 2.26t9b (600k)
BAND = {'r': 0.11, 'v1c': 0.004, 'd0': 0.0064}             # EPS09s Hessian on 2.4M (symmetrised), EPS09S_BAND_DM25.md

INK, MUTED, GRID = '#1f2328', '#59636e', '#d8dee4'
C_DO, C_CT = '#2a6fdb', '#d4572a'
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False})
fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.8), constrained_layout=True)
x = np.linspace(0, 1, 50)
for k, (a, key, lab) in enumerate(((ax[0], 1, r'$c$-quark $v_1$ (0.5 < |y| < 2.5)'),
                                   (ax[1], 2, r'D$^0$ d$v_1$/d$y$ (10–80 %)'))):
    r0 = np.array([p[0][0] for p in DO]); y0 = np.array([p[key][0] for p in DO])
    fit = np.polyfit(r0, y0, 1)
    a.plot(x, np.polyval(fit, x), color=C_DO, lw=1.5, alpha=0.6, zorder=1)
    a.errorbar(r0, y0, xerr=[p[0][1] for p in DO], yerr=[p[key][1] for p in DO], fmt='o', ms=7, color=C_DO,
               mec='white', mew=1.5, elinewidth=1.5, capsize=0, label='Duke–Owens (none, HIJ×0.4, full HIJ)', zorder=3)
    bx = [0, BAND['r']]; by = [0, BAND['v1c' if key == 1 else 'd0']]
    a.errorbar([p[0][0] for p in CT], [p[key][0] for p in CT],
               xerr=[np.hypot(p[0][1], b) for p, b in zip(CT, bx)],
               yerr=[np.hypot(p[key][1], b) for p, b in zip(CT, by)], fmt='s', ms=7, color=C_CT, mfc='white',
               mew=2, elinewidth=1.5, capsize=0, label='CTEQ6.1M (free, ×EPS09s)', zorder=4)
    a.errorbar([PUB[0][0]], [PUB[key][0]], xerr=[PUB[0][1]], yerr=[PUB[key][1]], fmt='^', ms=7, color='#59636e',
               mfc='white', mew=1.8, elinewidth=1.5, capsize=0, label='public AMPT 2.26t9b', zorder=4)
    if key == 2:
        a.axhspan(-0.080 - 0.023, -0.080 + 0.023, color=GRID, alpha=0.7, lw=0, zorder=0)
        a.axhline(-0.080, color=MUTED, lw=1, zorder=0)
        a.text(0.02, -0.080, 'STAR −0.080 ± 0.023', color=MUTED, va='bottom', fontsize=9)
    a.axhline(0, color=GRID, lw=1, zorder=0)
    a.set_xlim(0, 1); a.set_xlabel(r'charm-to-bulk tilt ratio $\varepsilon_1(c)/\varepsilon_1(u)$')
    a.set_ylabel(lab); a.grid(axis='y', color=GRID, lw=0.6)
    a.text(0.02, 0.97, '(a)' if key == 1 else '(b)', transform=a.transAxes, va='top', fontweight='bold', color=INK)
ax[0].legend(frameon=False, loc='lower right', fontsize=9)
fa = np.polyfit([p[0][0] for p in DO], [p[1][0] for p in DO], 1)
ax[0].text(0.98, 0.30, f'fit: $v_1(c)$ = {fa[1]:.4f} + {fa[0]:.4f} r', transform=ax[0].transAxes, ha='right',
           color=MUTED, fontsize=9)
for ext in ('pdf', 'png'):
    fig.savefig(f'{OUT}.{ext}', dpi=200)
print('wrote', OUT)
