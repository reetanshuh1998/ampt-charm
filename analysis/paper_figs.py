#!/usr/bin/env python3
"""paper_figs.py OUTDIR : figures of the open-charm AMPT paper from the final samples (dm25, observation/results/json/), STAR
HEPData and the digitised published curves (data/models/, provenance and calibration in each file).

Footing (docs/CLAIMS_AUDIT.md): every model point is the same estimator and kinematic window as the STAR point it is
drawn against (D0 v2 = scalar product in STAR's eta-sub geometry; v1 vs the reaction plane, |y| < 0.8, pT > 1.5 GeV,
10-80 %; R_AA from per-c-quark yields AA / p+p, |y| < 1; Lambda_c/D0 in STAR's pT bins).  Transport output: markers
with jackknife errors; published curves: lines.  EPS09s band: Hessian 90 % CL (observation/results/tables/EPS09S_BAND_DM25.md).
Writes OUTDIR/fig_*.pdf and .png.
"""
import csv
import glob
import json
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

D = str(__import__('pathlib').Path(__file__).resolve().parents[1])   # repository root
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)

# reference categorical order (dataviz palette), data in black; published models in fixed hues, dashed
INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e6e6e3'
C_FW, C_PUB = '#2a78d6', '#eb6834'
CM = ['#1baf7a', '#eda100', '#e87ba4', '#6250d6', '#008300', '#e34948', '#8a8a85']
DASH = [(0, (5, 2)), (0, (1, 1.2)), (0, (6, 2, 1, 2)), (0, (3, 1, 1, 1, 1, 1)), (0, (8, 2.5)), (0, (2, 2)), '-']
plt.rcParams.update({'font.size': 9, 'font.family': 'serif', 'mathtext.fontset': 'stix', 'axes.spines.top': True,
                     'axes.spines.right': True, 'axes.edgecolor': INK, 'axes.linewidth': 0.7, 'xtick.direction': 'in',
                     'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True, 'legend.frameon': False,
                     'legend.fontsize': 7, 'lines.linewidth': 1.4, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02})
W1, W2 = 3.4, 7.0           # PRD single / double column (in)


def J(f):
    return json.load(open(f'{D}/{f}'))


def hep(pattern, block=0):
    fn = glob.glob(f'{D}/data/hepdata/{pattern}')[0]
    blocks, cur = [], []
    for row in csv.reader(open(fn)):
        try:
            if not row:
                raise ValueError
            cur.append([float(x) for x in row])
        except ValueError:
            if cur:
                blocks.append(np.array(cur)); cur = []
    if cur:
        blocks.append(np.array(cur))
    return blocks[block]


R = dict(fw=J('observation/results/json/dm25_main.json'), pp=J('observation/results/json/dm25_ppmain.json'), pub=J('observation/results/json/dm25_pub.json'),
         pppub=J('observation/results/json/dm25_pppub.json'), c1=J('observation/results/json/dm25_c1.json'), ppc1=J('observation/results/json/dm25_ppc1.json'))
FAN, MV2, MV1 = J('data/models/fan2026_improved_ampt.json'), J('data/models/star_d0v2_fig4_models.json'), J('data/models/star_d0v1_fig3_models.json')
MLC, MRAA, OLI, MDS = (J('data/models/star_lc_models.json'), J('data/models/star_raa_models.json'),
                       J('data/models/oliva2021_d0v1.json'), J('data/models/star_ds_models.json'))
LAB = dict(fw='This work', pub='Public AMPT', c1='This work, $c_1$ (Jacobi)')
EPS = {}
for line in open(f'{D}/observation/results/tables/EPS09S_BAND_DM25.md'):
    m = re.match(r'\| (R_AA \S+ % pT [\d.]+|D0 dv1/dy[^|]*|v1\(c\)|tilt ratio[^|]*) \| ([-\d.]+) \| ([\d.]+) \| ([\d.]+) \|', line)
    if m:
        EPS[m.group(1).strip()] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))


def val(r, k):
    return np.array(r[k]['val'], float), np.array(r[k]['err'], float)


def data(ax, x, y, st, sy, xw, label='STAR', ms=3.6):
    for xi, yi, si in zip(x, y, sy):
        ax.add_patch(Rectangle((xi - xw, yi - si), 2 * xw, 2 * si, fill=False, ec=INK, lw=0.6, zorder=5))
    ax.errorbar(x, y, st, fmt='o', color=INK, ms=ms, lw=0.9, capsize=0, label=label, zorder=6)


def model(ax, x, v, e, k, dx=0.0, label=None, ms=3.8):
    ok = np.isfinite(v) & np.isfinite(e)
    mk = dict(fw='s', pub='^', c1='D')[k]
    col = dict(fw=C_FW, pub=C_PUB, c1=C_FW)[k]
    ax.errorbar(np.asarray(x)[ok] + dx, v[ok], e[ok], fmt=mk, color=col, mfc=col if k == 'fw' else 'white',
                ms=ms, lw=1.0, capsize=0, label=label or LAB[k], zorder=4)


def curve(ax, x, y, i, label, lw=1.1):
    ax.plot(x, y, color=CM[i % len(CM)], ls=DASH[i % len(DASH)], lw=lw, label=label, zorder=2)


def tag(ax, s, x=0.04, y=0.95, **kw):
    ax.text(x, y, s, transform=ax.transAxes, va='top', **kw)


def save(fig, nm):
    for ext in ('pdf', 'png'):
        fig.savefig(f'{OUT}/{nm}.{ext}', dpi=250)
    plt.close(fig)
    print('wrote', nm)


# =========================================================== Fig. 1  p+p baseline
st = hep('ins1711377/*/D^0spectrainppcollisions.csv')
pt_s, y_s = st[:, 0], st[:, 1]
e_st, e_sy = st[:, 2], 0.5 * (st[:, 4] + np.abs(st[:, 5]))
fig, ax = plt.subplots(2, 1, figsize=(W1, 3.6), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 1], hspace=0))
for k, pk, dx in (('fw', 'pp', -0.05), ('pub', 'pppub', 0.05)):
    r = R[pk]
    ed = np.array(r['edges']['PT_RAA']); c = 0.5 * (ed[1:] + ed[:-1]); w = np.diff(ed)
    v, e = val(r, 'dN_D0_mb_perc')
    dens, ed_ = v / w, e / w
    ok = dens > 0
    m = np.exp(np.interp(pt_s, c[ok], np.log(dens[ok]))) / (2 * np.pi * pt_s) / 2.0
    me = m * np.interp(pt_s, c[ok], ed_[ok] / dens[ok])
    wt = 1 / (np.hypot(e_st, e_sy) ** 2 + (me * y_s / m) ** 2)
    a = np.sum(wt * y_s * m) / np.sum(wt * m * m)
    chi2 = np.sum(wt * (y_s - a * m) ** 2)
    # full model spectrum, same normalisation, at the model bin centres
    full = a * dens / (2 * np.pi * c) / 2.0
    ax[0].plot(c[ok & (c < 8)], full[ok & (c < 8)], color=C_FW if k == 'fw' else C_PUB, ls='-' if k == 'fw' else '--',
               lw=1.2, label=f'{LAB[k]}, $\\chi^2$ = {chi2:.1f}')
    ax[1].errorbar(pt_s + dx, a * m / y_s, np.hypot(e_st, e_sy) / y_s, fmt='s' if k == 'fw' else '^',
                   color=C_FW if k == 'fw' else C_PUB, mfc=C_FW if k == 'fw' else 'white', ms=3.5, lw=0.9, capsize=0)
data(ax[0], pt_s, y_s, e_st, e_sy, 0.12, label='STAR p+p, D$^0$ ($|y|<1$)')
ax[0].set(yscale='log', ylabel=r'$\frac{1}{2\pi p_T}\frac{d^2N}{dp_Tdy}$ (a.u.)', xlim=(0, 7))
ax[0].legend(loc='lower left')
tag(ax[0], 'p+p, $\\sqrt{s}$ = 200 GeV\nnormalisation fitted\n$\\chi^2$: 6 points, 5 d.o.f.', x=0.45, y=0.95)
ax[1].axhline(1, color=MUTED, lw=0.7)
ax[1].set(xlabel='$p_T$ (GeV/$c$)', ylabel='model/STAR', ylim=(0, 2.4))
save(fig, 'fig_pp')

# =========================================================== Fig. 2  R_AA
fig, ax = plt.subplots(1, 3, figsize=(W2, 2.6), sharey=True, gridspec_kw=dict(wspace=0))
for j, cl in enumerate(('0-10', '10-40', '40-80')):
    a_ = ax[j]
    t = hep('ins1711377/*/RAA.csv', j)
    sel = t[:, 3] > 0.9
    f = FAN['RAA_D0'][cl]
    curve(a_, f['pt'], f['raa'], 0, 'Improved AMPT (Fan et al.)')
    if cl == '0-10':
        for i, (nm, c) in enumerate(MRAA['RAA_D0_0_10'].items()):
            curve(a_, c['pt'], c['raa'], 3 + i, nm)
    for k, pk, dx in (('fw', 'pp', 0.0), ('pub', 'pppub', 0.07)):
        A_, Ae = val(R[k], f'dN_D0_{cl}_perc'); P_, Pe = val(R[pk], 'dN_D0_mb_perc')
        ed = np.array(R[k]['edges']['PT_RAA']); x = 0.5 * (ed[1:] + ed[:-1])
        x = np.array([t[np.argmin(np.abs(t[:, 1] - lo)), 3] if np.min(np.abs(t[:, 1] - lo)) < 2e-3 else xc
                      for lo, xc in zip(ed[:-1], x)])
        with np.errstate(invalid='ignore', divide='ignore'):
            rr, re_ = A_ / P_, A_ / P_ * np.hypot(Ae / A_, Pe / P_)
        s = x < 8
        if k == 'fw':
            band = np.array([EPS.get(f'R_AA {cl} % pT {xx:.2f}', (np.nan, np.nan, np.nan)) for xx in x[s]])
            okb = np.isfinite(band[:, 0])
            a_.fill_between(x[s][okb], rr[s][okb] - band[okb, 2], rr[s][okb] + band[okb, 1], color=C_FW, alpha=0.18,
                            lw=0, label='This work, EPS09s 90 % CL', zorder=1)
        model(a_, x[s], rr[s], re_[s], k, dx)
    data(a_, t[:, 3], t[:, 4], t[:, 5], 0.5 * (np.abs(t[:, 7]) + np.abs(t[:, 8])), 0.1,
         label='STAR D$^0$' if j == 0 else None)
    a_.axhline(1, color=MUTED, lw=0.6)
    a_.set(xlabel='$p_T$ (GeV/$c$)', xlim=(0, 7.6), ylim=(0, 1.75))
    tag(a_, f'{cl.replace("-", "–")} %', x=0.72)
ax[0].set_ylabel('$R_{AA}$ (D$^0$, $|y|<1$)')
h, l = [], []
for a_ in ax:
    hh, ll = a_.get_legend_handles_labels()
    for x_, y_ in zip(hh, ll):
        if y_ not in l:
            h.append(x_); l.append(y_)
fig.legend(h, l, loc='lower center', ncol=4, fontsize=6.6, bbox_to_anchor=(0.5, -0.2))
save(fig, 'fig_raa')

# =========================================================== Fig. 3  Lambda_c / D0
fig, ax = plt.subplots(1, 2, figsize=(W2, 2.9), gridspec_kw=dict(wspace=0.25))
a_ = ax[0]
sel_models = [('Catania, coal.+frag. (10-80 %)', 'Catania (coal.+frag.)'), ('Tsinghua (10-80 %)', 'Tsinghua'),
              ('Ko et al.: di-quark (0-5 %)', 'Ko et al. (di-quark, 0–5 %)'), ('Cao et al. (10-80 %)', 'Cao et al.'),
              ('PYTHIA 8, CR', 'PYTHIA 8 CR (p+p)')]
r_ = MLC['fig2b']['Rapp et al. (0-20 %)']
a_.fill_between(r_['x'], r_['lo'], r_['hi'], color=CM[6], alpha=0.25, lw=0, label='He–Rapp (0–20 %)')
curve(a_, FAN['LcD0']['AuAu_10-80']['pt'], FAN['LcD0']['AuAu_10-80']['ratio'], 0, 'Improved AMPT (Fan et al.)')
for i, (k, lab) in enumerate(sel_models):
    curve(a_, MLC['fig2b'][k]['x'], MLC['fig2b'][k]['y'], i + 1, lab, lw=1.0)
ed = np.array(R['fw']['edges']['PT_LC']); xl = np.array([2.92, 4.07, 5.94])
for k, dx in (('fw', 0.0), ('pub', 0.12), ('c1', -0.12)):
    A_, Ae = val(R[k], 'dN_Lc_lcbins_10-80_perc'); B_, Be = val(R[k], 'dN_D0_lcbins_10-80_perc')
    model(a_, xl, A_ / B_, A_ / B_ * np.hypot(Ae / A_, Be / B_), k, dx)
t = hep('ins1762441/*/Figure2.csv')
data(a_, t[:, 0], t[:, 1], t[:, 2], t[:, 4], 0.1, label='STAR 10–80 %')
a_.set(xlabel='$p_T$ (GeV/$c$)', ylabel=r'$\Lambda_c^+/D^0$', xlim=(0, 8), ylim=(0, 2.6))
a_.legend(loc='upper right', fontsize=5.8, ncol=2)
a_ = ax[1]
c = MLC['fig3']['Catania (coal.+frag.)']
curve(a_, c['x'], c['y'], 1, 'Catania (coal.+frag.)')
t = hep('ins1762441/*/Figure3.csv')
for k, dx in (('fw', 0.0), ('pub', 5.0), ('c1', -5.0)):
    q, qe = [], []
    for cl in ('50-80', '20-50', '0-20'):
        A_, Ae = val(R[k], f'dN_Lc_3to6_{cl}_perc'); B_, Be = val(R[k], f'dN_D0_3to6_{cl}_perc')
        q.append(A_[0] / B_[0]); qe.append(A_[0] / B_[0] * np.hypot(Ae[0] / A_[0], Be[0] / B_[0]))
    model(a_, t[:, 0], np.array(q), np.array(qe), k, dx)
data(a_, t[:, 0], t[:, 1], t[:, 2], t[:, 4], 6, label='STAR, $3<p_T<6$ GeV/$c$')
a_.set(xlabel=r'$\langle N_{\rm part}\rangle$', ylabel=r'$\Lambda_c^+/D^0$', xlim=(0, 330), ylim=(0, 1.9))
a_.legend(loc='upper left', fontsize=6.2)
save(fig, 'fig_lc')

# =========================================================== Fig. 4  hadrochemistry (|y| < 1, per c quark)
SPEC = [('D0', 'D$^0$'), ('Dp', 'D$^+$'), ('Ds', 'D$_s^+$'), ('Lc', r'$\Lambda_c^+$')]


def fractions(r, cl):
    n = np.array([val(r, f'N_y1{s}_{cl}_perc')[0][0] for s, _ in SPEC])
    return n / n.sum()


STARF = np.array([41, 18, 15, 78], float)
STARFE = np.hypot([1, 1, 1, 13], [5, 3, 5, 28])
fig, a_ = plt.subplots(figsize=(W1, 2.6))
xs = np.arange(4)
groups = [('p+p', R['pp'], 'mb', '#c9c8c2'), ('Au+Au 40–80 %', R['fw'], '40-80', '#9ec5ef'),
          ('Au+Au 10–40 %', R['fw'], '10-40', '#5a9be3'), ('Au+Au 0–10 %', R['fw'], '0-10', '#1f5fb0')]
for i, (lab, r, cl, col) in enumerate(groups):
    a_.bar(xs + (i - 1.5) * 0.17, fractions(r, cl), 0.15, color=col, label=f'This work, {lab}', zorder=2)
sf = STARF / STARF.sum()
a_.errorbar(xs + 0.4, sf, sf * np.sqrt((STARFE / STARF) ** 2), fmt='o', color=INK, ms=3.6, lw=0.9, capsize=0,
            label='STAR 10–40 % (prelim.)', zorder=5)
a_.set(xticks=xs, xticklabels=[s for _, s in SPEC], ylabel='fraction ($|y|<1$, per $c$ quark)', ylim=(0, 0.85))
a_.legend(loc='upper left', fontsize=6.2)
save(fig, 'fig_hadrochem')

# =========================================================== Fig. 5  D0 v2 and NCQ
fig, ax = plt.subplots(1, 3, figsize=(W2, 2.6), gridspec_kw=dict(wspace=0.36))
for j, (cl, fn) in enumerate((('10-40', 3), ('0-80', 4))):
    a_ = ax[j]
    t = hep(f'ins1510298/*/Figure{fn}.csv')
    if cl == '0-80':
        nm = {'TAMU_cdiff': 'TAMU', 'SUBATECH': 'SUBATECH', 'hydro3D': '3D visc. hydro', 'PHSD': 'PHSD', 'LBT': 'LBT',
              'Duke': 'Duke'}
        for i, k in enumerate(nm):
            c = MV2['curves'][k]
            curve(a_, c['pt'], c['v2'], i, nm[k], lw=0.9)
    ed = np.array(R['fw']['edges']['PT_V2'])
    for k, dx in (('fw', 0.0), ('pub', 0.08)):
        v, e = val(R[k], f'v2SPS_D0_{cl}')
        model(a_, t[:, 0], v, e, k, dx)
    data(a_, t[:, 0], t[:, 1], t[:, 2], t[:, 4], 0.09, label='STAR D$^0$')
    a_.axhline(0, color=MUTED, lw=0.6)
    a_.set(xlabel='$p_T$ (GeV/$c$)', xlim=(0, 7), ylim=(-0.06, 0.27))
    tag(a_, f'{cl.replace("-", "–")} %', x=0.04)
h, l = [], []
for a_ in ax[:2]:
    hh, ll = a_.get_legend_handles_labels()
    for x_, y_ in zip(hh, ll):
        if y_ not in l:
            h.append(x_); l.append(y_)
fig.legend(h, l, loc='lower center', ncol=5, fontsize=6.4, bbox_to_anchor=(0.36, -0.24))
ax[0].set_ylabel('$v_2$ (D$^0$, $|y|<1$)')
a_ = ax[2]
xk = np.array(R['fw']['edges']['KET']); xc = 0.5 * (xk[1:] + xk[:-1])
for i, (s, lab) in enumerate((('pi', r'$\pi$'), ('K', 'K'), ('p', 'p'))):
    v, e = val(R['fw'], f'v2ncq_{s}_10-40')
    a_.errorbar(xc, v, e, fmt='o-', ms=2.6, lw=0.8, color=CM[i + 3], label=lab)
v, e = val(R['fw'], 'v2ncq_D0_10-40')
a_.errorbar(xc, v, e, fmt='s-', ms=3.8, lw=1.2, color=C_FW, label='D$^0$')
a_.set(xlabel=r'$(m_T-m_0)/n_q$ (GeV)', ylabel=r'$v_2/n_q$ (reaction plane)', xlim=(0, 1.6), ylim=(-0.005, 0.08))
tag(a_, 'This work, 10–40 %', x=0.04)
a_.legend(loc='lower right', fontsize=6)
save(fig, 'fig_v2')

# =========================================================== Fig. 6  D0 v1 and Delta v1
fig, ax = plt.subplots(1, 2, figsize=(W2, 2.8), gridspec_kw=dict(wspace=0.25))
yc = np.array([-0.6, -0.2, 0.2, 0.6])
a_ = ax[0]
curve(a_, MV1['v1_D0avg']['hydro_EM_Chatterjee_Bozek']['y'], MV1['v1_D0avg']['hydro_EM_Chatterjee_Bozek']['v1'], 0,
      'Hydro+EM (Chatterjee–Bożek)')
curve(a_, MV1['v1_D0avg']['AMPT_Nasim_Singha']['y'], MV1['v1_D0avg']['AMPT_Nasim_Singha']['v1'], 2, 'AMPT (Nasim–Singha)')
c = OLI['top']['eta_m 1.1']
curve(a_, c['y'], c['v1'], 3, 'Catania (Oliva et al.), $b$ = 9 fm')
for k, dx in (('fw', 0.0), ('pub', 0.03)):
    v, e = val(R[k], 'v1_D0avg_10-80')
    model(a_, yc, v, e, k, dx)
t = hep('ins1733225/*/Figure3a*.csv')
data(a_, t[:, 0], t[:, 1], t[:, 2], t[:, 4], 0.025, label='STAR D$^0$+$\\bar{\\rm D}^0$')
a_.axhline(0, color=MUTED, lw=0.6)
a_.set(xlabel='$y$', ylabel=r'$v_1$', xlim=(-1, 1), ylim=(-0.08, 0.11))
tag(a_, 'Au+Au 10–80 %, $p_T>1.5$ GeV/$c$', x=0.04, y=0.08)
a_.legend(loc='upper right', fontsize=6)
a_ = ax[1]
curve(a_, MV1['dv1_D0_minus_D0bar']['hydro_EM_Chatterjee_Bozek']['y'], MV1['dv1_D0_minus_D0bar']['hydro_EM_Chatterjee_Bozek']['dv1'],
      0, 'Hydro+EM (Chatterjee–Bożek)')
curve(a_, MV1['dv1_D0_minus_D0bar']['EM_Das_et_al']['y'], MV1['dv1_D0_minus_D0bar']['EM_Das_et_al']['dv1'], 4, 'EM only (Das et al.)')
for k, dx in (('fw', 0.0), ('pub', 0.03)):
    v0, e0 = val(R[k], 'v1_D0_10-80'); v1_, e1 = val(R[k], 'v1_D0bar_10-80')
    model(a_, yc, v0 - v1_, np.hypot(e0, e1), k, dx)
t = hep('ins1733225/*/Figure3b*.csv')
data(a_, t[:, 0], t[:, 1], t[:, 2], t[:, 4], 0.025, label=r'STAR D$^0-\bar{\rm D}^0$')
a_.axhline(0, color=MUTED, lw=0.6)
a_.set(xlabel='$y$', ylabel=r'$\Delta v_1$ (D$^0-\bar{\rm D}^0$)', xlim=(-1, 1), ylim=(-0.09, 0.09))
a_.legend(loc='upper right', fontsize=6)
save(fig, 'fig_v1')

# =========================================================== Fig. 7  Ds/D0
fig, ax = plt.subplots(1, 2, figsize=(W2, 2.6), gridspec_kw=dict(wspace=0.25))
t = hep('ins1843268/*/Figure2a(10-40%).csv')
a_ = ax[0]
for i, (k, c) in enumerate(MDS['fig3b'].items()):
    if '20-40' in k or 'PYTHIA' in k:
        curve(a_, c['x'], c['y'], 6 if 'PYTHIA' in k else 1, k.replace(', 20-40 %', ' (20–40 %)').replace('PYTHIA 8 p+p', 'PYTHIA 8 (p+p)'))
groups = [[0, 1], [2], [3], [4], [5]]
for k, pk, dx in (('fw', 'fw', 0.0), ('pub', 'pub', 0.08)):
    A_, Ae = val(R[k], 'dN_Ds_dsbins_10-40_perc'); B_, Be = val(R[k], 'dN_D0_dsbins_10-40_perc')
    a = np.array([A_[g].sum() for g in groups]); ae = np.array([np.sqrt((Ae[g] ** 2).sum()) for g in groups])
    b = np.array([B_[g].sum() for g in groups]); be = np.array([np.sqrt((Be[g] ** 2).sum()) for g in groups])
    model(a_, t[:, 0], a / b, a / b * np.hypot(ae / a, be / b), k, dx)
data(a_, t[:, 0], t[:, 1], t[:, 2], 0.5 * (np.abs(t[:, 4]) + np.abs(t[:, 5])), 0.1, label='STAR 10–40 %')
a_.set(xlabel='$p_T$ (GeV/$c$)', ylabel=r'$D_s^\pm/D^0$', xlim=(0, 8.5), ylim=(0, 0.9))
a_.legend(loc='upper right', fontsize=6)
a_ = ax[1]
t = hep('ins1843268/*/Figure4.csv')
c = MDS['fig4']['Tsinghua (seq. coal.)']
curve(a_, c['x'], c['y'], 1, 'Tsinghua (seq. coal.)')
for k, dx in (('fw', 0.0), ('pub', 5.0)):
    q, qe = [], []
    for cl in ('d60-80', 'd40-60', 'd20-40', 'd10-20', '0-10'):
        A_, Ae = val(R[k], f'dN_Ds_dsbins_{cl}_perc'); B_, Be = val(R[k], f'dN_D0_dsbins_{cl}_perc')
        a, ae = A_[1:5].sum(), np.sqrt((Ae[1:5] ** 2).sum()); b, be = B_[1:5].sum(), np.sqrt((Be[1:5] ** 2).sum())
        q.append(a / b); qe.append(a / b * np.hypot(ae / a, be / b))
    model(a_, t[:, 0], np.array(q), np.array(qe), k, dx)
data(a_, t[:, 0], t[:, 1], t[:, 2], 0.5 * (np.abs(t[:, 4]) + np.abs(t[:, 5])), 6, label='STAR, $1.5<p_T<5$ GeV/$c$')
a_.set(xlabel=r'$\langle N_{\rm part}\rangle$', ylabel=r'$D_s^\pm/D^0$', xlim=(0, 340), ylim=(0, 0.9))
a_.legend(loc='upper left', fontsize=6)
save(fig, 'fig_ds')
