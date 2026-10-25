#!/usr/bin/env python3
"""compare_bulk.py OUTDIR ARM.json [ARM.json ...]

Bulk (light-hadron) flow of each arm against the published 200 GeV Au+Au data, with
the model in the experiment's own kinematic windows (analyze_arm.py v3 keys):

  directed flow, slopes dv1/dy at midrapidity
    STAR PRL 108 202301 (HEPData ins1081744), 10-70 %:
        pi 0.15-0.75 GeV  -0.27 +- 0.01 +- 0.01 %     K0S       -0.17 +- 0.02 +- 0.04 %
        p  0.40-1.00 GeV  -0.15 +- 0.05 +- 0.08 %     pbar      -0.46 +- 0.06 +- 0.04 %
        plus a 10 % common (event-plane) uncertainty
        and the pi / p / pbar slope versus centrality (their Fig. 3)
    STAR PRL 123 162301 (D0 paper), 10-80 %, K+- pT > 0.2 GeV: -0.0030 +- 0.0001 +- 0.0002
  elliptic flow v2(pT), |y| < 1
    STAR PRC 77 054901 (ins777954) K0S v2{EP}, 10-40 % and 0-80 %
    PHENIX PRC 85 064914 (ins1093596) pi, K, p v2, 10-40 %, pT > 1 GeV

The model v1 uses the true reaction plane (STAR: spectator ZDC-SMD plane).  v2: with
BULK_V2=ep (2026-10-22) the model is measured like each experiment -- STAR K0S with the
eta-sub event plane (analyze_arm v5 v2EPS_*), PHENIX with the 1<|eta|<2.8 plane (v2EPF_*);
the default stays v2{RP} (v2_*), which is 25-30 % below v2{EP} at 10-40 % in AMPT.  Kaon v1: our K+- is compared with K0S in
10-70 % (same quark content, STAR's K+- has 6x the error).
Writes OUTDIR/BULK_FLOW.md, OUTDIR/fig_bulk_v1.png, OUTDIR/fig_bulk_v2.png.
"""
import json
import os
import sys

import numpy as np

HD = str(__import__('pathlib').Path(__file__).resolve().parents[1] / 'data' / 'hepdata')
OUT = sys.argv[1]
ARMS = sys.argv[2:]
os.makedirs(OUT, exist_ok=True)

V2M = os.environ.get('BULK_V2', 'rp')        # rp | ep | sp (sp: scalar product, analyze_arm v6)
EP = V2M != 'rp'
KS, KF = dict(rp=('v2_', 'v2_'), ep=('v2EPS_', 'v2EPF_'), sp=('v2SPS_', 'v2SPF_'))[V2M]   # STAR-, PHENIX-geometry keys
COMMON = 0.10                   # STAR first-order event-plane common uncertainty (relative)
# (value, stat, sys) in absolute units
V1_REF = {
    'pi_10-70': (-0.0027, 0.0001, 0.0001, 'STAR PRL 108 202301, pi 0.15-0.75 GeV'),
    'K_10-70': (-0.0017, 0.0002, 0.0004, 'STAR PRL 108 202301, K0S'),
    'p_10-70': (-0.0015, 0.0005, 0.0008, 'STAR PRL 108 202301, p 0.4-1.0 GeV'),
    'pbar_10-70': (-0.0046, 0.0006, 0.0004, 'STAR PRL 108 202301, pbar 0.4-1.0 GeV'),
    'K_10-80': (-0.0030, 0.0001, 0.0002, 'STAR PRL 123 162301, K+- pT>0.2 GeV'),
}
D0_REF = (-0.080, 0.017, 0.016)        # STAR PRL 123 162301, D0+D0bar, pT > 1.5 GeV, 10-80 %
CENT = ['c0-5', 'c5-10', 'c10-20', 'c20-30', 'c30-40', 'c40-50', 'c50-60', 'c60-70', 'c70-80']


def hep_blocks(fn):
    """HEPData csv -> list of (meta lines, header, rows)"""
    blocks, meta, hdr, rows = [], [], None, []
    for line in open(fn):
        line = line.rstrip('\n')
        if line.startswith('#:'):
            if hdr is not None:
                blocks.append((meta, hdr, rows)); meta, hdr, rows = [], None, []
            meta.append(line)
        elif not line.strip():
            if hdr is not None:
                blocks.append((meta, hdr, rows)); meta, hdr, rows = [], None, []
        elif hdr is None:
            hdr = line
        else:
            rows.append([float(x) if x.strip().strip('-') else np.nan for x in line.split(',')])
    if hdr is not None:
        blocks.append((meta, hdr, rows))
    return blocks


def star_v1_cent():
    out = {}
    for meta, hdr, rows in hep_blocks(f'{HD}/ins1081744/HEPData-ins1081744-v1-csv/Figure3.csv'):
        sp = [m for m in meta if 'Particle' in m][0]
        key = 'pi' if 'pi' in sp else ('pbar' if 'bar' in sp else 'p')
        r = np.array(rows)
        out[key] = (r[:, 0], r[:, 1] / 100, r[:, 2] / 100)      # centrality, slope, stat
    return out


def k0s_v2(cls):
    for meta, hdr, rows in hep_blocks(f'{HD}/ins777954/HEPData-ins777954-v1-csv/Figure7(K_{{S}}^{{0}}).csv'):
        if any(m.endswith('$' + cls + '\\%$') for m in meta):
            r = np.array(rows)
            r = r[np.all(np.isfinite(r), 1)]
            return r[:, 0], r[:, 1], np.hypot(r[:, 2], r[:, 4])
    raise KeyError(cls)


def phenix_v2(i):
    fn = f'{HD}/ins1093596/HEPData-ins1093596-v1-csv/Figure10b-{i}.csv'
    for meta, hdr, rows in hep_blocks(fn):
        if 'v_2' in hdr:
            r = np.array(rows)
            return r[:, 0], r[:, 3], np.hypot(r[:, 4], r[:, 6])
    raise KeyError(fn)


def get(r, k):
    v = r.get(k)
    return (np.nan, np.nan) if v is None else (v[0], v[1])


def v2_chi2(r, key, data, ptmin=0.0, ptmax=3.0):
    e = np.array(r['edges']['PT_L'])
    c = 0.5 * (e[1:] + e[:-1])
    if key not in r:
        return np.nan, 0
    mv, me = np.array(r[key]['val']), np.array(r[key]['err'])
    x, y, s = data
    ok = (c >= max(ptmin, x.min())) & (c <= min(ptmax, x.max())) & np.isfinite(mv) & (me > 0)
    d = np.interp(c[ok], x, y)
    ds = np.interp(c[ok], x, s)
    return float(np.sum((mv[ok] - d) ** 2 / (me[ok] ** 2 + ds ** 2))), int(ok.sum())


R = []
for f in ARMS:
    r = json.load(open(f))
    r['_name'] = os.path.basename(f).replace('.json', '')
    R.append(r)
SC = star_v1_cent()

lines = ['# Bulk flow against published 200 GeV Au+Au data', '',
         'Model: reaction-plane v1/v2 in STAR windows (analyze_arm.py v3). Errors: model jackknife; '
         'data stat (+) sys (+) 10 % common for v1.', '',
         '## Directed-flow slopes dv1/dy', '',
         '| arm | π 10–70 | K 10–70 | p 10–70 | p̄ 10–70 | K 10–80 | χ²(5) | π vs cent χ²(9) | partons 10–70 | D0 10–80 |',
         '|---|---|---|---|---|---|---|---|---|---|']
ref = '| **STAR** | ' + ' | '.join(f'{V1_REF[k][0]:+.4f}' for k in ('pi_10-70', 'K_10-70', 'p_10-70', 'pbar_10-70', 'K_10-80')) \
      + f' | — | — | — | {D0_REF[0]:+.3f} |'
lines.append(ref)
summary = {}
for r in R:
    cells, chi = [], 0.0
    for k in ('pi_10-70', 'K_10-70', 'p_10-70', 'pbar_10-70', 'K_10-80'):
        m, me = get(r, 'slopeS_' + k)
        d, ds, dy, _ = V1_REF[k]
        sd2 = ds ** 2 + dy ** 2 + (COMMON * d) ** 2
        chi += (m - d) ** 2 / (me ** 2 + sd2)
        cells.append(f'{m:+.4f} ± {me:.4f}')
    cc, x, y, s = 0.0, *SC['pi']
    for i, k in enumerate(CENT):
        m, me = get(r, f'slopeS_pi_{k}')
        cc += (m - y[i]) ** 2 / (me ** 2 + s[i] ** 2 + (COMMON * y[i]) ** 2)
    lq = get(r, 'slopeS_lq_10-70')
    dd = get(r, 'slope_D0avg_10-80')
    summary[r['_name']] = dict(chi_v1=chi, chi_cent=cc)
    lines.append(f"| {r['_name']} | " + ' | '.join(cells) + f' | {chi:.1f} | {cc:.1f} | '
                 f'{lq[0]:+.4f} ± {lq[1]:.4f} | {dd[0]:+.3f} ± {dd[1]:.3f} |')

lines += ['', '## Symmetry check: rapidity-even v1 relative to the reaction plane (10–70 %)', '',
          'Must vanish for Au+Au with the true reaction plane; STAR (spectator plane) and the slopes see only the odd part.', '',
          '| arm | π | K | p | p̄ |', '|---|---|---|---|---|']
for r in R:
    cells = []
    for sp in ('pi', 'K', 'p', 'pbar'):
        k = f'v1s_{sp}_10-70'
        if k not in r:
            cells.append('—'); continue
        v, e = np.array(r[k]['val']), np.array(r[k]['err'])
        ev, ee = 0.5 * (v + v[::-1])[:3], 0.5 * np.hypot(e, e[::-1])[:3]
        w = 1 / ee ** 2
        cells.append(f'{np.sum(w * ev) / np.sum(w):+.4f} ± {1 / np.sqrt(np.sum(w)):.4f}')
    lines.append(f"| {r['_name']} | " + ' | '.join(cells) + ' |')

lines += ['', '## Elliptic flow v2(pT), |y| < 1 (χ² against interpolated data, number of points in brackets; model estimator: ' + dict(rp='true reaction plane', ep='event plane as in each experiment', sp='scalar product in each experiment\'s geometry')[V2M] + ')', '',
          '| arm | K vs STAR K0S 10–40 | K vs STAR K0S 0–80 | π vs PHENIX 10–40 (pT>1) | K vs PHENIX 10–40 | p vs PHENIX 10–40 |',
          '|---|---|---|---|---|---|']
D = dict(k1040=k0s_v2('10-40'), k080=k0s_v2('0-80'), ppi=phenix_v2(1), pk=phenix_v2(2), pp=phenix_v2(3))
for r in R:
    c = [v2_chi2(r, KS + 'K_10-40', D['k1040'], 0.3, 3.0), v2_chi2(r, KS + 'K_0-80', D['k080'], 0.3, 3.0),
         v2_chi2(r, KF + 'pi_10-40', D['ppi'], 1.0, 3.0), v2_chi2(r, KF + 'K_10-40', D['pk'], 1.0, 3.0),
         v2_chi2(r, KF + 'p_10-40', D['pp'], 1.0, 3.0)]
    summary[r['_name']]['chi_v2'] = sum(x[0] for x in c)
    summary[r['_name']]['n_v2'] = sum(x[1] for x in c)
    lines.append(f"| {r['_name']} | " + ' | '.join(f'{a:.1f} ({n})' for a, n in c) + ' |')

lines += ['', '## Bulk chemistry at |y| < 0.5 (analysis v3 yields)', '',
          'STAR K−/π− (PRC 79 034909, HEPData ins793126 Fig. 32): 0.152 ± 0.019 (N_part 299–351), 0.142 ± 0.014 (N_part 78).', '',
          '| arm | Nch(\\|η\\|<0.5) 0–10 % | K−/π− 0–10 % | K−/π− 40–80 % | p̄/p 0–10 % | p/π+ 0–10 % |', '|---|---|---|---|---|---|']
for r in R:
    def y_(nm, k):
        d = r.get(f'yld_{nm}_{k}')
        return d['val'][0] if d else np.nan
    n = r.get('nch5_0-10', {}).get('val', [np.nan])[0]
    lines.append(f"| {r['_name']} | {n:.0f} | {y_('Km', '0-10') / y_('pim', '0-10'):.3f} | {y_('Km', '40-80') / y_('pim', '40-80'):.3f} | "
                 f"{y_('pbar', '0-10') / y_('p', '0-10'):.3f} | {y_('p', '0-10') / y_('pip', '0-10'):.3f} |")

lines += ['', '## Parton-phase end (src_v9 Z record)', '', '| arm | events with Z | partons frozen only by the HQFRZG guard |', '|---|---|---|']
for r in R:
    c = r['counts']
    if c.get('Z_events', 0):
        lines.append(f"| {r['_name']} | {c['Z_events']} | {c['Z_nfrz0']} of {c['Z_nfrzt']} "
                     f"({100 * c['Z_nfrz0'] / max(c['Z_nfrzt'], 1):.3f} %) |")
open(f'{OUT}/BULK_FLOW.md', 'w').write('\n'.join(lines) + '\n')
json.dump(summary, open(f'{OUT}/bulk_summary.json', 'w'), indent=1)
print('\n'.join(lines))

# ------------------------------------------------------------------ figures
import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

cols = plt.rcParams['axes.prop_cycle'].by_key()['color']
fig, ax = plt.subplots(2, 3, figsize=(13, 7.5))
ys = np.array(R[0]['edges']['Y_S']); yc = 0.5 * (ys[1:] + ys[:-1])
f1 = {}
for meta, hdr, rows in hep_blocks(f'{HD}/ins1081744/HEPData-ins1081744-v1-csv/Figure1(left).csv') + \
        hep_blocks(f'{HD}/ins1081744/HEPData-ins1081744-v1-csv/Figure1(right).csv'):
    sp = [m for m in meta if 'Particle' in m][0].split(',')[-1]
    f1[sp] = np.array(rows)
panels = [('pi', '$\\pi^{\\pm}$', '$\\pi^{\\pm}$'), ('K', 'K (model K$^\\pm$, data K$^0_S$)', '$K_{S}^{0}$'),
          ('p', 'p', '$p$'), ('pbar', '$\\bar p$', '$\\bar{p}$')]
for ia, (sp, title, dkey) in enumerate(panels):
    a = ax.flat[ia]
    d = f1[dkey]
    a.errorbar(d[:, 0], d[:, 1] / 100, d[:, 2] / 100, fmt='ks', label='STAR 10–70 %', zorder=5)
    for ir, r in enumerate(R):
        k = f'v1s_{sp}_10-70'
        if k in r:   # rapidity-odd part, as STAR's spectator-plane v1 measures
            v, e = np.array(r[k]['val']), np.array(r[k]['err'])
            a.errorbar(yc + 0.01 * ir, 0.5 * (v - v[::-1]), 0.5 * np.hypot(e, e[::-1]), fmt='o-', ms=3,
                       color=cols[ir % 10], label=r['_name'])
    a.axhline(0, color='0.6', lw=0.8); a.set_title(title); a.set_xlabel('y'); a.set_ylabel('$v_1^{odd}$')
a = ax.flat[4]
x, y, s = SC['pi']
a.errorbar(x, y, s, fmt='ks', label='STAR π')
cx = np.array([2.5, 7.5, 15, 25, 35, 45, 55, 65, 75])
for ir, r in enumerate(R):
    v = [get(r, f'slopeS_pi_{k}') for k in CENT]
    a.errorbar(cx + 0.5 * ir, [t[0] for t in v], [t[1] for t in v], fmt='o-', ms=3, color=cols[ir % 10])
a.axhline(0, color='0.6', lw=0.8); a.set_xlabel('centrality (%)'); a.set_ylabel('d$v_1$/dy  π'); a.set_title('pion slope vs centrality')
a = ax.flat[5]
yq = np.array(R[0]['edges']['Y_Q']); qc = 0.5 * (yq[1:] + yq[:-1])
for ir, r in enumerate(R):
    if 'v1_lq_all' in r:
        a.errorbar(qc, r['v1_lq_all']['val'], r['v1_lq_all']['err'], fmt='o-', ms=3, color=cols[ir % 10], label=r['_name'])
a.axhline(0, color='0.6', lw=0.8); a.set_xlim(-3.5, 3.5); a.set_ylim(-0.03, 0.03)
a.set_xlabel('y'); a.set_ylabel('$v_1$ light partons'); a.set_title('u+d at parton freeze-out (0–80 %)')
ax.flat[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(f'{OUT}/fig_bulk_v1.png', dpi=130)

fig, ax = plt.subplots(1, 3, figsize=(13, 4))
e = np.array(R[0]['edges']['PT_L']); c = 0.5 * (e[1:] + e[:-1])
for a, (key, dat, lab) in zip(ax, ((KS + 'K_10-40', D['k1040'], 'K, 10–40 % (STAR K$^0_S$)'),
                                   (KF + 'pi_10-40', D['ppi'], 'π, 10–40 % (PHENIX)'),
                                   (KF + 'p_10-40', D['pp'], 'p, 10–40 % (PHENIX)'))):
    a.errorbar(dat[0], dat[1], dat[2], fmt='ks', ms=3, label='data')
    for ir, r in enumerate(R):
        if key in r:
            a.errorbar(c, r[key]['val'], r[key]['err'], fmt='o-', ms=3, color=cols[ir % 10], label=r['_name'])
    a.set_xlim(0, 3.2); a.set_ylim(0, 0.3); a.set_xlabel('$p_T$ (GeV)'); a.set_ylabel('$v_2$'); a.set_title(lab)
ax[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(f'{OUT}/fig_bulk_v2.png', dpi=130)
print('figures written to', OUT)
