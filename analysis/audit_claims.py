#!/usr/bin/env python3
"""audit_claims.py OUT.md : apples-to-apples comparison of this work, public AMPT and every published model curve
we hold, against STAR, for the claims of the paper (2026-10-25, docs/CLAIMS_AUDIT.md).

Independent of compare_final.py / pilot_v14.py: own HEPData reader, own chi2.  One convention for every entity:
  chi2 = sum (m - d)^2 / sd^2,  sd^2 = stat^2 + sys^2 (sys symmetrised as the mean of |+| and |-|);
  D0 v2 0-80 % adds STAR's one-sided non-flow^2 where m < d (STAR PRL 118 Table I).
  Published curves carry no model error -> plain chi2.  Our transport output carries a statistical error ->
  chi2_inf = sum [(m - d)^2 - sm^2] / sd^2 (its expectation equals the chi2 of the infinite-statistics model),
  so both columns estimate the same quantity.
  Curves are evaluated two ways: at the STAR abscissa ('pt') and averaged over the STAR bin ('bin'); our output is
  already a bin average.  A curve that does not cover a bin is flagged (*).
D0 dv1/dy: straight line through the origin over the four STAR y bins (STAR's own fit), the same STAR-stat
weights for every entity; curves evaluated at y = +-0.2, +-0.6 and averaged over the 0.4-wide bins.
"""
import csv
import glob
import json
import sys

import numpy as np

D = str(__import__('pathlib').Path(__file__).resolve().parents[1])   # repository root
OUT = sys.argv[1]


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


def J(f):
    return json.load(open(f'{D}/{f}'))


OURS, OPP, PUB, PPP = (J('observation/results/json/dm25_main.json'), J('observation/results/json/dm25_ppmain.json'),
                       J('observation/results/json/dm25_pub.json'), J('observation/results/json/dm25_pppub.json'))
FAN, MV2, MV1 = J('data/models/fan2026_improved_ampt.json'), J('data/models/star_d0v2_fig4_models.json'), J('data/models/star_d0v1_fig3_models.json')
MLC, MRAA, OLI = J('data/models/star_lc_models.json'), J('data/models/star_raa_models.json'), J('data/models/oliva2021_d0v1.json')
L = ['# Claims audit: every entity on the same footing (ana/audit_claims.py)', '', __doc__.split('\n', 1)[1].strip(), '']


def val(r, k):
    return np.array(r[k]['val'], float), np.array(r[k]['err'], float)


def c2_curve(m, d, sd):
    return float(np.sum((m - d) ** 2 / sd ** 2))


def c2_inf(m, sm, d, sd):
    return float(np.sum(((m - d) ** 2 - sm ** 2) / sd ** 2)), float(np.sqrt(np.sum((4 * sm ** 2 * (m - d) ** 2 + 2 * sm ** 4) / sd ** 4)))


def curve_at(x, y, xs, edges=None):
    x, y = np.asarray(x, float), np.asarray(y, float)
    o = np.argsort(x); x, y = x[o], y[o]
    if edges is None:
        return np.interp(xs, x, y), np.zeros(len(xs), bool)
    m = np.array([np.mean(np.interp(np.linspace(a, b, 401), x, y)) for a, b in edges])
    cover = np.array([(a >= x[0] - 1e-6) and (b <= x[-1] + 1e-6) for a, b in edges])
    return m, ~cover


def row(obs, name, c, n, note=''):
    return f'| {obs} | {name} | {c} | {n} | {note} |'


# ----------------------------------------------------------------- data tables
V2 = {cl: hep(f'ins1510298/*/Figure{i}.csv') for cl, i in (('10-40', 3), ('0-80', 4))}
V2EDGES = list(zip([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0], [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 8.0]))
RAA = {cl: hep('ins1711377/*/RAA.csv', i) for i, cl in enumerate(('0-10', '10-40', '40-80'))}
LCPT, LCC = hep('ins1762441/*/Figure2.csv'), hep('ins1762441/*/Figure3.csv')
LCEDGES = [(2.5, 3.5), (3.5, 5.0), (5.0, 8.0)]
V1 = hep('ins1733225/*/Figure3a*.csv')
T = {}           # (observable) -> list of (entity, chi2 or (chi2inf, err), ndf, note)

# ----------------------------------------------------------------- D0 v2
for cl, t in V2.items():
    d, st, sy, nf = t[:, 1], t[:, 2], 0.5 * (np.abs(t[:, 4]) + np.abs(t[:, 5])), np.abs(t[:, 7])
    key = f'D0 v2 {cl} %'
    for nm, r in (('this work', OURS), ('public AMPT', PUB)):
        m, sm = val(r, f'v2SPS_D0_{cl}')
        sd = np.sqrt(st ** 2 + sy ** 2 + (np.where(m < d, nf ** 2, 0) if cl == '0-80' else 0))
        T.setdefault(key, []).append((nm, c2_inf(m, sm, d, sd), len(d), 'v2{SP}, STAR eta-sub geometry'))
        if cl == '10-40':   # sensitivity: STAR's one-sided non-flow also in 10-40 %
            sd2 = np.sqrt(st ** 2 + sy ** 2 + np.where(m < d, nf ** 2, 0))
            T.setdefault(key + ' (non-flow incl.)', []).append((nm, c2_inf(m, sm, d, sd2), len(d), 'sensitivity'))
    if cl == '0-80':
        for mk, c in MV2['curves'].items():
            for how, ed in (('pt', None), ('bin', V2EDGES)):
                m, unc = curve_at(c['pt'], c['v2'], t[:, 0], ed)
                ok = ~unc if ed is not None else (t[:, 0] >= min(c['pt'])) & (t[:, 0] <= max(c['pt']))
                sd = np.sqrt(st ** 2 + sy ** 2 + np.where(m < d, nf ** 2, 0))
                T.setdefault(key, []).append((f'{mk} [{how}]', c2_curve(m[ok], d[ok], sd[ok]), int(ok.sum()),
                                              f'STAR Table I: {MV2["chi2_star"][mk][0]} for {MV2["chi2_star"][mk][1]} d.o.f.'))

# ----------------------------------------------------------------- R_AA (8 STAR points, 0.9 < pT < 8)
for cl, t in RAA.items():
    sel = (t[:, 3] > 0.9) & (t[:, 3] < 8)
    t = t[sel]
    d, sd = t[:, 4], np.hypot(t[:, 5], 0.5 * (np.abs(t[:, 7]) + np.abs(t[:, 8])))
    ed = list(zip(t[:, 1], t[:, 2]))
    key = f'D0 R_AA {cl} %'
    for nm, r, pp in (('this work', OURS, OPP), ('public AMPT', PUB, PPP)):
        a, ae = val(r, f'dN_D0_{cl}_perc'); p, pe = val(pp, 'dN_D0_mb_perc')
        # model bins matched to the STAR points by bin edges (not by index): PT_RAA == STAR's binning
        e = np.array(r['edges']['PT_RAA']); ib = [int(np.argmin(np.abs(e[:-1] - lo))) for lo in t[:, 1]]
        assert np.allclose(e[ib], t[:, 1], atol=2e-3) and np.allclose(e[np.array(ib) + 1], t[:, 2], atol=2e-3)
        m, sm = (a / p)[ib], (a / p * np.hypot(ae / a, pe / p))[ib]
        T.setdefault(key, []).append((nm, c2_inf(m, sm, d, sd), len(d), 'per-c-quark yields AA/pp'))
    curves = {'Fan et al. 2026 (improved AMPT)': (FAN['RAA_D0'][cl]['pt'], FAN['RAA_D0'][cl]['raa'])}
    if cl == '0-10':
        for k, c in MRAA['RAA_D0_0_10'].items():
            curves[k] = (c['pt'], c['raa'])
    for k, (x, y) in curves.items():
        for how, e in (('pt', None), ('bin', ed)):
            m, unc = curve_at(x, y, t[:, 3], e)
            T.setdefault(key, []).append((f'{k} [{how}]', c2_curve(m, d, sd), len(d), '*' * int(unc.any())))

# ----------------------------------------------------------------- Lambda_c / D0
d, sd = LCPT[:, 1], np.hypot(LCPT[:, 2], LCPT[:, 4])
for nm, r in (('this work', OURS), ('public AMPT', PUB)):
    a, ae = val(r, 'dN_Lc_lcbins_10-80_perc'); b, be = val(r, 'dN_D0_lcbins_10-80_perc')
    T.setdefault('Λc/D0 vs pT (10-80 %)', []).append((nm, c2_inf(a / b, a / b * np.hypot(ae / a, be / b), d, sd), 3, ''))
curves = {'Fan et al. 2026 (improved AMPT)': (FAN['LcD0']['AuAu_10-80']['pt'], FAN['LcD0']['AuAu_10-80']['ratio'])}
curves.update({k: (c['x'], c['y']) for k, c in MLC['fig2b'].items()})
for k, (x, y) in curves.items():
    for how, e in (('pt', None), ('bin', LCEDGES)):
        m, unc = curve_at(x, y, LCPT[:, 0], e)
        T['Λc/D0 vs pT (10-80 %)'].append((f'{k} [{how}]', c2_curve(m, d, sd), 3, '*' * int(unc.any())))
d, sd = LCC[:, 1], np.hypot(LCC[:, 2], LCC[:, 4])
for nm, r in (('this work', OURS), ('public AMPT', PUB)):
    q, qe = [], []
    for cl in ('50-80', '20-50', '0-20'):
        a, ae = val(r, f'dN_Lc_3to6_{cl}_perc'); b, be = val(r, f'dN_D0_3to6_{cl}_perc')
        q.append(a[0] / b[0]); qe.append(a[0] / b[0] * np.hypot(ae[0] / a[0], be[0] / b[0]))
    T.setdefault('Λc/D0 vs centrality (3-6 GeV)', []).append((nm, c2_inf(np.array(q), np.array(qe), d, sd), 3, ''))
for k in ('Catania (coal.)', 'Catania (coal.+frag.)'):
    c = MLC['fig3'][k]
    m, _ = curve_at(c['x'], c['y'], LCC[:, 0])
    T['Λc/D0 vs centrality (3-6 GeV)'].append((k, c2_curve(m, d, sd), 3, 'curve drawn at the 3 STAR Npart'))

# ----------------------------------------------------------------- Ds / D0 (not in the 46-point scoreboard)
MDS = J('data/models/star_ds_models.json')
DSG = {'0-10': [[1, 2], [3], [4], [5]], '10-40': [[0, 1], [2], [3], [4], [5]], '40-80': [[1, 2], [3], [4], [5]]}


def ds_ratio(r, cl, groups):
    a, ae = val(r, f'dN_Ds_dsbins_{cl}_perc'); b, be = val(r, f'dN_D0_dsbins_{cl}_perc')
    A_ = np.array([a[g].sum() for g in groups]); Ae = np.array([np.sqrt((ae[g] ** 2).sum()) for g in groups])
    B_ = np.array([b[g].sum() for g in groups]); Be = np.array([np.sqrt((be[g] ** 2).sum()) for g in groups])
    return A_ / B_, A_ / B_ * np.hypot(Ae / A_, Be / B_)


for cl, groups in DSG.items():
    t = hep(f'ins1843268/*/Figure2a({cl}%).csv')
    d, sd = t[:, 1], np.hypot(t[:, 2], 0.5 * (np.abs(t[:, 4]) + np.abs(t[:, 5])))
    key = f'Ds/D0 vs pT ({cl} %)'
    for nm, r in (('this work', OURS), ('public AMPT', PUB)):
        q, qe = ds_ratio(r, cl, groups)
        T.setdefault(key, []).append((nm, c2_inf(q, qe, d, sd), len(d), 'error ignores Ds-D0 correlation (upper bound)'))
    mods = MDS['fig3a_0_10'] if cl == '0-10' else {k: c for k, c in MDS['fig3b'].items()
                                                    if ('40-80' in k) == (cl == '40-80') and 'PYTHIA' not in k}
    for k, c in mods.items():
        if 'PYTHIA' in k:
            continue
        m, _ = curve_at(c['x'], c['y'], t[:, 0])
        T[key].append((k + (' [20-40 % curve]' if cl == '10-40' else ''), c2_curve(m, d, sd), len(d), 'curve'))
t = hep('ins1843268/*/Figure4.csv')
d, sd = t[:, 1], np.hypot(t[:, 2], 0.5 * (np.abs(t[:, 4]) + np.abs(t[:, 5])))
for nm, r in (('this work', OURS), ('public AMPT', PUB)):
    qq = [ds_ratio(r, cl, [[1, 2, 3, 4]]) for cl in ('d60-80', 'd40-60', 'd20-40', 'd10-20', '0-10')]
    T.setdefault('Ds/D0 vs centrality (1.5-5 GeV)', []).append(
        (nm, c2_inf(np.array([x[0][0] for x in qq]), np.array([x[1][0] for x in qq]), d, sd), len(d), ''))
c = MDS['fig4']['Tsinghua (seq. coal.)']
T['Ds/D0 vs centrality (1.5-5 GeV)'].append(('Tsinghua (seq. coal.)', c2_curve(curve_at(c['x'], c['y'], t[:, 0])[0], d, sd), len(d), 'curve'))

# ----------------------------------------------------------------- sums: the 46-point scoreboard and with Ds/D0
SB = ['D0 v2 10-40 %', 'D0 v2 0-80 %', 'D0 R_AA 0-10 %', 'D0 R_AA 10-40 %', 'D0 R_AA 40-80 %', 'Λc/D0 vs pT (10-80 %)',
      'Λc/D0 vs centrality (3-6 GeV)']
DS = ['Ds/D0 vs pT (0-10 %)', 'Ds/D0 vs pT (10-40 %)', 'Ds/D0 vs pT (40-80 %)', 'Ds/D0 vs centrality (1.5-5 GeV)']
sums = []
for nm in ('this work', 'public AMPT'):
    def tot(keys):
        rows_ = [next(x for x in T[k] if x[0] == nm) for k in keys]
        return sum(x[1][0] for x in rows_), float(np.sqrt(sum(x[1][1] ** 2 for x in rows_))), sum(x[2] for x in rows_)
    sums.append((nm, tot(SB), tot(DS), tot(SB + DS)))
L += ['## Sums', '', '| entity | 46-point scoreboard | Ds/D0 terms | all |', '|---|---|---|---|']
for nm, a, b, c in sums:
    L.append(f'| {nm} | {a[0]:.1f} ± {a[1]:.1f} ({a[2]}) | {b[0]:.1f} ± {b[1]:.1f} ({b[2]}) | {c[0]:.1f} ± {c[1]:.1f} ({c[2]}) |')
L.append('')

# ----------------------------------------------------------------- tables
for obs, rows in T.items():
    L += [f'## {obs}', '', '| entity | χ² (ours/public: χ²∞ ± err) | ndf | note |', '|---|---|---|---|']
    for nm, c, n, note in rows:
        cs = f'{c[0]:.1f} ± {c[1]:.1f}' if isinstance(c, tuple) else f'{c:.1f}'
        L.append(f'| {nm} | {cs} | {n} | {note} |')
    L.append('')

# ----------------------------------------------------------------- D0 dv1/dy
yc = np.array([-0.6, -0.2, 0.2, 0.6]); w = 1 / V1[:, 2] ** 2
edges = [(-0.8, -0.4), (-0.4, 0.0), (0.0, 0.4), (0.4, 0.8)]


def slope(v, ww=w):
    return float(np.sum(ww * yc * v) / np.sum(ww * yc * yc))


k_star = slope(V1[:, 1]); e_star = float(np.sqrt(1 / np.sum(w * yc * yc)))
sd_star = float(np.hypot(0.017, 0.016))
L += ['## D0+D0bar dv1/dy, 10-80 %, pT > 1.5 GeV/c', '',
      f'STAR points refitted through the origin with stat weights: {k_star:+.4f} ± {e_star:.4f} '
      f'(published −0.080 ± 0.017 stat; total with sys {sd_star:.3f}).', '',
      '| entity | slope [at y] | slope [bin avg] | (STAR − model)/σ | note |', '|---|---|---|---|---|']
ent = []
for nm, r in (('this work', OURS), ('public AMPT', PUB)):
    v, e = val(r, 'v1_D0avg_10-80')
    k_own = r['slope_D0avg_10-80']
    ent.append((nm, slope(v), slope(v), k_own[1], f'production slope (model weights) {k_own[0]:+.4f} ± {k_own[1]:.4f}'))
for mk, lab in (('hydro_EM_Chatterjee_Bozek', 'hydro+EM (Chatterjee–Bożek), STAR Fig. 3'),
                ('AMPT_Nasim_Singha', 'AMPT (Nasim–Singha), STAR Fig. 3')):
    c = MV1['v1_D0avg'][mk]
    ent.append((lab, slope(curve_at(c['y'], c['v1'], yc)[0]), slope(curve_at(c['y'], c['v1'], yc, edges)[0]), 0.0, 'curve'))
for pan, k in (('top', 'eta_m 1.1'), ('top', 'eta_m 1.2'), ('bottom', 'charm tilted like the bulk')):
    c = OLI[pan][k]
    ent.append((f'Catania (Oliva et al. 2021), {k}', slope(curve_at(c['y'], c['v1'], yc)[0]),
                slope(curve_at(c['y'], c['v1'], yc, edges)[0]), 0.0, 'curve, b = 9 fm'))
for nm, a, b, e, note in ent:
    L.append(f'| {nm} | {a:+.4f} | {b:+.4f} | {(a - (-0.080)) / np.hypot(sd_star, e):+.1f} | {note} |')
L.append('')
open(OUT, 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
