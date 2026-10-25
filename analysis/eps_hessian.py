#!/usr/bin/env python3
"""eps_hessian.py OUT.md PP.json SET01.json[.gz] ... SET31.json[.gz] : EPS09s NLO error band of the paper observables.

The 31 inputs are ONE sample analysed 31 times, each with the per-event weight of EPS09s NLO error set K
(analysis/analyze_arm.py --shw tab:K, weights from analysis/eps09s/eps09s_errw.py; README Sec. 2.3).  Set 1 is the
central set (= the in-code weight of the production).  For each observable X the EPS09 Hessian prescription
(Eskola, Paukkunen, Salgado, JHEP 04 (2009) 065, Eq. 5) gives the 90 % CL band
  dX+ = sqrt(sum_k max(X_{2k} - X_1, X_{2k+1} - X_1, 0)^2),  dX- = sqrt(sum_k max(X_1 - X_{2k}, X_1 - X_{2k+1}, 0)^2),
k = 1..15 (sets 2..31); the last column names the eigen-direction that dominates.  Observables and chi2 conventions
are those of analysis/audit_claims.py (STAR tables, model at the data abscissae, D0 v2 from the scalar product in
STAR's eta-sub geometry, STAR's one-sided non-flow in 0-80 %).  The p+p reference is unweighted, so it is the same
for every set."""
import csv
import glob
import gzip
import json
import os
import sys
from pathlib import Path

import numpy as np

np.seterr(all='ignore')                                                  # tiny test samples give empty bins
D = Path(__file__).resolve().parents[1]                                  # repository root
OUT, PP, SETS = sys.argv[1], sys.argv[2], sys.argv[3:]
assert len(SETS) == 31, f'need the 31 EPS09s sets, got {len(SETS)}'


def load(fn):
    return json.load(gzip.open(fn, 'rt') if fn.endswith('.gz') else open(fn))


def hd(pattern, i=0):
    """block i of a HEPData CSV table"""
    out, cur = [], []
    for row in csv.reader(open(glob.glob(f'{D}/data/hepdata/{pattern}')[0])):
        try:
            if not row or row[0].startswith('#'):
                raise ValueError
            cur.append([float(x.replace('−', '-')) for x in row])
        except ValueError:
            if cur:
                out.append(np.array(cur))
            cur = []
    if cur:
        out.append(np.array(cur))
    return out[i]


S_V2 = {'10-40': hd('ins1510298/*/Figure3.csv'), '0-80': hd('ins1510298/*/Figure4.csv')}
S_RAA = {k: hd('ins1711377/*/RAA.csv', i) for i, k in enumerate(('0-10', '10-40', '40-80'))}
S_LCPT = hd('ins1762441/*/Figure2.csv')
S_LCC = hd('ins1762441/*/Figure3.csv')
RPP = load(PP)


def g(r, k):
    return np.array(r[k]['val'], float), np.array(r[k]['err'], float)


def centres(r, key):
    e = np.array(r['edges'][key])
    return 0.5 * (e[1:] + e[:-1])


def at(xm, vm, em, xd):
    """model value at the data abscissae (nearest bin; the analysis bins are matched to the data)"""
    i = np.array([np.argmin(np.abs(xm - x)) for x in xd])
    return vm[i], em[i]


def ratio(a, ae, b, be):
    with np.errstate(invalid='ignore', divide='ignore'):
        q = a / b
        return q, q * np.hypot(ae / a, be / b)


def odd(r, key, lo=0.5, hi=2.5):
    yq = np.array(r['edges']['Y_Q']); c = 0.5 * (yq[1:] + yq[:-1])
    v, e = g(r, key)
    s = (np.abs(c) > lo) & (np.abs(c) < hi) & np.isfinite(v) & (e > 0)
    w = 1 / e[s] ** 2
    return np.sum(w * np.sign(c[s]) * v[s]) / np.sum(w)


def c2inf(m, sm, d, sd2):
    ok = np.isfinite(m) & np.isfinite(d)
    sm = np.where(np.isfinite(sm), sm, 0)
    return np.sum(((m - d) ** 2 - sm ** 2)[ok] / sd2[ok])


def observables(r):
    o = {}
    c = r['counts']
    o['⟨w⟩'] = c.get('sumw_shw', 0) / c['sumw']
    o['tilt ratio ε1(c)/ε1(u)'] = odd(r, 'eps1_cq_all') / odd(r, 'eps1_uq_all')
    o['v1(c)'] = odd(r, 'v1_cq_all')
    o['D0 dv1/dy 10–80 %'] = r.get('slope_D0avg_10-80', [np.nan])[0]     # absent when a small sample has no fit
    o['Δv1 D0−D̄0 10–80 %'] = r.get('slope_dD0_10-80', [np.nan])[0]
    chi = 0.0
    for cl in ('10-40', '0-80'):
        t = S_V2[cl]
        mv, me = at(centres(r, 'PT_V2'), *g(r, f'v2SPS_D0_{cl}'), t[:, 0])
        sd2 = t[:, 2] ** 2 + t[:, 4] ** 2
        if cl == '0-80':
            sd2 = sd2 + np.where(mv < t[:, 1], t[:, 7] ** 2, 0)
        chi += c2inf(mv, me, t[:, 1], sd2)
    for cl in ('0-10', '10-40', '40-80'):
        t = S_RAA[cl]
        sy = 0.5 * (np.abs(t[:, 7]) + np.abs(t[:, 8]))
        sel = (t[:, 3] > 0.9) & (t[:, 3] < 8)
        rr, re_ = ratio(*g(r, f'dN_D0_{cl}_perc'), *g(RPP, 'dN_D0_mb_perc'))
        mv, me = at(centres(r, 'PT_RAA'), rr, re_, t[:, 3])
        for p, m in zip(t[sel, 3], mv[sel]):
            o[f'R_AA {cl} % pT {p:.2f}'] = m
        o[f'χ²∞ R_AA {cl} %'] = c2inf(mv[sel], me[sel], t[sel, 4], t[sel, 5] ** 2 + sy[sel] ** 2)
        chi += o[f'χ²∞ R_AA {cl} %']
    t = S_LCPT
    q, qe = ratio(*g(r, 'dN_Lc_lcbins_10-80_perc'), *g(r, 'dN_D0_lcbins_10-80_perc'))
    for p, m in zip(t[:, 0], q):
        o[f'Λc/D0 pT {p:.2f}'] = m
    chi += c2inf(q, qe, t[:, 1], t[:, 2] ** 2 + t[:, 4] ** 2)
    t = S_LCC
    qs = [ratio(*(x[:1] for x in g(r, f'dN_Lc_3to6_{cl}_perc')), *(x[:1] for x in g(r, f'dN_D0_3to6_{cl}_perc')))
          for cl in ('50-80', '20-50', '0-20')]
    for cl, (a, _) in zip(('50-80', '20-50', '0-20'), qs):
        o[f'Λc/D0 {cl} %'] = a[0]
    chi += c2inf(np.array([a[0] for a, _ in qs]), np.array([b[0] for _, b in qs]), t[:, 1], t[:, 2] ** 2 + t[:, 4] ** 2)
    o['χ²∞ sum (46)'] = chi
    return o


R1 = load(SETS[0])
X = [observables(R1)] + [observables(load(f)) for f in SETS[1:]]
name = os.path.basename(SETS[0]).split('__')[0]
nev = R1['counts']['events']
L = [f'# EPS09s NLO error band (Hessian, 90 % CL) on {name} '
     f'({f"{nev / 1e6:.1f}M" if nev >= 1e5 else f"{nev} events,"} Au+Au; '
     f'p+p ref {os.path.basename(RPP["arm"])})', '',
     'Same events for every set (offline re-weighting, gluon legs; quark legs, ≈ 4 % of legs, at the central factor). '
     'Central = set 1 = the in-code weight. Columns: central, +, −, and the eigen-pair (k = 1..15) with the largest '
     'contribution. analysis/eps_hessian.py.', '',
     '| observable | central | +δ | −δ | dominant pair |', '|---|---|---|---|---|']
for k in X[0]:
    x0 = X[0][k]
    up = np.array([max(X[2 * i + 1][k] - x0, X[2 * i + 2][k] - x0, 0) for i in range(15)])
    dn = np.array([max(x0 - X[2 * i + 1][k], x0 - X[2 * i + 2][k], 0) for i in range(15)])
    dom = int(np.argmax(up ** 2 + dn ** 2)) + 1
    L.append(f'| {k} | {x0:.4f} | {np.sqrt(np.sum(up ** 2)):.4f} | {np.sqrt(np.sum(dn ** 2)):.4f} | {dom} '
             f'(sets {2 * dom}/{2 * dom + 1}) |')
open(OUT, 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
