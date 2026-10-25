#!/usr/bin/env python3
"""indep_check.py AADIR PPDIR OUT.json : independent re-implementation of the headline observables (claims audit,
2026-10-25).  Written from the file-format and method documentation only (docs/CHARM_DAT_FORMAT.md,
docs/FLOW_METHOD.md, STAR papers); it imports nothing from dmlib.py or analyze_arm.py, so a bookkeeping bug there
would show up as a disagreement on the same events.

Observables (all D0 = D0 + D0bar after feed-down of D*0, D*+ ; |y| < 1 unless stated):
  v1(y) of D0+D0bar, 10-80 %, pT > 1.5 GeV, y bins of 0.4 in |y| < 0.8, vs the true RP (+x); slope through 0
  Delta v1 slope (D0 - D0bar)
  v2{RP} and v2{SP} (STAR eta-sub geometry, resolution per 5 % bin) of D0, 10-40 % and 0-80 %
  per-c-quark D0 yields in STAR's R_AA bins (0-10, 10-40, 40-80 %) and in p+p -> R_AA
  Lambda_c/D0 in 10-80 % (STAR pT bins) and in 3 < pT < 6 GeV for 0-20, 20-50, 50-80 %
Weights: charm observables x NCOLT (E field 15) x nPDF factor (E field 17); the charm count x NCOLT only.
Centrality: Nch(|eta|<0.5) percentile within the sample, scaled by (bmax/14.9 fm)^2.
Errors: jackknife over files.  Feed-down random numbers follow the documented production convention (one
numpy default_rng per file, seed 20261017 + 7919 x file index, three draws per decaying hadron: branch, cos theta,
phi), so the same decays are drawn and any difference from analyze_arm.py is a real discrepancy.  With
INDEP_RNG=own a private generator is used instead (statistical agreement only).
"""
import gzip
import json
import math
import os
import sys

import warnings

import numpy as np

AADIR, PPDIR, OUT = sys.argv[1:4]
warnings.filterwarnings('ignore', category=RuntimeWarning)   # empty bins in small samples give NaN
BMAX = 14.5
FRAC = min(1.0, (BMAX / 14.9) ** 2)
rng = np.random.default_rng(12345)
NOE = [0]

MASS = {421: 1.86484, 411: 1.86966, 431: 1.96835, 4122: 2.28646, 111: 0.13498, 211: 0.13957, 22: 0.0}
# strong/EM feed-down (PDG): D*0 -> D0 pi0 / D0 gamma ; D*+ -> D0 pi+ (67.7 %), else D+ ;
# Sigma_c(2455), Sigma_c(2520) -> Lambda_c pi (100 %).  Cumulative branching fractions.
DSTAR = {423: [(0.647, 421, 111), (1.0, 421, 22)], 413: [(0.677, 421, 211), (0.984, 411, 111), (1.0, 411, 22)],
         433: [(0.936, 431, 22), (0.994, 431, 111), (1.0, 431, 22)],
         4112: [(1.0, 4122, 211)], 4212: [(1.0, 4122, 111)], 4222: [(1.0, 4122, 211)],
         4114: [(1.0, 4122, 211)], 4214: [(1.0, 4122, 111)], 4224: [(1.0, 4122, 211)]}


def decay_daughter(kf, p, m, rng):
    """D* -> D + X, isotropic in the rest frame; returns (kf_D, px, py, pz, m_D)"""
    r = rng.uniform()
    for cum, d1, d2 in DSTAR[abs(kf)]:
        if r <= cum:
            break
    m1, m2 = MASS[d1], MASS[d2]
    q = math.sqrt(max((m * m - (m1 + m2) ** 2) * (m * m - (m1 - m2) ** 2), 0)) / (2 * m)
    ct = rng.uniform(-1, 1); st = math.sqrt(1 - ct * ct); ph = rng.uniform(0, 2 * math.pi)
    k = np.array([q * st * math.cos(ph), q * st * math.sin(ph), q * ct]); e1 = math.sqrt(q * q + m1 * m1)
    P = np.array(p); E = math.sqrt(P @ P + m * m); bv = P / E; b2 = bv @ bv
    if b2 > 0:
        g = 1 / math.sqrt(1 - b2); bk = bv @ k
        k = k + ((g - 1) * bk / b2 + g * e1) * bv
    return (int(math.copysign(d1, kf)), k[0], k[1], k[2], m1)


def is_charm_hadron(kf):
    a = abs(kf)
    if a < 1000:
        return (a // 100) % 10 == 4 or (a // 10) % 10 == 4
    return '4' in str(a)[:-1]


def c_count(kf):
    """number of c plus cbar quarks carried by a hadron"""
    a = abs(kf); s = str(a)
    if a < 1000:
        h, t = (a // 100) % 10, (a // 10) % 10
        return (h == 4) + (t == 4)
    return s[:-1].count('4')


def read_events(red, charm):
    E = {}
    if charm:
        for line in gzip.open(charm, 'rt'):
            if line.startswith('E '):
                f = line.split()
                E[int(f[1])] = f          # last block of a repeated event wins
    evs, cur = [], None
    for line in gzip.open(red, 'rt'):
        f = line.split()
        if len(f) == 13:
            cur = dict(iev=int(f[0]), b=float(f[3]), nch=int(f[11]), had=[], Q=None)
            evs.append(cur)
        elif len(f) == 9 and is_charm_hadron(int(f[0])):
            cur['had'].append((int(f[0]),) + tuple(float(x) for x in f[1:4]) + (float(f[4]),))
        elif f and f[0] == 'Q':
            cur['Q'] = np.array(f[2:14], float).reshape(4, 3)
    out = []
    for e in evs:
        f = E.get(e['iev'])
        if f is not None and ((len(f) > 16 and int(f[16]) > 0) or not np.isfinite(float(f[13]))):
            continue                       # events with non-finite partons (documented drop rule)
        e['ncoll'] = float(f[15]) if f is not None and len(f) > 15 else None
        if e['ncoll'] is None and charm:
            NOE[0] += 1                    # Au+Au event without an E record (reported, weight 1)
        e['wpdf'] = float(f[17]) if f is not None and len(f) > 17 else 1.0
        out.append(e)
    return out


def files(d, tag):
    r = sorted(x for x in os.listdir(d) if x.startswith('red_'))
    return [(os.path.join(d, x), os.path.join(d, x.replace('red_', 'charm_'))) for x in r]


Y_V1 = [-0.8, -0.4, 0.0, 0.4, 0.8]
PT_V2 = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 8.0]
PT_RAA = [0.0, 0.485, 0.971, 1.469, 1.951, 2.469, 2.951, 3.769, 4.951, 5.789, 7.391, 9.889]
PT_LC = [2.5, 3.5, 5.0, 8.0]


def bin_of(x, e):
    for i in range(len(e) - 1):
        if e[i] <= x < e[i + 1]:
            return i
    return -1


def rapidity(px, py, pz, m):
    en = math.sqrt(px * px + py * py + pz * pz + m * m)
    return 0.5 * math.log((en + pz) / (en - pz))


class Sum:
    """per-file sums S[name][file] = (num array, den array)"""
    def __init__(self, nf):
        self.nf, self.d = nf, {}

    def add(self, name, nb, ifile, b, x, w):
        if b < 0:
            return
        if name not in self.d:
            self.d[name] = (np.zeros((self.nf, nb)), np.zeros((self.nf, nb)))
        self.d[name][0][ifile, b] += w * x
        self.d[name][1][ifile, b] += w


def ratio_jk(num, den):
    N, D = num.sum(0), den.sum(0)
    nf = num.shape[0]
    with np.errstate(invalid='ignore', divide='ignore'):
        v = N / D
        jk = (N - num) / (D - den)
    return v, np.sqrt((nf - 1) / nf * np.nansum((jk - v) ** 2, 0)), jk


def analyse(d, sysname):
    fl = files(d, sysname)
    nf = len(fl)
    S = Sum(nf)
    allev = [read_events(r, c) for r, c in fl]
    nch = np.sort(np.concatenate([[e['nch'] for e in ev] for ev in allev]))[::-1]
    for i, ev in enumerate(allev):
        frng = rng if os.environ.get('INDEP_RNG') == 'own' else np.random.default_rng(20261017 + 7919 * i)
        for e in ev:
            if sysname == 'pp':
                cl = ['mb']; cen = None
            else:
                # percentile: midpoint of the tied range of this Nch in the descending list
                lo = np.searchsorted(-nch, -e['nch'], 'left'); hi = np.searchsorted(-nch, -e['nch'], 'right')
                cen = 100 * FRAC * 0.5 * (lo + hi) / len(nch)
                cl = [k for k, (a, b) in dict([('0-10', (0, 10)), ('10-40', (10, 40)), ('40-80', (40, 80)),
                      ('10-80', (10, 80)), ('0-80', (0, 80)), ('0-20', (0, 20)), ('20-50', (20, 50)),
                      ('50-80', (50, 80))]).items() if a <= cen < b]
            wn = e['ncoll'] if (sysname == 'auau' and e['ncoll'] is not None) else 1.0
            wq = wn * e['wpdf']
            ncq = sum(c_count(h[0]) for h in e['had'])
            for k in cl:
                S.add(f'nc_{k}', 1, i, 0, ncq, wn)
            # feed-down
            H = []
            for h in e['had']:
                if abs(h[0]) in DSTAR:
                    H.append(decay_daughter(h[0], h[1:4], h[4], frng))
                else:
                    H.append(h)
            # SP / EP ingredients (STAR eta-sub: A = -1<eta<-0.05, B = 0.05<eta<1, pT-weighted Q)
            sp = None
            if sysname == 'auau' and cen is not None and cen < 80 and e['Q'] is not None and np.all(e['Q'][:, 0] > 0):
                qa = (e['Q'][0, 1] + 1j * e['Q'][0, 2]) / e['Q'][0, 0]
                qb = (e['Q'][1, 1] + 1j * e['Q'][1, 2]) / e['Q'][1, 0]
                jc = int(min(cen, 79.999) // 5)
                sp = (qa, qb, jc)
                for k in cl:
                    S.add(f'res_{k}', 16, i, jc, (qa * np.conj(qb)).real, 1.0)
            for kf, px, py, pz, m in H:
                pt = math.hypot(px, py)
                y = rapidity(px, py, pz, m)
                phi = math.atan2(py, px)
                d0 = abs(kf) == 421
                for k in cl:
                    if d0 and pt > 1.5:
                        b = bin_of(y, Y_V1)
                        S.add(f'v1_D0avg_{k}', 4, i, b, math.cos(phi), wq)
                        S.add(f'v1_{"D0" if kf > 0 else "D0bar"}_{k}', 4, i, b, math.cos(phi), wq)
                    if abs(y) >= 1:
                        continue
                    if d0:
                        b = bin_of(pt, PT_V2)
                        S.add(f'v2RP_{k}', 8, i, b, math.cos(2 * phi), wq)
                        if sp is not None and b >= 0:
                            p = math.sqrt(px * px + py * py + pz * pz)
                            eta = 0.5 * math.log((p + pz) / (p - pz))
                            qo = sp[0] if eta > 0 else sp[1]          # opposite hemisphere
                            S.add(f'v2SP_{k}', 8 * 16, i, b * 16 + sp[2], (np.exp(2j * phi) * np.conj(qo)).real, wq)
                        S.add(f'dN_D0_{k}', 11, i, bin_of(pt, PT_RAA), 1.0, wq)
                        S.add(f'dN_D0lc_{k}', 3, i, bin_of(pt, PT_LC), 1.0, wq)
                        if 3 < pt < 6:
                            S.add(f'dN_D0_3to6_{k}', 1, i, 0, 1.0, wq)
                    if abs(kf) == 4122:
                        S.add(f'dN_Lclc_{k}', 3, i, bin_of(pt, PT_LC), 1.0, wq)
                        if 3 < pt < 6:
                            S.add(f'dN_Lc_3to6_{k}', 1, i, 0, 1.0, wq)
    return S, nf, sum(len(ev) for ev in allev)


def per_c(S, nm, k, nb):
    """yield per charm quark in bins, jackknife (numerator sums of w, denominator sums of w*nc); NaN when the
    sample has no entry at all (small test productions)"""
    if f'nc_{k}' not in S.d:
        nan = np.full(nb, np.nan)
        return nan, nan, None
    num = S.d[nm][1] if nm in S.d else np.zeros((S.nf, nb))     # no such hadron in the sample: yield 0                          # sum of w over hadrons in the bin
    cq = S.d[f'nc_{k}'][0][:, 0]              # sum of wn * nc
    return ratio_jk(num, np.repeat(cq[:, None], num.shape[1], 1))


res = {}
SA, nfa, nea = analyse(AADIR, 'auau')
SP_, nfp, nep = analyse(PPDIR, 'pp')
res['events'] = dict(auau=nea, pp=nep, files_auau=nfa, files_pp=nfp, events_without_E_record=NOE[0])
yc = np.array([-0.6, -0.2, 0.2, 0.6])


def slope(v, e):
    w = 1 / e ** 2
    return float(np.sum(w * yc * v) / np.sum(w * yc * yc))


for nm in ('v1_D0avg_10-80', 'v1_D0_10-80', 'v1_D0bar_10-80'):
    if nm not in SA.d:
        continue
    v, e, jk = ratio_jk(*SA.d[nm])
    res[nm] = dict(val=v.tolist(), err=e.tolist())
    k = slope(v, e); ks = np.array([slope(j, e) for j in jk])
    res['slope_' + nm[3:]] = [k, float(np.sqrt((nfa - 1) / nfa * np.sum((ks - k) ** 2)))]
v0, e0, j0 = ratio_jk(*SA.d['v1_D0_10-80']) if 'v1_D0_10-80' in SA.d else (np.full(4, np.nan),) * 3
v1, e1, j1 = ratio_jk(*SA.d['v1_D0bar_10-80']) if 'v1_D0bar_10-80' in SA.d else (np.full(4, np.nan),) * 3
ee = np.hypot(e0, e1); k = slope(v0 - v1, ee)
ks = np.array([slope(a - b, ee) for a, b in zip(j0, j1)])
res['slope_dD0_10-80'] = [k, float(np.sqrt((nfa - 1) / nfa * np.sum((ks - k) ** 2)))]
for k in ('10-40', '0-80'):
    if f'v2RP_{k}' not in SA.d or f'v2SP_{k}' not in SA.d:
        continue
    v, e, _ = ratio_jk(*SA.d[f'v2RP_{k}'])
    res[f'v2_D0_{k}'] = dict(val=v.tolist(), err=e.tolist())
    # SP: per 5 % bin j, <Re u q*>_j / sqrt<Re qA qB*>_j, combined as sum_j N_j/R_j / sum_j D_j
    num, den = SA.d[f'v2SP_{k}']; rn, rd = SA.d[f'res_{k}']

    def spv(num, den, rn, rd):
        R2 = rn.sum(0) / rd.sum(0)
        R = np.where(R2 > 0, np.sqrt(np.abs(R2)), np.nan)
        N = num.sum(0).reshape(8, 16); D = den.sum(0).reshape(8, 16)
        ok = np.isfinite(R)
        return np.sum(np.where(ok, N / R, 0), 1) / np.sum(np.where(ok, D, 0), 1)
    v = spv(num, den, rn, rd)
    jk = np.array([spv(np.delete(num, i, 0), np.delete(den, i, 0), np.delete(rn, i, 0), np.delete(rd, i, 0))
                   for i in range(nfa)])
    res[f'v2SPS_D0_{k}'] = dict(val=v.tolist(), err=np.sqrt((nfa - 1) / nfa * np.nansum((jk - v) ** 2, 0)).tolist())
pv, pe, _ = per_c(SP_, 'dN_D0_mb', 'mb', 11)
res['pp_dN_D0_perc'] = dict(val=pv.tolist(), err=pe.tolist())
for k in ('0-10', '10-40', '40-80'):
    av, ae, _ = per_c(SA, f'dN_D0_{k}', k, 11)
    res[f'dN_D0_{k}_perc'] = dict(val=av.tolist(), err=ae.tolist())
    with np.errstate(invalid='ignore', divide='ignore'):
        res[f'RAA_{k}'] = dict(val=(av / pv).tolist(), err=(av / pv * np.hypot(ae / av, pe / pv)).tolist())
lv, le, _ = per_c(SA, 'dN_Lclc_10-80', '10-80', 3); dv, de, _ = per_c(SA, 'dN_D0lc_10-80', '10-80', 3)
res['LcD0_pt_10-80'] = dict(val=(lv / dv).tolist(), err=(lv / dv * np.hypot(le / lv, de / dv)).tolist())
for k in ('50-80', '20-50', '0-20'):
    lv, le, _ = per_c(SA, f'dN_Lc_3to6_{k}', k, 1); dv, de, _ = per_c(SA, f'dN_D0_3to6_{k}', k, 1)
    res[f'LcD0_3to6_{k}'] = [float(lv[0] / dv[0]), float(lv[0] / dv[0] * math.hypot(le[0] / lv[0], de[0] / dv[0]))]
json.dump(res, open(OUT, 'w'), indent=1)
print(json.dumps(res['events']), 'slope D0avg', res['slope_D0avg_10-80'])
