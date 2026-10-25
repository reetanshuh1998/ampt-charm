#!/usr/bin/env python3
"""analyze_arm.py ARM_DIR OUT.json [--sys auau|pp] [--bmax B] [--bcent] [--notrig] [-j N]

v5 (2026-10-22): event-plane v2 in the experiments' own methods (v2EPS_*: STAR eta-sub;
v2EPF_*: PHENIX 1<|eta|<2.8), resolution-corrected per 5 % centrality bin (R_epres_*).
v9 (2026-10-23): + N_y1{D0,Dp,Ds,Lc,Xic,Omc} midrapidity ground-state yields (hadrochemistry).
v8 (2026-10-23): --shw none|hij:S re-weights the shadowing offline from the per-event '# hqshw' inputs in
hq_*.log.gz (x1, x2, s_A, s_B; HIJING 1.383 factor, HQSHAD-scaled by S; self-check against the logged w).
v7 (2026-10-23): src_v15 shadowing weight (E field 18) on charm observables, not on the c-quark count.
v6 (same day): scalar product v2SPS_* / v2SPF_* (Luzum-Ollitrault, rms limit, resolution-free)
and the OPV resolution parameter chi_s and alpha per 5 % bin in R_epres_*.
v2FAC_D0_*: factorised D0 v2{EP} = sum_j kappa_j N_j / sum_j D_j with kappa_j = v2{EP}/v2{RP} of
charged hadrons (1-3 GeV) per 5 % bin; v2EPS_D0_* is its noisier direct cross-check.
Every earlier key is unchanged (no random numbers are drawn by the new code).

v3 (2026-10-21): STAR-window bulk v1 (v1s_*, slopeS_*), light-parton v1, Z-record gate,
per-event yields at |y|<0.5 (yld_*) and Nch(|eta|<0.5) (nch5_*) per class.

One pass over an arm (red_* + charm_* files) -> every framework observable and
gate, with leave-one-file-out jackknife errors (docs/VALIDATION.md).

v2 (2026-10-18): binary-collision weighting.  A triggered arm (AMPT_HQTRIG)
carries one hard c-cbar pair per event at every b, while charm production
scales with N_coll(b).  Every CHARM observable of a triggered Au+Au arm is
therefore weighted by <N_coll>(b) (ana/ncoll_b_auau200.txt, HIJING's own NCOLT),
so a centrality class averages D mesons as the data (and an untriggered run) do.
Light hadrons are never weighted.  --notrig (dm1_pub) switches the weight off.
Files are processed in parallel (-j).  v1 of this script: analyze_arm_v1.py.bak.

Centrality (auau): Nch(|eta|<0.5) percentile inside the arm times
(BMAX/14.9 fm)^2; --bcent uses b (c = (b/14.9)^2).
"""
import argparse
import gzip
import os
import json
import math
import pickle
import subprocess
import sys
from multiprocessing import Pool

import numpy as np

HERE = __import__('os').path.dirname(__import__('os').path.abspath(__file__))   # ana/ locally, ana_v4/ on KEKCC
sys.path.insert(0, HERE)
import dmlib as L  # noqa: E402

Y_V1 = np.array([-0.8, -0.4, 0.0, 0.4, 0.8])
Y_V1F = np.arange(-2.0, 2.01, 0.4)
Y_Q = np.arange(-4.0, 4.01, 0.5)
PT_V2 = np.array([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 8.0])
PT_RAA = np.array([0.0, 0.485, 0.971, 1.469, 1.951, 2.469, 2.951, 3.769, 4.951, 5.789, 7.391, 9.889])
PT_LC = np.array([2.5, 3.5, 5.0, 8.0])
# Ds/D0 (STAR PRL 127 092301): fine bins, merged to STAR's 1.0-2.0 / 1.5-2.5 / 2.5-3.5 / 3.5-5 / 5-8 in compare_final
PT_DS = np.array([1.0, 1.5, 2.0, 2.5, 3.5, 5.0, 8.0])
PT_L = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0])
PT_V1 = np.array([0.0, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0])
KET = np.array([0.0, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0])
# v3 (2026-10-21) bulk directed flow with STAR's own windows: PRL 108 202301 (10-70 %,
# pi 0.15-0.75, K 0.2-0.6, p/pbar 0.4-1.0 GeV) and the kaons of PRL 123 162301 (10-80 %,
# pT > 0.2 GeV); slope = linear fit through the origin over |y| < 0.9.
Y_S = np.array([-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9])
V1S_SP = (('pi', (211, -211), 0.15, 0.75), ('K12', (321, -321), 0.2, 0.6), ('K', (321, -321), 0.2, 5.0),
          ('p', (2212,), 0.4, 1.0), ('pbar', (-2212,), 0.4, 1.0))
CLS_V1 = {'10-70': (10, 70), '10-80': (10, 80), 'c0-5': (0, 5), 'c5-10': (5, 10), 'c10-20': (10, 20),
          'c20-30': (20, 30), 'c30-40': (30, 40), 'c40-50': (40, 50), 'c50-60': (50, 60),
          'c60-70': (60, 70), 'c70-80': (70, 80)}
CLS_AA = {'0-10': (0, 10), '10-40': (10, 40), '40-80': (40, 80), '10-80': (10, 80), '0-80': (0, 80),
          '0-20': (0, 20), '20-50': (20, 50), '50-80': (50, 80), 'mb': (-1e9, 1e9)}
# extra classes of STAR's Ds/D0 vs centrality (Fig. 4); only nc_ and the Ds/D0 bins are filled for them
CLS_DS = {'d10-20': (10, 20), 'd20-40': (20, 40), 'd40-60': (40, 60), 'd60-80': (60, 80)}
_nc = np.loadtxt(f'{HERE}/ncoll_b_auau200.txt')
# v5 (2026-10-22) event-plane v2, the experiments' own estimator (the v2_* keys stay v2{RP}).
#  S: STAR eta-sub (PRL 118 212301 D0; PRC 77 054901 K0S): Psi_2 from charged tracks
#     0.2 < pT < 2 GeV, pT-weighted, in |eta| < 1 on the hemisphere OPPOSITE the particle,
#     beyond |eta| > 0.05; POI |y| < 1; resolution sqrt<cos 2(Psi_A - Psi_B)>.
#  F: PHENIX RXN-like (PRC 85 064914): Psi_2 from 1.0 < |eta| < 2.8 (north, south);
#     POI |eta| < 0.35, averaged over the two planes; resolution sqrt<cos 2(Psi_N - Psi_S)>.
# The resolution correction is applied per 5 % centrality bin (CB_EP) before a class is
# combined, as STAR does "for each centrality interval".  Files from reduce_dmeson_v2.awk
# carry the planes' Q-vectors from every charged hadron; for older files the plane uses
# the kept hadrons (all K+-, p, pbar; 25 % of the pions at |y| < 1.5; 10 % of pi/K beyond)
# with their prescale weights, so there the model resolution is below the data's.
CB_EP = np.arange(0.0, 80.01, 5.0)
NCB = len(CB_EP) - 1
# v6 (2026-10-22) scalar product (Luzum & Ollitrault PRC 87 044907): v2{SP} = <Re(u q_opp*)> /
# sqrt<Re(q_A q_B*)> with normalised sub-event flow vectors q, same geometry as the event planes;
# resolution-independent (rms limit).  2-D accumulator prefix -> (output prefix, resolution key):
TWO_D = {'v2S_': ('v2EPS_', 'epres_S'), 'v2F_': ('v2EPF_', 'epres_F'),
         'v2U_': ('v2SPS_', 'spres_S'), 'v2V_': ('v2SPF_', 'spres_F')}


def opv_alpha(rsub):
    """Ollitrault, Poskanzer & Voloshin PRC 80 014904: sub-event resolution R_sub = (sqrt(pi)/2)
    chi e^(-chi^2/2) [I0 + I1](chi^2/2) -> chi_s; v{subEP}^2 = <v>^2 + (alpha - 1) sigma_v^2 with
    alpha = 2 - 4 i1^2 / (i0 + i1)^2, i = I(chi_s^2/2) (their Eqs. 22-23)."""
    def i1(x):
        return sum((x / 2) ** (2 * k + 1) / (math.factorial(k) * math.factorial(k + 1)) for k in range(30))

    def rofchi(c):
        x = c * c / 2
        return np.sqrt(np.pi) / 2 * c * np.exp(-x) * (np.i0(x) + i1(x))
    out = []
    for r in rsub:
        if not np.isfinite(r) or r <= 0 or r >= 0.99:
            out.append((np.nan, np.nan)); continue
        lo, hi = 0.0, 10.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if rofchi(mid) < r else (lo, mid)
        c = 0.5 * (lo + hi); x = c * c / 2
        out.append((c, 2 - 4 * i1(x) ** 2 / (np.i0(x) + i1(x)) ** 2))
    return out


def NB(e):
    return len(e) - 1


def ncoll(b):
    return float(np.interp(b, _nc[:, 0], _nc[:, 1]))


def binof(x, edges):
    i = np.searchsorted(edges, x, side='right') - 1
    return np.where((i >= 0) & (i < len(edges) - 1), i, -1)


def pseudorap(px, py, pz):
    p = np.sqrt(px * px + py * py + pz * pz)
    with np.errstate(divide='ignore', invalid='ignore'):
        return 0.5 * np.log((p + pz) / (p - pz))


def planes(e):
    """second-order event planes of one event: angles A (eta<0), B (eta>0) [STAR], N, S
    [PHENIX], from the Q line if the reduced file has one, else from the kept light hadrons;
    None where a side is empty."""
    if 'Q' in e:   # reduce_dmeson_v2.awk files: every charged hadron, as the detectors see them
        Q = e['Q']
        if np.any(Q[:, 0] <= 0):
            return None
        out = dict(zip('ABNS', (np.angle(q[1] + 1j * q[2]) / 2 for q in Q)))
        out.update({'q' + k: (q[1] + 1j * q[2]) / q[0] for k, q in zip('ABNS', Q)})
        return out
    Lh = e['light']
    if not Lh.size:
        return None
    pt = np.hypot(Lh[:, 1], Lh[:, 2])
    et = pseudorap(Lh[:, 1], Lh[:, 2], Lh[:, 3])
    q = Lh[:, 5] * np.exp(2j * np.arctan2(Lh[:, 2], Lh[:, 1]))
    ok = (pt > 0.2) & (pt < 2.0)
    sides = dict(A=ok & (et < -0.05) & (et > -1.0), B=ok & (et > 0.05) & (et < 1.0),
                 N=ok & (et > 1.0) & (et < 2.8), S=ok & (et < -1.0) & (et > -2.8))
    out = {}
    for nm, s in sides.items():
        if not s.any():
            return None
        # STAR weights its TPC plane by pT; the PHENIX RXN plane is unweighted
        wq_ = (pt[s] if nm in ('A', 'B') else 1.0) * Lh[s, 5]
        Qc = np.sum((pt[s] if nm in ('A', 'B') else 1.0) * q[s])
        out[nm] = np.angle(Qc) / 2
        out['q' + nm] = Qc / np.sum(wq_)     # normalised flow vector (scalar product)
    return out


class Loc:
    """per-file accumulators: name -> [sum w x, sum w, count]"""

    def __init__(self):
        self.d = {}

    def add(self, name, nbin, ib, x, w):
        if name not in self.d:
            self.d[name] = [np.zeros(nbin), np.zeros(nbin), np.zeros(nbin)]
        a = self.d[name]
        ib = np.atleast_1d(ib)
        x = np.broadcast_to(np.atleast_1d(np.asarray(x, float)), ib.shape)
        w = np.broadcast_to(np.atleast_1d(np.asarray(w, float)), ib.shape)
        np.add.at(a[0], ib, w * x)
        np.add.at(a[1], ib, w)
        np.add.at(a[2], ib, 1.0)


def hij_rrx(x, s, h6=0.1, A=197.0):
    # HIJING 1.383 shadowing factor exactly as PYSTFU / HQRRX (charm_mod.f) apply it
    bbr2 = np.minimum(1.0, s * s / 1.44 / A ** 0.66666)
    abx = A ** 0.33333333 - 1.0
    apx = h6 * 4.0 / 3.0 * abx * np.sqrt(1.0 - bbr2)
    aax = 1.192 * np.log(A) ** 0.1666666
    return (aax * (x ** 3 - 1.2 * x ** 2 + 0.21 * x) + 1.0
            - (apx - 1.079 * abx * np.sqrt(x) / np.log(A + 1.0)) * np.exp(-x ** 2 / 0.01))


def shw_weights(rfn, mode):
    """--shw: per-event shadowing weight {iev: w} from hq_*.log.gz (src_v15 '# hqshw' lines), or None.
    'none' -> 1 (also charm-less events, see stale); 'tab:K' -> EPS09s error set K (hqerr_ tables); 'hij:S' -> prod_b [1 + S (R_HIJING(x_b, s_b) - 1)].  Self-check: if the run used HQSHAD = S0
    (logged r1, r2), hij:S0 reproduces the logged w; asserted for S = 0.4 runs."""
    if mode is None:
        return None
    if mode.startswith('tab:'):
        # EPS09s error set K from ana/eps09s_tool/eps09s_errw.py tables (STEP_LOG 99); K = 1 = the in-code weight
        k = int(mode.split(':')[1])
        return {int(l.split()[0]): float(l.split()[k]) for l in
                gzip.open(rfn.replace('/red_', '/hqerr_').replace('.dat.gz', '.txt.gz'), 'rt')}
    hq = rfn.replace('/red_', '/hq_').replace('.dat.gz', '.log.gz')
    out, px = {}, None
    for l in gzip.open(hq, 'rt'):
        if '# hqshw' not in l:
            continue
        f = l.split()
        iev, x1, x2, sa, sb, r1, r2 = int(f[3]), float(f[8]), float(f[9]), float(f[13]), float(f[14]), float(f[16]), float(f[17])
        stale = r1 == 1.0 and r2 == 1.0 and ((x1, x2) == px or f[5] == '0')
        px = (x1, x2)
        if mode == 'none' or stale:
            # stale: HQSHWT never ran (failed charm trigger, nc0 = 0; r = 1 and x repeats the previous event's, or
            # kf = x = 0 on a chunk's first event; 0.6 % of events, STEP_LOG 97) -> w = 1 as in the in-code weighting; no charm, so no charm key is affected
            out[iev] = 1.0
            continue
        S = float(mode.split(':')[1])
        q = [1 + S * (hij_rrx(x, sv) - 1) for x, sv in ((x1, sa), (x2, sb))]
        out[iev] = q[0] * q[1]
        if S == 0.4:
            # s is logged to 1e-3 fm; near the rim (s ~ 1.2 A^1/3) sqrt(1 - s^2/R^2) makes R steep in s, so the check
            # accepts any s within the rounding (STEP_LOG 97)
            for xv, sv, rv in ((x1, sa, r1), (x2, sb, r2)):
                ok = min(abs(1 + S * (hij_rrx(xv, sv + ds) - 1) - rv) for ds in np.linspace(-5e-4, 5e-4, 21))
                assert ok < 2e-3, (hq, iev, q, r1, r2)
    return out


def work(args):
    ifile, rfn, cfn, opt = args
    rng = np.random.default_rng(opt['seed'] + 7919 * ifile)
    reds = L.read_red(rfn)
    chs = L.read_charm(cfn) if cfn else {}
    shw = shw_weights(rfn, opt.get('shw'))
    A = Loc()
    cnt = dict(dropped_nonfinite=0, events=0, events_charm=0, netc_nonzero_events=0, H=0, final_charm_hadrons=0,
               H_final_mismatch_events=0, I=0, F=0, nmatch=0, nfail=0, ccbar_hadrons=0, sumw=0.0, sumw_shw=0.0,
               Z_events=0, Z_nfrz0=0, Z_nfrzt=0, H_dmax_overflow=0)
    extra = dict(dE=[], dmax=[], ncol=[], sig=[], evb=[])
    srt, frac = opt['srt'], opt['frac']
    for e in reds:
        ce = chs.get(e['iev'])
        if ce is not None:
            E = ce['E']
            # stock HIJING occasionally emits a NaN-momentum hadron whose partons
            # poison ZPC (STEP_LOG 35): drop the event entirely.  src_v5 counts
            # non-finite partons in E[16]; older campaigns show it as NaN dE (E[13]).
            bad = (len(E) > 16 and int(E[16]) > 0) or not np.isfinite(float(E[13]))
            if bad:
                cnt['dropped_nonfinite'] += 1
                continue
        cnt['events'] += 1
        b = e['b']
        if opt['sys'] == 'pp':
            cen, incls = -1, ['mb']
        else:
            if opt['bcent']:
                cen = 100 * (b / L.B_MAX_INEL) ** 2
            else:
                hi = np.searchsorted(-srt, -e['nch5'], side='left')
                lo = np.searchsorted(-srt, -e['nch5'], side='right')
                cen = 100 * frac * 0.5 * (hi + lo) / len(srt)
            incls = [k for k, (c0, c1) in CLS_AA.items() if c0 <= cen < c1]
        v1cls = [] if opt['sys'] == 'pp' else [k for k, (c0, c1) in CLS_V1.items() if c0 <= cen < c1]
        dscls = [] if opt['sys'] == 'pp' else [k for k, (c0, c1) in CLS_DS.items() if c0 <= cen < c1]
        # event planes (v5): jc = 5 % centrality bin of the resolution correction
        ep = planes(e) if (opt['sys'] == 'auau' and 0 <= cen < 80) else None
        jc = int(min(cen, 79.999) // 5) if ep else -1
        if ep:
            A.add('epres_S', NCB, jc, np.cos(2 * (ep['A'] - ep['B'])), 1.0)
            A.add('epres_F', NCB, jc, np.cos(2 * (ep['N'] - ep['S'])), 1.0)
            A.add('spres_S', NCB, jc, np.real(ep['qA'] * np.conj(ep['qB'])), 1.0)
            A.add('spres_F', NCB, jc, np.real(ep['qN'] * np.conj(ep['qS'])), 1.0)
        wq = 1.0                                                          # charm weight
        if opt['trig'] and opt['sys'] == 'auau':
            # exact NCOLT of the event (E record field 15, src_v2), else <N_coll>(b)
            wq = float(ce['E'][15]) if (ce is not None and len(ce['E']) > 15) else ncoll(b)
        wn = wq                                                           # weight of the charm count
        if ce is not None and len(ce['E']) > 17:
            # src_v15 AMPT_HQSHW=1: nuclear modification of the pair (E field 18) multiplies every charm
            # observable but not the c-quark count, so per-c-quark yields (R_AA) keep it
            wq *= float(ce['E'][17]) if shw is None else shw[int(ce['E'][1])]
            cnt['sumw_shw'] += wq
        cnt['sumw'] += wn
        extra['evb'].append((e['iev'], round(b, 4)))
        had = [(h[0], h[1], h[2], h[3], h[4]) for h in e['charm']]
        cn = [L.charm_number(h[0]) for h in had]
        nc = sum(abs(c) for c in cn) + 2 * sum(1 for c in cn if c == 0)
        cnt['final_charm_hadrons'] += len(had)
        cnt['ccbar_hadrons'] += sum(1 for c in cn if c == 0)
        cnt['netc_nonzero_events'] += int(sum(cn) != 0)
        cnt['events_charm'] += int(len(had) > 0)
        for k in dscls:
            A.add(f'nc_{k}', 1, 0, nc, wn)
        for k in incls:
            A.add(f'nc_{k}', 1, 0, nc, wn)
            A.add(f'nch5_{k}', 1, 0, float(e['nch5']), 1.0)
        if ce is not None:
            E = ce['E']   # E iev b NP NT NATT npar nc0 nhqs nhqj1 nhqovf nmatch nfail dE sigma
            cnt['I'] += len(ce['I'])
            cnt['F'] += len(ce['F'])
            cnt['H'] += len(ce['H'])
            cnt['nmatch'] += int(E[11])
            cnt['nfail'] += int(E[12])
            for z in ce.get('Z', []):          # src_v9 parton-phase end record
                cnt['Z_events'] += 1
                cnt['Z_nfrzt'] += int(z[3])
                cnt['Z_nfrz0'] += int(z[4])
            extra['dE'].append(float(E[13]))
            if len(E) > 14:
                extra['sig'].append(float(E[14]))
            cnt['H_final_mismatch_events'] += int(len(ce['H']) != len(had))
            # a Fortran field overflow ('*****', seen with the src_v12 catch-all bug) is kept
            # as NaN and counted instead of aborting the file
            for h in ce['H']:
                ovf = '*' in h[14]
                cnt['H_dmax_overflow'] += int(ovf)
                extra['dmax'].append(float('nan') if ovf else float(h[14]))
            for f in ce['F']:
                extra['ncol'].append(int(f[11]))
                A.add('ncol_by_cent10', 8, int(min(max(cen, 0), 79.999) // 10), float(f[11]), 1.0)
            if 0 <= cen < 80 or opt['sys'] == 'pp':
                # quark level (R1): v1 at freeze-out (F, G), eps1 at formation (I, L)
                for tag, recs, wt in (('cq', ce['F'], wq), ('uq', [g for g in ce['G'] if abs(int(g[2])) == 1], 1.0)):
                    if recs:
                        q = np.array([[float(x) for x in r[3:7]] for r in recs])
                        pt = np.hypot(q[:, 0], q[:, 1])
                        y = L.rap(*q.T)
                        c1 = q[:, 0] / np.maximum(pt, 1e-9)
                        by = binof(y, Y_Q)
                        for wn, (p0, p1) in (('all', (0, 5)), ('lo', (0, 1)), ('hi', (1, 5))):
                            s = (by >= 0) & (pt > p0) & (pt < p1)
                            A.add(f'v1_{tag}_{wn}', NB(Y_Q), by[s], c1[s], wt)
                lq = [g for g in ce['G'] if abs(int(g[2])) in (1, 2)]
                if lq:
                    q = np.array([[float(x) for x in r[3:7]] for r in lq])
                    pt = np.hypot(q[:, 0], q[:, 1])
                    y = L.rap(*q.T)
                    by = binof(y, Y_S)
                    s = (by >= 0) & (pt < 5)
                    for k in v1cls:
                        A.add(f'v1s_lq_{k}', NB(Y_S), by[s], (q[s, 0] / np.maximum(pt[s], 1e-9)), 1.0)
                    by = binof(y, Y_Q)
                    s = (by >= 0) & (pt < 5)
                    A.add('v1_lq_all', NB(Y_Q), by[s], (q[s, 0] / np.maximum(pt[s], 1e-9)), 1.0)
                for tag, recs, wt in (('cq', ce['I'], wq), ('uq', [g for g in ce['L'] if abs(int(g[2])) == 1], 1.0)):
                    if recs:
                        q = np.array([[float(x) for x in r[3:11]] for r in recs])
                        en = np.sqrt(np.sum(q[:, :4] ** 2, 1))
                        xf = q[:, 4] + q[:, 0] / en * q[:, 7]
                        yf = q[:, 5] + q[:, 1] / en * q[:, 7]
                        pt = np.hypot(q[:, 0], q[:, 1])
                        ys = L.rap(*q[:, :4].T)
                        c1 = xf / np.maximum(np.hypot(xf, yf), 1e-9)
                        by = binof(ys, Y_Q)
                        for wn, (p0, p1) in (('all', (0, 5)), ('lo', (0, 1)), ('hi', (1, 5))):
                            s = (by >= 0) & (pt > p0) & (pt < p1)
                            A.add(f'eps1_{tag}_{wn}', NB(Y_Q), by[s], c1[s], wt)
        # -------- charm hadrons after feed-down (weight wq)
        fd = L.feed_down(had, rng)
        if fd:
            H = np.array([h[:5] for h in fd], float)
            kf = H[:, 0].astype(int)
            pt = np.hypot(H[:, 1], H[:, 2])
            y = L.rap(*H[:, 1:5].T)
            phi = np.arctan2(H[:, 2], H[:, 1])
            c1, c2 = np.cos(phi), np.cos(2 * phi)
            d0 = np.abs(kf) == 421
            mid = np.abs(y) < 1
            mTd = np.sqrt(pt ** 2 + H[:, 4] ** 2) - H[:, 4]
            for k in incls:
                for nm, sel in (('D0', kf == 421), ('D0bar', kf == -421), ('D0avg', d0)):
                    s = sel & (pt > 1.5)
                    for yy, tg in ((Y_V1, ''), (Y_V1F, 'F')):
                        by = binof(y, yy)
                        ss = s & (by >= 0)
                        A.add(f'v1{tg}_{nm}_{k}', NB(yy), by[ss], c1[ss], wq)
                    bp = binof(pt, PT_V1)
                    ss = sel & (np.abs(y) < 0.8) & (bp >= 0)
                    A.add(f'v1odd_pt_{nm}_{k}', NB(PT_V1), bp[ss], (np.sign(y) * c1)[ss], wq)
                bp = binof(pt, PT_V2)
                s = d0 & mid & (bp >= 0)
                A.add(f'v2_D0_{k}', NB(PT_V2), bp[s], c2[s], wq)
                if ep:
                    pso = np.where(pseudorap(*H[:, 1:4].T) > 0, ep['A'], ep['B'])
                    A.add(f'v2S_D0_{k}', NB(PT_V2) * NCB, bp[s] * NCB + jc, np.cos(2 * (phi - pso))[s], wq)
                    qo = np.where(pseudorap(*H[:, 1:4].T) > 0, ep['qA'], ep['qB'])
                    A.add(f'v2U_D0_{k}', NB(PT_V2) * NCB, bp[s] * NCB + jc, np.real(np.exp(2j * phi) * np.conj(qo))[s], wq)
                    A.add(f'v2R_D0_{k}', NB(PT_V2) * NCB, bp[s] * NCB + jc, c2[s], wq)
                bk = binof(mTd / 2, KET)
                s = d0 & mid & (bk >= 0)
                A.add(f'v2ncq_D0_{k}', NB(KET), bk[s], c2[s] / 2, wq)
                br = binof(pt, PT_RAA)
                s = d0 & mid & (br >= 0)
                A.add(f'dN_D0_{k}', NB(PT_RAA), br[s], 1.0, wq)
                for nm, sel in (('Lc', np.abs(kf) == 4122), ('D0', d0)):
                    bl = binof(pt, PT_LC)
                    s = sel & mid & (bl >= 0)
                    A.add(f'dN_{nm}_lcbins_{k}', NB(PT_LC), bl[s], 1.0, wq)
                    s = sel & mid & (pt > 3) & (pt < 6)
                    A.add(f'dN_{nm}_3to6_{k}', 1, np.zeros(s.sum(), int), 1.0, wq)
                for nm, code in (('D0', 421), ('Dp', 411), ('Ds', 431), ('Lc', 4122)):
                    A.add(f'N_{nm}_{k}', 1, 0, float(np.sum(np.abs(kf) == code)), wq)
                # v9: midrapidity (|y| < 1) ground-state yields for STAR's charm hadrochemistry (dsigma/dy table)
                for nm, codes in (('D0', (421,)), ('Dp', (411,)), ('Ds', (431,)), ('Lc', (4122,)),
                                  ('Xic', (4132, 4232, 4312, 4322, 4314, 4324)), ('Omc', (4332, 4334))):
                    A.add(f'N_y1{nm}_{k}', 1, 0, float(np.sum(np.isin(np.abs(kf), codes) & mid)), wq)
            for k in incls + dscls:
                bd = binof(pt, PT_DS)
                for nm, sel in (('Ds', np.abs(kf) == 431), ('D0', d0)):
                    s = sel & mid & (bd >= 0)
                    A.add(f'dN_{nm}_dsbins_{k}', NB(PT_DS), bd[s], 1.0, wq)
            for k in incls:
                for nm, code in (('Dst0', 423), ('Dstp', 413), ('D0prim', 421), ('Dpprim', 411), ('jpsi', 443), ('etac', 441)):
                    A.add(f'N_{nm}_{k}', 1, 0, float(sum(1 for h in had if abs(h[0]) == code)), wq)
                s = d0 & mid
                A.add(f'meanpt_D0_{k}', 1, np.zeros(s.sum(), int), pt[s], wq)
        # -------- light hadrons (never weighted)
        Lh = e['light']
        if Lh.size:
            kf = Lh[:, 0].astype(int)
            pt = np.hypot(Lh[:, 1], Lh[:, 2])
            y = L.rap(*Lh[:, 1:5].T)
            phi = np.arctan2(Lh[:, 2], Lh[:, 1])
            w = Lh[:, 5]
            c1, c2 = np.cos(phi), np.cos(2 * phi)
            mT = np.sqrt(pt ** 2 + Lh[:, 4] ** 2) - Lh[:, 4]
            if ep:
                etl = pseudorap(Lh[:, 1], Lh[:, 2], Lh[:, 3])
                # a POI is never in its own plane: S uses the opposite hemisphere, F (|eta| < 0.35) both sides
                cS = np.cos(2 * (phi - np.where(etl > 0, ep['A'], ep['B'])))
                cF = 0.5 * (np.cos(2 * (phi - ep['N'])) + np.cos(2 * (phi - ep['S'])))
                u2 = np.exp(2j * phi)
                cU = np.real(u2 * np.conj(np.where(etl > 0, ep['qA'], ep['qB'])))
                cV = 0.5 * np.real(u2 * np.conj(ep['qN'] + ep['qS']))
                # kappa_j = v2{EP}/v2{RP} of charged hadrons, |y| < 1, 1 < pT < 3 GeV (factorised D0 v2{EP})
                s = (np.abs(y) < 1) & (pt > 1.0) & (pt < 3.0) & np.isin(np.abs(kf), (211, 321, 2212))
                A.add('kapRP', NCB, np.full(s.sum(), jc), c2[s], w[s])
                A.add('kapS', NCB, np.full(s.sum(), jc), cS[s], w[s])
            for k in incls:
                for nm, sel in (('K', np.abs(kf) == 321), ('Kp', kf == 321), ('Km', kf == -321),
                                ('pi', np.abs(kf) == 211), ('p', kf == 2212), ('pbar', kf == -2212)):
                    for yy, tg in ((Y_V1, ''), (Y_V1F, 'F'), (Y_Q, 'Q')):
                        by = binof(y, yy)
                        s = sel & (pt > 0.2) & (pt < 5) & (by >= 0)
                        A.add(f'v1{tg}_{nm}_{k}', NB(yy), by[s], c1[s], w[s])
                mid = np.abs(y) < 1
                for nm, sel, nq in (('pi', np.abs(kf) == 211, 2), ('K', np.abs(kf) == 321, 2), ('p', np.abs(kf) == 2212, 3)):
                    bp = binof(pt, PT_L)
                    s = sel & mid & (bp >= 0)
                    A.add(f'v2_{nm}_{k}', NB(PT_L), bp[s], c2[s], w[s])
                    if ep:
                        A.add(f'v2S_{nm}_{k}', NB(PT_L) * NCB, bp[s] * NCB + jc, cS[s], w[s])
                        s = sel & (np.abs(etl) < 0.35) & (bp >= 0)
                        A.add(f'v2F_{nm}_{k}', NB(PT_L) * NCB, bp[s] * NCB + jc, cF[s], w[s])
                        A.add(f'v2V_{nm}_{k}', NB(PT_L) * NCB, bp[s] * NCB + jc, cV[s], w[s])
                        s = sel & mid & (bp >= 0)
                        A.add(f'v2U_{nm}_{k}', NB(PT_L) * NCB, bp[s] * NCB + jc, cU[s], w[s])
                    bk = binof(mT / nq, KET)
                    s = sel & mid & (bk >= 0)
                    A.add(f'v2ncq_{nm}_{k}', NB(KET), bk[s], c2[s] / nq, w[s])
            # yields at |y| < 0.5 per event (bulk chemistry check; STAR PRC 79 034909 K-/pi-)
            mid5 = np.abs(y) < 0.5
            for k in incls:
                for nm, code in (('pim', -211), ('pip', 211), ('Km', -321), ('Kp', 321), ('p', 2212), ('pbar', -2212)):
                    A.add(f'yld_{nm}_{k}', 1, 0, float(np.sum(w[mid5 & (kf == code)])), 1.0)
            by = binof(y, Y_S)
            for nm, codes, p0, p1 in V1S_SP:
                s = np.isin(kf, codes) & (pt > p0) & (pt < p1) & (by >= 0)
                for k in v1cls:
                    A.add(f'v1s_{nm}_{k}', NB(Y_S), by[s], c1[s], w[s])
    return ifile, A.d, cnt, extra


def header_scan(fn):
    out = subprocess.run(f"zcat {fn} | awk 'NF==13{{print $12, $4}}'", shell=True,
                         capture_output=True, text=True).stdout
    return np.array([list(map(float, l.split())) for l in out.splitlines()]).reshape(-1, 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('arm')
    ap.add_argument('out')
    ap.add_argument('--sys', default='auau')
    ap.add_argument('--bmax', type=float, default=14.5)
    ap.add_argument('--bcent', action='store_true')
    ap.add_argument('--notrig', action='store_true', help='untriggered arm: no N_coll weight')
    ap.add_argument('-j', type=int, default=6)
    ap.add_argument('--maxfiles', type=int, default=0)
    ap.add_argument('--seed', type=int, default=20261017)
    ap.add_argument('--shw', help="offline shadowing weight: none | hij:S (needs hq_*.log.gz with '# hqshw')")
    # split mode (2026-10-21, single-core farm jobs): --part K --nparts N analyses files i with
    # i % N == K and pickles their per-file accumulators to OUT; --merge OUT.part* finishes.
    # File seeds depend only on the global index i, so merge(parts) == one full run.
    ap.add_argument('--part', type=int)
    ap.add_argument('--nparts', type=int)
    ap.add_argument('--merge', nargs='+')
    a = ap.parse_args()
    files = L.arm_files(a.arm)
    if a.maxfiles:
        files = files[:a.maxfiles]
    NF = NFILES = len(files)
    print(f'{a.arm}: {NF} files, trig-weight={not a.notrig}', flush=True)
    if a.merge:
        # streaming merge (analysis v13, STEP_LOG 107): one part pickle at a time, per-file accumulators summed into
        # NB <= 2000 contiguous jackknife blocks (delete-a-block jackknife); NB = NF (identical to the old merge) for
        # arms with <= 2000 files.  The 7982-file dm25_main merge needed > 4 GB when every file was kept.
        NB = min(NF, int(os.environ.get('ANA_JKMAX', '2000')))
        acc, cnt, ex, evb, seen, nevb = {}, None, None, {}, [], 0
        for f in a.merge:
            for i, d, c, e in pickle.load(open(f, 'rb')):
                b = i * NB // NF
                seen.append(i)
                for nm, (u, v_, w) in d.items():
                    if nm not in acc:
                        acc[nm] = tuple(np.zeros((NB, len(u))) for _ in range(3))
                    acc[nm][0][b] += u
                    acc[nm][1][b] += v_
                    acc[nm][2][b] += w
                cnt = dict(c) if cnt is None else {k: cnt[k] + c[k] for k in cnt}
                if ex is None:
                    ex = {k: [] for k in e if k != 'evb'}
                for k in ex:                  # per-event/per-hadron lists as compact arrays (v13b)
                    ex[k].append(np.asarray(e[k], dtype=np.float64))
                if evb is not None:
                    evb[str(files[i][0])] = e['evb']
                    nevb += len(e['evb'])
                    if nevb > 3_000_000:
                        evb = None            # pairing record skipped above 3M events (12M-event p+p merges)
        assert sorted(seen) == list(range(NF)), 'merge: parts do not cover every file once'
        ex = {k: np.concatenate(v) if v else np.zeros(0) for k, v in ex.items()}
        NF = NB
    else:
        with Pool(a.j) as P:
            hs = P.map(header_scan, [r for _, r, _ in files])
        nch = np.concatenate([h[:, 0] for h in hs])
        opt = dict(sys=a.sys, bcent=a.bcent, trig=not a.notrig, seed=a.seed, shw=a.shw,
                   srt=np.sort(nch)[::-1], frac=min(1.0, (a.bmax / L.B_MAX_INEL) ** 2))
        sel = [i for i in range(NF) if a.part is None or i % a.nparts == a.part]
        with Pool(a.j) as P:
            parts = P.map(work, [(i, files[i][1], files[i][2], opt) for i in sel], chunksize=1)
        if a.part is not None:
            pickle.dump(parts, open(a.out, 'wb'))
            print(f'part {a.part}/{a.nparts}: {len(sel)} files -> {a.out}')
            return
    if not a.merge:
        names = sorted(set().union(*[p[1].keys() for p in parts]))
        acc = {}
        for nm in names:
            nbin = next(len(p[1][nm][0]) for p in parts if nm in p[1])
            num, den, n = np.zeros((NF, nbin)), np.zeros((NF, nbin)), np.zeros((NF, nbin))
            for i, d, _, _ in parts:
                if nm in d:
                    num[i], den[i], n[i] = d[nm]
            acc[nm] = (num, den, n)
        cnt = {k: sum(p[2][k] for p in parts) for k in parts[0][2]}
        ex = {k: sum((p[3][k] for p in parts), []) for k in parts[0][3]}
        evb = {str(files[p[0]][0]): p[3]['evb'] for p in parts}
    acc = {nm: acc[nm] for nm in sorted(acc)}

    def ratio(num, den):
        N, D = num.sum(0), den.sum(0)
        with np.errstate(invalid='ignore', divide='ignore'):
            v = N / D
            jk = (N[None] - num) / (D[None] - den)
        return v, np.sqrt((NF - 1) / NF * np.nansum((jk - v[None]) ** 2, 0)), jk

    res = dict(arm=a.arm, sys=a.sys, nfiles=NFILES, jk_units=NF, trig_weight=not a.notrig, counts=cnt,
               edges=dict(Y_V1=Y_V1.tolist(), Y_V1F=Y_V1F.tolist(), Y_Q=Y_Q.tolist(), Y_S=Y_S.tolist(), PT_V2=PT_V2.tolist(),
                          PT_RAA=PT_RAA.tolist(), PT_LC=PT_LC.tolist(), PT_DS=PT_DS.tolist(), PT_L=PT_L.tolist(), KET=KET.tolist(),
                          PT_V1=PT_V1.tolist()),
               file_idx=[f[0] for f in files])
    JK = {}

    def epv2(num, den, rn, rd):
        """sum_j (N_j / R_j) / sum_j D_j over the 5 % bins j, R_j = sqrt<cos 2 dPsi>_j;
        a bin whose <cos 2 dPsi> is not positive cannot be corrected and is left out."""
        nb = num.shape[1] // NCB
        with np.errstate(invalid='ignore', divide='ignore'):
            r2 = rn.sum(0) / rd.sum(0)
            R = np.where(r2 > 0, np.sqrt(np.abs(r2)), np.nan)
            N, D = num.sum(0).reshape(nb, NCB), den.sum(0).reshape(nb, NCB)
            okb = np.isfinite(R)
            v = np.sum(np.where(okb, N / R, 0), 1) / np.sum(np.where(okb, D, 0), 1)
            jks = []
            for i in range(NF):
                r2i = (rn.sum(0) - rn[i]) / (rd.sum(0) - rd[i])
                Ri = np.where(okb & (r2i > 0), np.sqrt(np.abs(r2i)), np.nan)
                Ni, Di = N - num[i].reshape(nb, NCB), D - den[i].reshape(nb, NCB)
                oki = np.isfinite(Ri)
                jks.append(np.sum(np.where(oki, Ni / Ri, 0), 1) / np.sum(np.where(oki, Di, 0), 1))
            jk = np.array(jks)
        return v, np.sqrt((NF - 1) / NF * np.nansum((jk - v[None]) ** 2, 0)), jk, D.sum(1)

    for nm, (num, den, n) in acc.items():
        if nm[:4] in TWO_D:
            out_, resk = TWO_D[nm[:4]]
            rn, rd, _ = acc[resk]
            v, e, jk, d = epv2(num, den, rn, rd)
            res[out_ + nm[4:]] = dict(val=v.tolist(), err=e.tolist(), sumw=d.tolist())
            JK[out_ + nm[4:]] = jk
    if 'kapRP' in acc:
        kn, kd, _ = acc['kapS']
        rn_, rd_, _ = acc['kapRP']
        sn, sd, _ = acc['epres_S']

        def kappa(dr=None):
            g = (lambda a: a.sum(0)) if dr is None else (lambda a: a.sum(0) - a[dr])
            with np.errstate(invalid='ignore', divide='ignore'):
                r2 = g(sn) / g(sd)
                R = np.where(r2 > 0, np.sqrt(np.abs(r2)), np.nan)
                return (g(kn) / g(kd)) / R / (g(rn_) / g(rd_))
        k0 = kappa()
        res['kappa_S'] = dict(cb=CB_EP.tolist(), val=k0.tolist())
        for nm, (num, den, n) in acc.items():
            if not nm.startswith('v2R_'):
                continue
            nb = num.shape[1] // NCB

            def fac(kp, dr=None):
                N = (num.sum(0) - (num[dr] if dr is not None else 0)).reshape(nb, NCB)
                D = (den.sum(0) - (den[dr] if dr is not None else 0)).reshape(nb, NCB)
                ok = np.isfinite(kp)
                with np.errstate(invalid='ignore', divide='ignore'):
                    return np.sum(np.where(ok, kp * N, 0), 1) / np.sum(np.where(ok, D, 0), 1)
            v = fac(k0)
            jk = np.array([fac(kappa(i), i) for i in range(NF)])
            res['v2FAC' + nm[3:]] = dict(val=v.tolist(), err=np.sqrt((NF - 1) / NF * np.nansum((jk - v[None]) ** 2, 0)).tolist())
    for nm in ('epres_S', 'epres_F'):
        if nm in acc:
            r2 = acc[nm][0].sum(0) / acc[nm][1].sum(0)
            Rs = np.where(r2 > 0, np.sqrt(np.abs(r2)), np.nan)
            ca = opv_alpha(Rs)
            res['R_' + nm] = dict(cb=CB_EP.tolist(), Rsub=Rs.tolist(), chi_s=[c for c, _ in ca], alpha=[a for _, a in ca])
    for nm, (num, den, n) in acc.items():
        if nm[:4] in TWO_D or nm.startswith('v2R_'):
            continue
        v, e, jk = ratio(num, den)
        res[nm] = dict(val=v.tolist(), err=e.tolist(), n=n.sum(0).tolist())
        JK[nm] = jk
    for nm, (num, den, n) in acc.items():   # yields per charm quark
        if not (nm.startswith('dN_') or nm.startswith('N_')):
            continue
        k = nm.split('_')[-1]
        if f'nc_{k}' not in acc:
            continue
        cq = acc[f'nc_{k}'][0][:, 0]
        v, e, _ = ratio(num, np.repeat(cq[:, None], num.shape[1], 1))
        res[nm + '_perc'] = dict(val=v.tolist(), err=e.tolist(), nc=float(acc[f'nc_{k}'][2].sum()))
    yc = 0.5 * (Y_V1[1:] + Y_V1[:-1])

    def slope(v, e):
        w = 1 / np.maximum(e, 1e-12) ** 2
        return float(np.sum(w * yc * v) / np.sum(w * yc * yc))

    for nm in list(acc):
        if nm.startswith('v1_') and len(res[nm]['val']) == 4:
            v, e = np.array(res[nm]['val']), np.array(res[nm]['err'])
            if np.any(~np.isfinite(v)) or np.any(e <= 0):
                continue
            k = slope(v, e)
            ks = np.array([slope(j, e) for j in JK[nm]])
            res['slope_' + nm[3:]] = [k, float(np.sqrt((NF - 1) / NF * np.nansum((ks - k) ** 2)))]
    ycs = 0.5 * (Y_S[1:] + Y_S[:-1])

    def slope_s(v, e):
        w = 1 / np.maximum(e, 1e-12) ** 2
        return float(np.sum(w * ycs * v) / np.sum(w * ycs * ycs))

    for nm in list(acc):
        if nm.startswith('v1s_'):
            v, e = np.array(res[nm]['val']), np.array(res[nm]['err'])
            if np.any(~np.isfinite(v)) or np.any(e <= 0):
                continue
            k = slope_s(v, e)
            ks = np.array([slope_s(j, e) for j in JK[nm]])
            res['slopeS_' + nm[4:]] = [k, float(np.sqrt((NF - 1) / NF * np.nansum((ks - k) ** 2)))]
    for k in CLS_AA:
        n0, n1 = f'v1_D0_{k}', f'v1_D0bar_{k}'
        if n0 in res and n1 in res:
            v = np.array(res[n0]['val']) - np.array(res[n1]['val'])
            e = np.hypot(res[n0]['err'], res[n1]['err'])
            if np.all(np.isfinite(v)) and np.all(e > 0):
                kk = slope(v, e)
                ks = np.array([slope(a_ - b_, e) for a_, b_ in zip(JK[n0], JK[n1])])
                res[f'slope_dD0_{k}'] = [kk, float(np.sqrt((NF - 1) / NF * np.nansum((ks - kk) ** 2)))]
    if len(ex['dE']):
        res['dE_direct'] = dict(mean=float(np.mean(ex['dE'])), rms=float(np.std(ex['dE'])),
                                max=float(np.max(np.abs(ex['dE']))))
    res['sigma_trig_mb'] = float(np.mean(ex['sig'])) if len(ex['sig']) else None
    if len(ex['dmax']):
        d = np.array(ex['dmax'])
        d = d[np.isfinite(d)]          # overflowed fields are counted in counts['H_dmax_overflow']
        res['coal_dmax_fm'] = dict(median=float(np.median(d)), p90=float(np.percentile(d, 90)),
                                   frac_gt2=float(np.mean(d > 2)),
                                   hist=np.histogram(d, bins=np.linspace(0, 5, 26))[0].tolist())
    if len(ex['ncol']):
        n = np.array(ex['ncol'])
        res['ncol_charm'] = dict(mean=float(n.mean()), frac0=float(np.mean(n == 0)))
    if evb is not None:
        res['event_b'] = evb   # gate G6 (pairing)
    json.dump(res, open(a.out, 'w'))
    print('wrote', a.out, cnt, flush=True)


if __name__ == '__main__':
    main()
