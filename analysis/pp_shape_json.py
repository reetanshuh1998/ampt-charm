#!/usr/bin/env python3
"""pp_shape_json.py LABEL=PP.json ... : V-D1 p+p D0 pT-shape test against STAR (PRC 99 034908, HEPData ins1711377)
from the analysis output instead of raw events (pp_shape.py loads all events and was OOM-killed on KEKCC).
Model: dN_D0_mb_perc = D0+D0bar per c quark per PT_RAA bin at |y| < 1 (feed-down applied, jackknife errors) ->
density / bin width, log-linearly interpolated to STAR's pT points, invariant yield 1/(2 pi pT) dN/dpT/dy.
Normalisation free (as pp_shape.py: the absolute triggered sigma_ccbar is not a model output here).
chi2 with STAR stat (+) average sys (+) model error, ndf = N - 1."""
import csv
import glob
import json
import sys

import numpy as np

D = str(__import__('pathlib').Path(__file__).resolve().parents[1])   # repository root
rows = [r for r in csv.reader(open(glob.glob(f'{D}/data/hepdata/ins1711377/*/D^0spectrainppcollisions.csv')[0]))
        if r and not r[0].startswith('#')][1:]
st = np.array([[float(x) for x in r[:6]] for r in rows])
pt_s, y_s = st[:, 0], st[:, 1]
e_s = np.hypot(st[:, 2], 0.5 * (st[:, 4] + np.abs(st[:, 5])))
print('STAR p+p D0 points (GeV):', ' '.join(f'{p:.2f}' for p in pt_s))
for spec in sys.argv[1:]:
    lab, fn = spec.split('=', 1)
    r = json.load(open(fn))
    ed = np.array(r['edges']['PT_RAA']); c = 0.5 * (ed[1:] + ed[:-1]); w = np.diff(ed)
    v = np.array(r['dN_D0_mb_perc']['val']); e = np.array(r['dN_D0_mb_perc']['err'])
    dens, edens = v / w, e / w
    ok = dens > 0
    m = np.exp(np.interp(pt_s, c[ok], np.log(dens[ok])))
    me = m * np.interp(pt_s, c[ok], edens[ok] / dens[ok])
    inv, einv = m / (2 * np.pi * pt_s) / 2.0, me / (2 * np.pi * pt_s) / 2.0
    wt = 1 / (e_s ** 2 + (einv * y_s / inv) ** 2)
    a = np.sum(wt * y_s * inv) / np.sum(wt * inv * inv)
    chi2 = np.sum(wt * (y_s - a * inv) ** 2)
    meanpt = np.sum(c * v) / np.sum(v)
    print(f'{lab:22s} <pT>(|y|<1) = {meanpt:.2f} GeV  chi2 = {chi2:5.1f} for {len(y_s) - 1} dof  model/STAR: '
          + ' '.join(f'{x:.2f}' for x in a * inv / y_s))
