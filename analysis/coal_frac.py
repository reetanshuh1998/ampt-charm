#!/usr/bin/env python3
"""Coalesced fraction of charm mesons vs pT (charm.dat H records, mode 1 = coalescence,
2 = fragmentation), |y| < 1.  Is the R_AA step at pT ~ 2.7 GeV the
m_T threshold (2.0 in p+p, 3.3 in Au+Au)?
usage: coal_frac.py LABEL FILE.gz [FILE.gz ...]"""
import sys, gzip, math
import numpy as np
edges = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8])
lab, files = sys.argv[1], sys.argv[2:]
n = np.zeros((len(edges) - 1, 3))
nd = np.zeros((len(edges) - 1, 3))
for fn in files:
    for l in gzip.open(fn, 'rt'):
        if l[0] != 'H':
            continue
        f = l.split()
        kf = abs(int(f[2]))
        if kf // 100 != 4 or kf >= 1000:
            continue
        px, py, pz, m = map(float, f[3:7])
        e = math.sqrt(px * px + py * py + pz * pz + m * m)
        if abs(0.5 * math.log((e + pz) / (e - pz))) > 1:
            continue
        b = np.searchsorted(edges, math.hypot(px, py)) - 1
        if 0 <= b < len(edges) - 1:
            n[b, int(f[-1])] += 1
            if kf in (421, 423):
                nd[b, int(f[-1])] += 1
print(f'# {lab}: {len(files)} files; charm mesons |y|<1; coalesced fraction = mode1/(mode1+mode2)')
print('# pT_lo pT_hi  N   f_coal(all D)  f_coal(D0,D*0)')
for i in range(len(edges) - 1):
    t = n[i, 1] + n[i, 2]; td = nd[i, 1] + nd[i, 2]
    print(f'{edges[i]:5.1f} {edges[i+1]:5.1f} {int(t):7d}  {n[i,1]/t if t else float("nan"):.3f}  {nd[i,1]/td if td else float("nan"):.3f}')
