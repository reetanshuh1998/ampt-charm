"""dmlib -- readers, feed-down decays and observables for the D-meson framework.

Inputs per arm directory (as written by scripts/run_job.sh):
  red_<arm>_<idx>.dat.gz    reduced ampt.dat (scripts/reduce_dmeson_v2.awk)
  charm_<arm>_<idx>.dat.gz  ana/charm.dat (docs/CHARM_DAT_FORMAT.md)

Conventions: reaction plane = +x
(AMPT puts the projectile at +x moving to +z, iphirp = 0), v_n = <cos n phi>.
Every event is one job-file member; errors are jackknife over job files,
and arms run with the same job index share HIJING seeds (paired).
"""
import glob
import gzip
import os
import re
from collections import defaultdict

import numpy as np

# ---------------------------------------------------------------- masses (GeV)
M = {111: 0.13498, 211: 0.13957, 22: 0.0,
     421: 1.86484, 411: 1.86966, 431: 1.96835,
     423: 2.00685, 413: 2.01026, 433: 2.1122,
     4122: 2.28646}
# two-body feed-down tables (PDG 2024 branching ratios), daughters (charm, other)
DECAYS = {
    423: [(0.647, 421, 111), (0.353, 421, 22)],
    413: [(0.677, 421, 211), (0.307, 411, 111), (0.016, 411, 22)],
    433: [(0.936, 431, 22), (0.058, 431, 111), (0.006, 431, 22)],
    4112: [(1.0, 4122, 211)], 4212: [(1.0, 4122, 111)], 4222: [(1.0, 4122, 211)],
    4114: [(1.0, 4122, 211)], 4214: [(1.0, 4122, 111)], 4224: [(1.0, 4122, 211)],
}

B_MAX_INEL = 14.9   # fm, Au+Au 200 GeV geometric cross section (STAR Glauber b(80%)=13.3)


def is_charm(kf):
    a = abs(kf)
    s = str(a)
    if a < 100:
        return False
    if a < 1000:
        return s[-2] == '4' or s[-3] == '4'
    return '4' in s[:-1]


def charm_number(kf):
    """+1 for a hadron with one c, -1 for one cbar, 0 for c-cbar / none."""
    a = abs(kf)
    s = str(a)
    sg = 1 if kf > 0 else -1
    if a < 1000:                      # meson q qbar: 4 in tens or hundreds
        h, t = int(s[-3]), int(s[-2])
        if h == 4 and t == 4:
            return 0
        if h == 4:                    # c qbar  (D+ 411, D0 421, Ds 431: c is the heavier)
            return sg
        return 0
    n = sum(1 for d in s[:-1] if d == '4')
    return sg * n


# ---------------------------------------------------------------- kinematics
def rap(px, py, pz, m):
    e = np.sqrt(px * px + py * py + pz * pz + m * m)
    return 0.5 * np.log((e + pz) / (e - pz))


def boost_from_rest(p_rest, parent_p, parent_m):
    """Lorentz-boost 4-vectors p_rest (N,4: E,px,py,pz) from the parent rest frame."""
    E = np.sqrt(np.sum(parent_p ** 2) + parent_m ** 2)
    b = parent_p / E
    b2 = np.dot(b, b)
    if b2 <= 0:
        return p_rest
    g = 1.0 / np.sqrt(1.0 - b2)
    bp = p_rest[:, 1:] @ b
    gam2 = (g - 1.0) / b2
    out = np.empty_like(p_rest)
    out[:, 0] = g * (p_rest[:, 0] + bp)
    out[:, 1:] = p_rest[:, 1:] + np.outer(gam2 * bp + g * p_rest[:, 0], b)
    return out


def two_body(parent_p, parent_m, m1, m2, rng):
    """isotropic two-body decay; returns momentum 3-vector of daughter 1"""
    pst = np.sqrt(max((parent_m ** 2 - (m1 + m2) ** 2) * (parent_m ** 2 - (m1 - m2) ** 2), 0.0)) / (2 * parent_m)
    c = rng.uniform(-1, 1)
    s = np.sqrt(1 - c * c)
    ph = rng.uniform(0, 2 * np.pi)
    p1 = np.array([[np.sqrt(pst ** 2 + m1 ** 2), pst * s * np.cos(ph), pst * s * np.sin(ph), pst * c]])
    return boost_from_rest(p1, np.asarray(parent_p, float), parent_m)[0, 1:]


def feed_down(hadrons, rng):
    """hadrons: list of (kf, px, py, pz, m). Returns list of (kf, px, py, pz, m, parent_kf)
    with D*, Ds*, Sigma_c decayed to their ground state (one generation is enough:
    every listed daughter is stable under the strong/EM interaction)."""
    out = []
    for kf, px, py, pz, m in hadrons:
        a = abs(kf)
        if a in DECAYS:
            r = rng.uniform()
            acc = 0.0
            for br, d1, d2 in DECAYS[a]:
                acc += br
                if r <= acc:
                    break
            p = two_body((px, py, pz), m, M[d1], M[d2], rng)
            out.append((int(np.sign(kf)) * d1, p[0], p[1], p[2], M[d1], kf))
        else:
            out.append((kf, px, py, pz, m, 0))
    return out


# ---------------------------------------------------------------- readers
def _open(fn):
    return gzip.open(fn, 'rt') if fn.endswith('.gz') else open(fn)


def read_red(fn):
    """reduced ampt.dat -> list of events: dict(b, nch5, nch10, charm=[(kf,px,py,pz,m,x,y,z,t)],
    light=np.array (kf,px,py,pz,m,w), and Q (4x3: sumw, Qx, Qy of windows A, B, N, S) when the
    file comes from reduce_dmeson_v2.awk)"""
    evs = []
    cur = None
    for line in _open(fn):
        f = line.split()
        if len(f) == 13:
            cur = dict(iev=int(f[0]), b=float(f[3]), nch5=int(f[11]), nch10=int(f[12]), charm=[], light=[])
            evs.append(cur)
        elif len(f) == 9:
            cur['charm'].append((int(f[0]),) + tuple(map(float, f[1:9])))
        elif len(f) == 6:
            cur['light'].append((int(f[0]),) + tuple(map(float, f[1:6])))
        elif f[0] == 'Q':   # reduce_dmeson_v2.awk: full-multiplicity Q2 vectors (A, B, N, S windows)
            cur['Q'] = np.array(f[2:14], float).reshape(4, 3)
    for e in evs:
        e['light'] = np.array(e['light'], float).reshape(-1, 6)
    return evs


def read_charm(fn):
    """charm.dat -> dict iev -> dict(E=[...], I=[...], L, F, G, H). Repeated HIJING
    attempts write a new E block for the same event; the LAST block is the kept event."""
    evs = {}
    cur = None
    for line in _open(fn):
        f = line.split()
        t = f[0]
        if t == 'E':
            cur = dict(E=f, I=[], L=[], F=[], G=[], H=[], Z=[])
            evs[int(f[1])] = cur
        elif cur is not None and t in cur:
            cur[t].append(f)
    return evs


def arm_files(d):
    red = sorted(glob.glob(os.path.join(d, 'red_*.dat.gz')))
    out = []
    for r in red:
        c = r.replace('/red_', '/charm_')
        idx = int(re.search(r'_(\d+)\.dat\.gz$', r).group(1))
        out.append((idx, r, c if os.path.exists(c) else None))
    return out


# ---------------------------------------------------------------- estimators
class Acc:
    """per-file accumulator of (sum w*x, sum w) in bins, for ratio-type observables
    (v_n, yields per charm quark) with leave-one-file-out jackknife errors."""

    def __init__(self, nbins, nfiles):
        self.num = np.zeros((nfiles, nbins))
        self.den = np.zeros((nfiles, nbins))

    def add(self, ifile, ibin, x, w=1.0):
        np.add.at(self.num[ifile], ibin, w * x)
        np.add.at(self.den[ifile], ibin, w)

    def result(self):
        N, D = self.num.sum(0), self.den.sum(0)
        with np.errstate(invalid='ignore', divide='ignore'):
            val = N / D
            nf = self.num.shape[0]
            jk = (N[None] - self.num) / (D[None] - self.den)
            err = np.sqrt((nf - 1) / nf * np.nansum((jk - val[None]) ** 2, 0))
        return val, err, D


def slope_through_origin(y, v, e):
    """weighted LSQ slope of v = k*y (STAR fit for D0 dv1/dy)."""
    w = 1 / e ** 2
    k = np.sum(w * y * v) / np.sum(w * y * y)
    return k, np.sqrt(1 / np.sum(w * y * y))


def slope_jackknife(acc, yc, sel):
    """dv1/dy through origin from an Acc in y bins, jackknife over files (uniform weights
    from the full-sample errors)."""
    val, err, _ = acc.result()
    k, _ = slope_through_origin(yc[sel], val[sel], err[sel])
    w = 1 / err[sel] ** 2
    nf = acc.num.shape[0]
    ks = []
    N, D = acc.num.sum(0), acc.den.sum(0)
    for i in range(nf):
        with np.errstate(invalid='ignore', divide='ignore'):
            vi = ((N - acc.num[i]) / (D - acc.den[i]))[sel]
        ks.append(np.sum(w * yc[sel] * vi) / np.sum(w * yc[sel] ** 2))
    ks = np.array(ks)
    return k, np.sqrt((nf - 1) / nf * np.sum((ks - k) ** 2))
