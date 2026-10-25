# Bulk flow against published 200 GeV Au+Au data

Model: reaction-plane v1/v2 in STAR windows (analyze_arm.py v3). Errors: model jackknife; data stat (+) sys (+) 10 % common for v1.

## Directed-flow slopes dv1/dy

| arm | π 10–70 | K 10–70 | p 10–70 | p̄ 10–70 | K 10–80 | χ²(5) | π vs cent χ²(9) | partons 10–70 | D0 10–80 |
|---|---|---|---|---|---|---|---|---|---|
| **STAR** | -0.0027 | -0.0017 | -0.0015 | -0.0046 | -0.0030 | — | — | — | -0.080 |
| dm25_main | -0.0015 ± 0.0002 | +0.0021 ± 0.0002 | +0.0055 ± 0.0006 | +0.0001 ± 0.0009 | +0.0021 ± 0.0002 | 277.7 | 49.9 | -0.0006 ± 0.0002 | -0.031 ± 0.005 |
| dm25_pub | +0.0019 ± 0.0003 | +0.0042 ± 0.0004 | +0.0034 ± 0.0009 | +0.0015 ± 0.0012 | +0.0041 ± 0.0004 | 429.0 | 196.9 | +0.0009 ± 0.0004 | +0.006 ± 0.010 |

## Symmetry check: rapidity-even v1 relative to the reaction plane (10–70 %)

Must vanish for Au+Au with the true reaction plane; STAR (spectator plane) and the slopes see only the odd part.

| arm | π | K | p | p̄ |
|---|---|---|---|---|
| dm25_main | +0.0001 ± 0.0001 | -0.0000 ± 0.0001 | +0.0006 ± 0.0003 | +0.0000 ± 0.0004 |
| dm25_pub | -0.0010 ± 0.0002 | -0.0011 ± 0.0002 | +0.0017 ± 0.0005 | +0.0150 ± 0.0006 |

## Elliptic flow v2(pT), |y| < 1 (χ²/N against interpolated data; model estimator: scalar product in each experiment's geometry)

| arm | K vs STAR K0S 10–40 | K vs STAR K0S 0–80 | π vs PHENIX 10–40 (pT>1) | K vs PHENIX 10–40 | p vs PHENIX 10–40 |
|---|---|---|---|---|---|
| dm25_main | 388.0/8 | 316.3/8 | 124.6/5 | 111.3/5 | 194.2/5 |
| dm25_pub | 259.7/8 | 198.5/8 | 59.7/5 | 46.5/5 | 26.4/5 |

## Bulk chemistry at |y| < 0.5 (analysis v3 yields)

STAR K−/π− (PRC 79 034909, HEPData ins793126 Fig. 32): 0.152 ± 0.019 (N_part 299–351), 0.142 ± 0.014 (N_part 78).

| arm | Nch(\|η\|<0.5) 0–10 % | K−/π− 0–10 % | K−/π− 40–80 % | p̄/p 0–10 % | p/π+ 0–10 % |
|---|---|---|---|---|---|
| dm25_main | 538 | 0.161 | 0.152 | 0.587 | 0.054 |
| dm25_pub | 542 | 0.159 | 0.146 | 0.654 | 0.093 |

## Parton-phase end (src_v9 Z record)

| arm | events with Z | partons frozen only by the HQFRZG guard |
|---|---|---|
| dm25_main | 2394600 | 5390083 of 7707250194 (0.070 %) |
