# Claims audit: every entity on the same footing (ana/audit_claims.py)

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

## Sums

| entity | 46-point scoreboard | Ds/D0 terms | all |
|---|---|---|---|
| this work | 91.5 ± 7.4 (46) | 196.0 ± 5.3 (18) | 287.4 ± 9.1 (64) |
| public AMPT | 173.5 ± 22.6 (46) | 95.3 ± 9.9 (18) | 268.7 ± 24.7 (64) |

## D0 v2 10-40 %

| entity | χ² (ours/public: χ²∞ ± err) | ndf | note |
|---|---|---|---|
| this work | 20.5 ± 6.3 | 8 | v2{SP}, STAR eta-sub geometry |
| public AMPT | 64.0 ± 21.6 | 8 | v2{SP}, STAR eta-sub geometry |

## D0 v2 10-40 % (non-flow incl.)

| entity | χ² (ours/public: χ²∞ ± err) | ndf | note |
|---|---|---|---|
| this work | 9.9 ± 3.6 | 8 | sensitivity |
| public AMPT | 29.5 ± 11.5 | 8 | sensitivity |

## D0 v2 0-80 %

| entity | χ² (ours/public: χ²∞ ± err) | ndf | note |
|---|---|---|---|
| this work | 14.5 ± 3.5 | 8 | v2{SP}, STAR eta-sub geometry |
| public AMPT | 15.8 ± 6.4 | 8 | v2{SP}, STAR eta-sub geometry |
| TAMU_cdiff [pt] | 8.0 | 7 | STAR Table I: 10.0/8 |
| TAMU_cdiff [bin] | 8.2 | 7 | STAR Table I: 10.0/8 |
| TAMU_nocdiff [pt] | 25.8 | 7 | STAR Table I: 29.5/8 |
| TAMU_nocdiff [bin] | 26.1 | 7 | STAR Table I: 29.5/8 |
| SUBATECH [pt] | 15.2 | 8 | STAR Table I: 15.2/8 |
| SUBATECH [bin] | 14.0 | 7 | STAR Table I: 15.2/8 |
| hydro3D [pt] | 3.6 | 6 | STAR Table I: 3.6/6 |
| hydro3D [bin] | 3.6 | 6 | STAR Table I: 3.6/6 |
| Duke [pt] | 35.8 | 8 | STAR Table I: 35.7/8 |
| Duke [bin] | 33.1 | 7 | STAR Table I: 35.7/8 |
| PHSD [pt] | 8.5 | 7 | STAR Table I: 8.7/7 |
| PHSD [bin] | 8.0 | 6 | STAR Table I: 8.7/7 |
| LBT [pt] | 11.0 | 8 | STAR Table I: 11.1/8 |
| LBT [bin] | 9.0 | 7 | STAR Table I: 11.1/8 |

## D0 R_AA 0-10 %

| entity | χ² (ours/public: χ²∞ ± err) | ndf | note |
|---|---|---|---|
| this work | 9.7 ± 1.3 | 8 | per-c-quark yields AA/pp |
| public AMPT | 15.8 ± 1.0 | 8 | per-c-quark yields AA/pp |
| Fan et al. 2026 (improved AMPT) [pt] | 5.5 | 8 |  |
| Fan et al. 2026 (improved AMPT) [bin] | 5.5 | 8 |  |
| LBT [pt] | 26.2 | 8 |  |
| LBT [bin] | 25.4 | 8 |  |
| Duke [pt] | 15.8 | 8 |  |
| Duke [bin] | 15.5 | 8 |  |

## D0 R_AA 10-40 %

| entity | χ² (ours/public: χ²∞ ± err) | ndf | note |
|---|---|---|---|
| this work | 10.4 ± 0.3 | 8 | per-c-quark yields AA/pp |
| public AMPT | 17.9 ± 0.8 | 8 | per-c-quark yields AA/pp |
| Fan et al. 2026 (improved AMPT) [pt] | 2.8 | 8 |  |
| Fan et al. 2026 (improved AMPT) [bin] | 2.8 | 8 |  |

## D0 R_AA 40-80 %

| entity | χ² (ours/public: χ²∞ ± err) | ndf | note |
|---|---|---|---|
| this work | 30.8 ± 0.5 | 8 | per-c-quark yields AA/pp |
| public AMPT | 29.0 ± 0.8 | 8 | per-c-quark yields AA/pp |
| Fan et al. 2026 (improved AMPT) [pt] | 0.9 | 8 |  |
| Fan et al. 2026 (improved AMPT) [bin] | 0.9 | 8 |  |

## Λc/D0 vs pT (10-80 %)

| entity | χ² (ours/public: χ²∞ ± err) | ndf | note |
|---|---|---|---|
| this work | 3.1 ± 0.2 | 3 |  |
| public AMPT | 15.2 ± 0.9 | 3 |  |
| Fan et al. 2026 (improved AMPT) [pt] | 1.7 | 3 |  |
| Fan et al. 2026 (improved AMPT) [bin] | 1.6 | 3 |  |
| Ko et al.: three quark (0-5 %) [pt] | 1.4 | 3 |  |
| Ko et al.: three quark (0-5 %) [bin] | 1.5 | 3 | * |
| Ko et al.: di-quark (0-5 %) [pt] | 0.7 | 3 |  |
| Ko et al.: di-quark (0-5 %) [bin] | 0.7 | 3 | * |
| Ko et al.: with flow (0-10 %) [pt] | 8.1 | 3 |  |
| Ko et al.: with flow (0-10 %) [bin] | 6.5 | 3 | * |
| Catania, coal.+frag. (10-80 %) [pt] | 1.8 | 3 |  |
| Catania, coal.+frag. (10-80 %) [bin] | 3.1 | 3 |  |
| Catania, coal. (10-80 %) [pt] | 66.2 | 3 |  |
| Catania, coal. (10-80 %) [bin] | 50.1 | 3 |  |
| Tsinghua (10-80 %) [pt] | 1.0 | 3 |  |
| Tsinghua (10-80 %) [bin] | 0.6 | 3 |  |
| Rapp et al. (0-20 %) [pt] | 4.1 | 3 |  |
| Rapp et al. (0-20 %) [bin] | 2.7 | 3 |  |
| Cao et al. (10-80 %) [pt] | 0.7 | 3 |  |
| Cao et al. (10-80 %) [bin] | 0.9 | 3 | * |
| PYTHIA 8 [pt] | 21.1 | 3 |  |
| PYTHIA 8 [bin] | 20.8 | 3 |  |
| PYTHIA 8, CR [pt] | 7.7 | 3 |  |
| PYTHIA 8, CR [bin] | 8.5 | 3 |  |

## Λc/D0 vs centrality (3-6 GeV)

| entity | χ² (ours/public: χ²∞ ± err) | ndf | note |
|---|---|---|---|
| this work | 2.5 ± 0.2 | 3 |  |
| public AMPT | 15.6 ± 0.5 | 3 |  |
| Catania (coal.) | 75.1 | 3 | curve drawn at the 3 STAR Npart |
| Catania (coal.+frag.) | 0.7 | 3 | curve drawn at the 3 STAR Npart |

## Ds/D0 vs pT (0-10 %)

| entity | χ² (ours/public: χ²∞ ± err) | ndf | note |
|---|---|---|---|
| this work | 17.7 ± 3.3 | 4 | error ignores Ds-D0 correlation (upper bound) |
| public AMPT | 19.9 ± 7.0 | 4 | error ignores Ds-D0 correlation (upper bound) |
| Tsinghua (seq. coal.) | 26.0 | 4 | curve |
| Catania (coal.) | 19.3 | 4 | curve |
| Catania (coal.+frag.) | 20.5 | 4 | curve |
| Cao, Ko | 36.0 | 4 | curve |
| He, Rapp (0-20 %) | 16.4 | 4 | curve |

## Ds/D0 vs pT (10-40 %)

| entity | χ² (ours/public: χ²∞ ± err) | ndf | note |
|---|---|---|---|
| this work | 68.9 ± 3.6 | 5 | error ignores Ds-D0 correlation (upper bound) |
| public AMPT | 31.8 ± 6.3 | 5 | error ignores Ds-D0 correlation (upper bound) |
| Tsinghua (seq. coal.), 20-40 % [20-40 % curve] | 8.0 | 5 | curve |

## Ds/D0 vs pT (40-80 %)

| entity | χ² (ours/public: χ²∞ ± err) | ndf | note |
|---|---|---|---|
| this work | 37.9 ± 2.0 | 4 | error ignores Ds-D0 correlation (upper bound) |
| public AMPT | 17.4 ± 2.8 | 4 | error ignores Ds-D0 correlation (upper bound) |
| Tsinghua (seq. coal.), 40-80 % | 33.7 | 4 | curve |

## Ds/D0 vs centrality (1.5-5 GeV)

| entity | χ² (ours/public: χ²∞ ± err) | ndf | note |
|---|---|---|---|
| this work | 71.5 ± 0.8 | 5 |  |
| public AMPT | 26.2 ± 1.2 | 5 |  |
| Tsinghua (seq. coal.) | 42.6 | 5 | curve |

## D0+D0bar dv1/dy, 10-80 %, pT > 1.5 GeV/c

STAR points refitted through the origin with stat weights: -0.0803 ± 0.0174 (published −0.080 ± 0.017 stat; total with sys 0.023).

| entity | slope [at y] | slope [bin avg] | (STAR − model)/σ | note |
|---|---|---|---|---|
| this work | -0.0306 | -0.0306 | +2.1 | production slope (model weights) -0.0312 ± 0.0052 |
| public AMPT | +0.0080 | +0.0080 | +3.5 | production slope (model weights) +0.0061 ± 0.0100 |
| hydro+EM (Chatterjee–Bożek), STAR Fig. 3 | -0.0200 | -0.0200 | +2.6 | curve |
| AMPT (Nasim–Singha), STAR Fig. 3 | -0.0147 | -0.0150 | +2.8 | curve |
| Catania (Oliva et al. 2021), eta_m 1.1 | -0.0574 | -0.0577 | +1.0 | curve, b = 9 fm |
| Catania (Oliva et al. 2021), eta_m 1.2 | -0.0530 | -0.0532 | +1.2 | curve, b = 9 fm |
| Catania (Oliva et al. 2021), charm tilted like the bulk | -0.0018 | -0.0017 | +3.4 | curve, b = 9 fm |

