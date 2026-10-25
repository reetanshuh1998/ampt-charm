# Claims audit before the paper draft (2026-10-25)

Question: for every observable where the paper would say "this work describes the data better than model X", is
the number right, and is the comparison like for like (same estimator, kinematic window, centrality, normalisation
and χ² definition)? This is the same exercise as the reaction-plane vs event-plane check (docs/FLOW_METHOD.md).

Tools:
- analysis/indep_check.py and analysis/indep_compare.py: independent re-analysis, raw events.
- analysis/audit_claims.py: independent χ² of every entity, with output in observation/results/tables/AUDIT_TABLES.md.
- data/models/extract_*.py: digitised model curves, each calibrated against STAR's own markers.

## 1. Is the analysis code right? (independent re-analysis on raw events)

**Sample.**
- 121 Au+Au files (36,300 events), spread over the whole dm25_main arm (every 66th file).
- 10 p+p files (300,000 events) of dm25_ppmain.

**The independent code (analysis/indep_check.py).**
- It was written from the format and method documentation, without importing dmlib.py or analyze_arm.py.
- It re-implements:
  - parsing;
  - the Nch-percentile centrality;
  - N_coll × nPDF event weights;
  - the per-c-quark normalisation;
  - D*, Ds*, Σc feed-down;
  - the v1 bins and the fit through the origin;
  - v2{RP};
  - v2{SP} in STAR's η-sub geometry, with the resolution taken per 5 % bin;
  - R_AA;
  - Λc/D0 vs pT and vs centrality.
- Feed-down draws follow the documented seed convention (one generator per file), so the same decays are drawn.

**Result.** Every observable agrees with analyze_arm.py on the same events to ≤ 6 × 10⁻¹⁸ (observation/results/tables/INDEP_CHECK.md).

**Second run, with a private generator.** All differences stay within the statistical errors (max |pull| 1.4 over 8-bin spectra).

**What this cannot catch.** A definition that is wrong in the specification itself. That is what §2 checks.

**χ² recomputation.** analysis/audit_claims.py uses its own HEPData reader and its own χ². It reproduces every term of the scoreboard exactly, for us and for public AMPT:
- D0 v2: 20.5 and 14.5;
- R_AA: 9.7, 10.4 and 30.8;
- Λc/D0: 3.1 and 2.5;
- public AMPT: 64.0, 15.8, 15.8, 17.9, 29.0, 15.2 and 15.6.

**A bug in the audit script itself.** Its first version matched the R_AA bins by index and was shifted by one bin, which gave 5.9 instead of 9.7. It now matches by bin edges, with an assertion. The production code was not affected.

## 2. Is the comparison like for like?

| Item | STAR | This work and public AMPT | Published curves | Verdict |
|---|---|---|---|---|
| D0 v2 estimator | η-sub EP, resolution per centrality bin | Scalar product, same geometry, per 5 % bin; SP ≈ EP at the rms limit, verified (docs/FLOW_METHOD.md §4) | Each group's own v2 (mostly smooth initial conditions, so RP = EP) | Each entity gives its best estimate of the measured quantity. Our χ² code reproduces STAR's Table I: SUBATECH 15.2/15.2, hydro 3.6/3.6, Duke 35.8/35.7, PHSD 8.5/8.7, LBT 11.0/11.1 |
| v2 0–80 % non-flow | One-sided band, Table I | Same | Same | Like for like |
| v2 10–40 % non-flow | Band given, not used in STAR's χ² | Not added (the scoreboard) | No curves | Including it gives ours 9.9 ± 3.6, public 29.5 ± 11.5; the ratio is unchanged |
| D0 v1 plane and sign | ZDC-SMD 1st-order spectator plane; projectile spectators v1 > 0 | True RP; HIJING puts the projectile at x = +b/2 (linana.f:3318, hijing1.383_ampt.f:784), so the projectile spectators lie at +x | As published | Sign verified directly; the spectator plane tracks the RP |
| D0 v1 window and fit | 10–80 %, pT > 1.5, 4 y bins, line through origin | Same bins and fit; the STAR refit gives −0.0803 (published −0.080) | Slope from curves at y = ±0.2, ±0.6 and bin-averaged (identical to 0.0003) | Like for like. Caveat: Oliva et al. use b = 9 fm, not 10–80 % |
| R_AA normalisation | dN_AA / (⟨N_bin⟩ dN_pp), Levy pp reference | Per-c-quark yields AA / pp. The two are equal when total charm scales with N_coll × nPDF, which the trigger + N_coll weight enforces. Au+Au events carry 12 % extra HIJING charm per event; the per-c normalisation cancels it, whereas an N_coll-only normalisation would bias R_AA by 12 % | Fan et al.: Eq. (5), the standard definition; LBT, Duke: as plotted by STAR | Like for like |
| R_AA points | 8 points, 0.9 < pT < 8 | Model bins are STAR's bins (asserted by edges) | Evaluated at ⟨pT⟩ and bin-averaged; the χ² differs by ≤ 0.8 | Like for like |
| Λc/D0 | 10–80 %, three pT bins; centrality at 3–6 GeV | Same bins and classes | Several curves are for other centralities (Ko 0–5 %, Rapp 0–20 %), as STAR notes | Our χ² reproduces STAR's quoted PYTHIA 20.7/3 (ours 20.8–21.1) and PYTHIA CR 8.2/3 (7.7–8.6). STAR's Catania (reduced 26.1) and Ko-with-flow (4.8) numbers use a finite-bin-width treatment we cannot replicate |
| χ² for transport output vs curves | — | χ²∞ (noise-subtracted), whose expectation is the χ² of the infinite-statistics model | Plain χ² (no model error) | Same quantity. Using an ordinary χ² with model error would favour small samples, which is why the bulk-v2 χ² table below is not used for ranking |

Digitisation:
- Λc/D0 vs pT: calibrated on the three STAR markers, residual 0.0011.
- LBT and Duke R_AA: eleven markers, residuals 0.005 in pT and 0.0007 in R_AA.
- Oliva et al. D0 v1: four markers per panel, residual 3 × 10⁻⁵.
- Fan et al.: earlier extraction, reproduced to 0.01. Their published curve has only 6 points, so it is interpolated linearly, as drawn.

## 3. What the numbers support (χ²; ours and public = χ²∞ ± stat)

| Observable (points) | This work | Public AMPT | Published models | Supportable statement |
|---|---|---|---|---|
| D0 v2 10–40 % (8) | 20.5 ± 6.3 | 64.0 ± 21.6 | none drawn | Better than public AMPT (2σ on the difference) |
| D0 v2 0–80 % (8) | 14.5 ± 3.5 | 15.8 ± 6.4 | hydro 3.6, TAMU 8.0, PHSD 8.0–8.5, LBT 9.0–11.0, SUBATECH 14.0–15.2, TAMU no-diff 26, Duke 33–36 | Comparable to public AMPT and SUBATECH; worse than hydro, TAMU, PHSD and LBT |
| R_AA 0–10 % (8) | 9.7 ± 1.3 | 15.8 ± 1.0 | Fan 5.5, Duke 15.5–15.8, LBT 25.4–26.2 | Better than public AMPT, Duke and LBT; worse than Fan |
| R_AA 10–40 % (8) | 10.4 ± 0.3 | 17.9 ± 0.8 | Fan 2.8 | Better than public AMPT; worse than Fan |
| R_AA 40–80 % (8) | 30.8 ± 0.5 | 29.0 ± 0.8 | Fan 0.9 | No better than public AMPT; much worse than Fan (the coalescence bump) |
| Λc/D0 vs pT (3) | 3.1 ± 0.2 | 15.2 ± 0.9 | Tsinghua 0.6–1.0, Ko di-quark 0.7, Cao 0.7–0.9, Ko three-quark 1.4, Fan 1.6, Catania c+f 1.8–3.1, Rapp 2.7–4.1, Ko flow 6.5–8.1, PYTHIA CR 7.7–8.6, PYTHIA 21, Catania coal. 50–66 | Far better than public AMPT; inside the band of coalescence models but at its upper end; worse than Fan |
| Λc/D0 vs centrality (3) | 2.5 ± 0.2 | 15.6 ± 0.5 | Catania c+f 0.7, Catania coal. 75 | Better than public AMPT; worse than Catania c+f. Our centrality trend is flat (0.60, 0.70, 0.64) where STAR rises |
| **Ds/D0 vs pT and centrality (18)** | **196.0 ± 5.3** | **95.3 ± 9.9** | 0–10 %: Rapp 16.4, Catania 19–21, Tsinghua 26, Cao–Ko 36; ours 17.7 | **Public AMPT is better by a factor of 2** (its old coalescence gives Ds/D0 ≈ 0.28; our hybrid ≈ 0.20). Every model misses the peripheral enhancement |
| D0 dv1/dy (STAR −0.080 ± 0.023) | −0.031 ± 0.005 (2.1σ; 2.0σ with the EPS09s band) | +0.006 ± 0.010 (3.5σ) | Catania (Oliva) −0.053 to −0.058 (1.0–1.2σ); hydro+EM −0.020 (2.6σ); AMPT Nasim–Singha −0.015 (2.8σ) | Better than public AMPT, the old AMPT and the hydro+EM curve shown by STAR. **Worse than Catania**, which fits a bulk tilt to charged-hadron v1 and keeps charm untilted |
| p+p D0 shape (5) | 3.5 | 6.7 | — | Both acceptable within STAR's p+p errors; "comparable or better", not decisive |
| Bulk π, K v2 | — | — | — | Ours and public agree to 3–5 %. Public has stronger proton mass ordering below 1.5 GeV (closer to PHENIX) |
| Bulk v1 slopes | π −0.0015, K +0.0021, p +0.0055 | π +0.0019, K +0.0042, p +0.0034 | STAR −0.0027 / −0.0017 / −0.0015 | Ours has the right sign for π (public does not); both have the wrong sign for K and p |

**Sums.**
- 46-point scoreboard (v2, R_AA, Λc/D0): this work 91.5 ± 7.4 vs public AMPT 173.5 ± 22.6.
- Ds/D0 (18 points): 196.0 ± 5.3 vs 95.3 ± 9.9.
- All 64 points: 287.4 ± 9.1 vs 268.7 ± 24.7. **Over all 64 points the two are statistically equal.**

## 4. Consequences for the paper's claims

**Keep (verified, like for like):**
- Against public AMPT, better on D0 v2 (10–40 %), R_AA (0–10, 10–40 %), Λc/D0 (pT and centrality) and D0 v1, and on the 46-point set by a factor of 1.9. State the observable set explicitly.
- D0 dv1/dy −0.031 ± 0.005: the largest of the AMPT calculations (Nasim–Singha −0.015, public +0.006) and larger than the hydro+EM curve shown by STAR. It comes with no tilt parameter.
- The tilt law: v1(c) is linear in ε1(c)/ε1(u) across shadowing, PDFs, r_BM and public AMPT, which falls on the line. The Catania limits (charm untilted −0.057, charm tilted like the bulk −0.002) are an independent confirmation of the same mechanism in a different transport model.
- The one-framework scope: no other single calculation covers v1, v2, R_AA in three centralities, Λc/D0 vs pT and centrality, Ds/D0 and hadrochemistry. This is a scope statement, not a "better fit" statement.

**Drop or qualify:**
- "Half the χ² of public AMPT" without qualification. It holds only on the 46-point set. With Ds/D0 the two models tie (287 ± 9 vs 269 ± 25). State both.
- "Best D0 v1 of any model". False: Catania (Oliva et al.) is closer, with a fitted bulk tilt. Our claim is the mechanism, not the value.
- "Ds/D0 is the AMPT-SM strangeness limit". Incomplete: public AMPT, with the same γ_s, gives 0.28. The hybrid coalescence/fragmentation lowers Ds/D0 to the fragmentation value (0.21). This is a limitation of our hadronisation, to be stated.
- Any statement that we beat the published models on R_AA or Λc/D0. Fan et al. (R_AA, Λc/D0 vs pT) and several coalescence models (Λc/D0) do better.
- p+p shape: say "comparable to or better than", not "better".
