# From the reaction plane to the event plane: how model v2 is compared with data

Written 2026-10-22. It describes the flow estimators of `analysis/analyze_arm.py` (introduced in its versions v5–v6; all paper results use them). Arm names such as dm8_p3r3 and dm15_aa are development samples of the same framework and are not shipped.

## 1. The problem

A transport model knows the true reaction plane (RP). In AMPT, the impact parameter lies along x, so Ψ_RP = 0.
An experiment never sees Ψ_RP. It reconstructs an event plane (EP) from the azimuths of tracks in its detector and corrects
for the finite resolution of that plane. These two procedures do not measure the same quantity:

- **Reaction plane:** v2{RP} = ⟨cos 2(φ − Ψ_RP)⟩. This is the mean flow projected onto a fixed geometric axis.
- **Event plane:** v2{EP} = ⟨cos 2(φ − Ψ_EP)⟩ / R. With an event-by-event fluctuating initial state, the event plane follows each
  event's own flow (participant) plane, and fluctuations and non-flow enter with a weight that depends on the
  resolution.

Comparing v2{RP} from a model with v2{EP} from data is therefore not like-for-like. In AMPT-SM Au+Au at 200 GeV,
κ = v2{EP}/v2{RP} is 1.25–1.41 at 10–40 %, about 1.5 for 0–80 % and about 3 at 0–5 % (dm15_aa, §4).

## 2. What the literature prescribes

| Point | Reference | Consequence for us |
|---|---|---|
| EP method: Q-vector Ψ_n = atan2(ΣwSin nφ, ΣwCos nφ)/n, sub-event resolution R_sub = √⟨cos n(Ψ_A − Ψ_B)⟩, full-event resolution by χ extrapolation, autocorrelation removed by excluding the POI from its own plane | Poskanzer & Voloshin, PRC 58 1671 (1998), [nucl-ex/9805001](https://arxiv.org/abs/nucl-ex/9805001) | Same construction in the model |
| v{EP} lies between ⟨v⟩ (high resolution) and √⟨v²⟩ (low resolution): v{subEP}² = ⟨v⟩² + (α − 1)σ_v², α = 2 − 4i₁²/(i₀+i₁)², i = I(χ_s²/2); non-flow δ enters with a similar weight | Ollitrault, Poskanzer & Voloshin, PRC 80 014904 (2009), [0904.2315](https://arxiv.org/abs/0904.2315), Eqs. 22–23, 31–32 | The model must reproduce the experiment's resolution, or show that both sit in the same limit |
| The EP result depends on the detector's resolution and acceptance. Use the scalar product v{SP} = ⟨Q u_A*⟩/√⟨Q_A Q_B*⟩, which always gives √⟨v²⟩. At RHIC, v2{EP} and v2{SP} differ by about 5 % or less | Luzum & Ollitrault, PRC 87 044907 (2013), [1209.2323](https://arxiv.org/abs/1209.2323) | Compute SP as a resolution-free cross-check |
| AMPT study: correct the resolution per event (or per narrow centrality bin) before averaging a wide class, i.e. ⟨v_obs/R⟩ rather than ⟨v_obs⟩/⟨R⟩; detector efficiency enters only through the multiplicity, and so the resolution | Nasim & Mohanty, IJMPE 24 1550027 (2015), [1402.6093](https://arxiv.org/abs/1402.6093) | 5 % bins inside every wide class |
| Fluctuations relative to the RP: v_PP² ≈ v_RP² + σ_v², v{2}² = ⟨v⟩² + σ_v² + δ, v{4} ≈ v_RP (Gaussian ε fluctuations) | Ollitrault, Poskanzer & Voloshin, Eqs. 13, 15, 49–50 | Explains κ > 1 even without non-flow |

**The experiments' own procedures.**

- **STAR D0 v2** (PRL 118 212301, [1701.06060](https://arxiv.org/abs/1701.06060)):
  - Ψ_2 is built from TPC tracks in the η hemisphere opposite the D0, outside |Δη| > 0.05.
  - The detector's azimuthal non-uniformity is corrected.
  - D0 yields are weighted by 1/R for each centrality interval.
  - A two-particle correlation method with the same η gap gives consistent results. That method sits at the √⟨v²⟩ limit.
  - Non-flow is estimated from p+p D–h correlations and shown as a separate one-sided band.
- **STAR K0S** (PRC 77 054901): v2{EP} from the TPC event plane.
- **PHENIX π/K/p** (PRC 85 064914, ppg123): reaction-plane detector (RXN) at 1.0 < |η| < 2.8, with a large gap.

## 3. Our implementation: experimental step and model counterpart

| Experimental step | Model (analyze_arm v5–v6) | Status |
|---|---|---|
| Centrality from mid-rapidity multiplicity | Nch(\|η\| < 0.5) percentile of the arm | as before |
| Ψ_2 from charged tracks, 0.2 < pT < 2 GeV, pT-weighted, \|η\| < 1, opposite hemisphere, \|η\| > 0.05 | Same cuts; plane from the reduced file's hadrons. Old files: all K±, p, p̄ plus 25 % of pions, prescale-weighted. reduce_dmeson_v2 files: every charged hadron (Q line) | done |
| PHENIX RXN plane 1 < \|η\| < 2.8 | N and S planes in that window, unweighted; POI \|η\| < 0.35 | done |
| Autocorrelation removal | The POI is never in its own plane (opposite hemisphere / gap) | done |
| Acceptance flattening and re-centring | Not needed: the model acceptance is azimuthally uniform | n/a |
| Resolution R_sub = √⟨cos 2(Ψ_A − Ψ_B)⟩ | Same, per 5 % centrality bin; the correction is applied per bin before a wide class is combined (Nasim–Mohanty) | done |
| Tracking efficiency (~85 %) | Affects only the multiplicity, hence χ_s. Both model and data sit near the rms limit (§4), so the effect is ≲ 1–2 % | quantified |
| Two-particle / correlation cross-check | Scalar product, same geometry: v2SPS_* (STAR), v2SPF_* (PHENIX) | done (v6) |
| Non-flow | The model's own non-flow at the same gap is included. STAR's one-sided band is kept as a systematic in the Table I χ². Model small gap (S) vs large gap (F) bounds AMPT's non-flow | done |
| N_coll weighting of D0 events | Each triggered event carries its exact NCOLT | as before |

The v1 comparison stays relative to the RP. STAR measures v1 against the first-order spectator (ZDC-SMD) plane, which
tracks the reaction plane and is free of mid-rapidity fluctuations and non-flow.

## 4. Verification

- All earlier output keys are byte-identical after v5 and v6 (dm8_p3r3, dm15_aa).
- Split/merge is byte-identical. Batch-farm and local runs agree to 3×10⁻¹³.
- For the D0, the direct EP and the factorised κ·v2{RP} agree: χ² 20.5/24 (dm8_p3r3).
- κ is species-independent (π, K, p within 5 %) and survives a Δη ≈ 1 gap (S ≈ F). So in AMPT it is fluctuation-driven, not non-flow.
- OPV α from the model resolution: R_sub = 0.42–0.50 at 5–50 % gives α = 1.97–1.98. A detector with R_sub = 0.65 would give
  α = 1.90. The resolution-mismatch bias on v2 is therefore ≈ (2 − α)σ_v²/(2v²) ≈ 1 % or less: the model and STAR both
  sit at the rms limit, where EP ≈ SP.
- **Direct check of EP against SP (dm8_p3r3, 89k events, analysis v6), pT 0.5–2.5 GeV:**

  | geometry | class | π SP/EP | K SP/EP | p SP/EP | SP/RP (π, K, p) |
  |---|---|---|---|---|---|
  | STAR η-sub | 10–40 % | 1.002 | 1.003 | 1.028 | 1.29, 1.32, 1.25 |
  | STAR η-sub | 0–80 % | 1.003 | 1.007 | 1.028 | 1.49, 1.50, 1.42 |
  | STAR η-sub | 40–80 % | 0.990 | 0.989 | 0.996 | 1.78, 1.65, 1.39 |
  | PHENIX 1<\|η\|<2.8 | 10–40 % | 1.008 | 0.986 | 1.018 | 1.27, 1.29, 1.35 |
  | PHENIX 1<\|η\|<2.8 | 0–80 % | 0.81 | 0.69 | 0.97 | 1.44, 1.53, 1.52 |

  In the STAR geometry EP and SP agree to ≤ 1 % for π and K and ≤ 3 % for protons in every class, as α ≈ 1.97–1.98
  predicts. In the PHENIX geometry they agree at 10–40 % (the class compared with PHENIX), but the old files' forward
  plane (10 % prescaled hadrons) has R_sub → 0–0.15 beyond 50 %, where the EP division is unstable; SP is not.
- **Decision:** the primary model estimator is the scalar product in each experiment's geometry (v2SPS_*, v2SPF_*),
  as recommended by Luzum & Ollitrault; it equals STAR's EP here and STAR's own two-particle cross-check, and its D0
  errors are ~10 % smaller than EP's. v2EPS_* is kept as the exact replica of STAR's procedure.

## 5. Remaining limits

- The D0 v2{EP} model errors are about 2.2/√2 × those of v2{RP} (resolution ≈ 0.45). Above 3 GeV they exceed STAR's at 300k events.
  Full-multiplicity Q vectors (reduce_dmeson_v2) and more events are needed. Compare arms of different size only with the
  noise-subtracted χ²∞ = Σ[(m − d)² − σ_m²]/σ_d².
- AMPT's non-flow (minijets, resonance decays) is not tuned to data. Its size at the STAR gap is bounded by S vs F.
