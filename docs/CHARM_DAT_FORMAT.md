# `ana/charm.dat` format (AMPT_CHARM_OUT=1)

The file is whitespace-separated, and each line starts with a record type.
Events appear in order.

- **Repeated events.** If HIJING repeats an event (for example NATT = 0),
  a new **E** block with the same event number is written. Keep the
  **last** block per event number, because it is the event that reaches
  `ampt.dat` (`analysis/dmlib.py: read_charm`).
- **Units.** GeV for momenta and masses; fm and fm/c for positions and
  times. The frame is the NN CMS, and the reaction plane is +x.

| Type | When written | Fields after the type letter |
|---|---|---|
| E | end of HTOP (before ZPC) | iev, b (fm), NP, NT (participants), NATT (hadrons before melting), npar (partons into ZPC), nc0 (c + c̄ partons), nhqs (charm saved from strings), nhqj1 (charm on projectile/target strings), nhqovf (save overflow, must be 0), nmatch, nfail (HQDIRECT), dE (GeV, HQDIRECT energy change), σ_trig = HINT1(61) (mb) |
| I | end of HTOP | iev, ityp (±4), px, py, pz, m, x, y, z, t_form, string id |
| L | end of HTOP, sampled with probability AMPT_LQFRAC | same as I, for u/d partons (ityp ±1, ±2) |
| F | start of PTOH (ZPC finished) | iev, ityp, px, py, pz, m, x, y, z, t (last interaction or formation point), ncol (ZPC scatterings of this c) |
| G | start of PTOH, sampled | same as F, for u/d partons (ncol column = 0) |
| H | end of PTOH (right after coalescence) | iev, KF, px, py, pz, m, x, y, z, t (hadron formation), pcx, pcy, pcz (charm constituent momentum), dmax (largest constituent separation in fm at the latest constituent freeze-out time, lab frame), nq (2 or 3) |
| Z | start of PTOH, only with ISOFT=5 or `AMPT_ZPC_TSTOP` | iev, ISOFT, partons in ZPC, partons frozen only by the HQFRZG guard (must be ≈ 0), ε_c (GeV/fm³, −1 = unset), lab stop time (fm/c, −1 = unset) |

**Initial positions.** In I and L, z = 0 and t = t_form. ZPC
free-streams a parton from (x, y, 0, 0) to its formation point, so the
spatial position at formation is x + (px/E) t_form, with
η_s = y (momentum rapidity). `analysis/analyze_arm.py` uses this to compute ε1.

**Reduced `ampt.dat`** (`red_*.dat.gz`, written by `scripts/reduce_dmeson_v2.awk`), one event at a
time:

- **Header:** the 11 stock header fields, then Nch(\|η\| < 0.5) and
  Nch(\|η\| < 1), 13 fields in all. Field 4 is b.
- **Q line** (right after the header): `Q nch A_w A_x A_y B_w B_x B_y N_w N_x N_y S_w S_x S_y`,
  14 fields: second-order Q vectors (Σw, Σw cos 2φ, Σw sin 2φ) of all charged hadrons with 0.2 < pT < 2 GeV,
  without prescale. A: −1 < η < −0.05 and B: 0.05 < η < 1, pT-weighted (STAR TPC η-sub planes);
  N: 1 < η < 2.8 and S: −2.8 < η < −1, unweighted (PHENIX RXN-like planes). nch counts the hadrons used.
- **Charm hadrons:** the stock line verbatim, 9 fields:
  KF px py pz m x y z t.
- **Light hadrons:** KF px py pz m w, 6 fields, where w is the inverse
  prescale:

  | Species | Rapidity window | Prescale | w |
  |---|---|---|---|
  | K±, p, p̄ | \|y\| < 1.5 | none | 1 |
  | π± | \|y\| < 1.5 | 0.25 | 4 |
  | π±, K± | 1.5 ≤ \|y\| < 4 | 0.10 | 10 |

  The prescale uses awk's `rand()`, seeded per job (`srand(seed)`), and never touches the AMPT random
  sequence. Only gawk reproduces the paper files' choice exactly (README Sec. 3).

**Other files of a job** (`scripts/run_job.sh`): `charm_*.dat.gz` is `ana/charm.dat` above; `hq_*.log.gz` is the
AMPT log, which with `AMPT_HQSHW=1` carries one `# hqshw` line per event (the two incoming partons' x, flavour,
Q², transverse distances from the nuclear centres and the logged shadowing factors), used by `analysis/eps09s/eps09s_errw.py`.
