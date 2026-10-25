# ampt-charm: open-charm hadrons in AMPT at RHIC

An extension of the public **AMPT** model (v1.26t9b/v2.26t9b, string melting) for open-charm hadrons in Au+Au and
p+p collisions at √s<sub>NN</sub> = 200 GeV. It contains everything needed to rebuild the program, regenerate the
events, and reproduce every number and figure of the accompanying paper (in preparation). That includes D⁰ directed
and elliptic flow, R<sub>AA</sub>, Λc/D⁰, Ds/D⁰ and the charm hadrochemistry, compared with STAR, public AMPT and other
transport models. **No batch farm is needed**: everything runs on a Linux laptop or in Docker.

- **What was changed in AMPT, where, and why:** [`docs/CHANGES.pdf`](docs/CHANGES.pdf). It covers every switch, routine
  and patched line, with the motivating papers.
- **The change set:** [`patches/ampt-v1.26t9b-v2.26t9b_charm.patch`](patches/ampt-v1.26t9b-v2.26t9b_charm.patch)
  (+170/−17 lines in seven stock files) and four new files in [`src_new/`](src_new/). The official AMPT source is
  downloaded from its authors' page and checksum-verified; it is not redistributed here.
- **Results:** [`observation/results/`](observation/results/) holds the analysis outputs of the full production
  (2.4M Au+Au, 11.9M p+p events), the χ² tables against STAR for every model, and the figures.

## 1. Quick start

### Option A: Docker (any OS, nothing else to install)
```bash
git clone https://github.com/<owner>/ampt-charm.git && cd ampt-charm
docker build -t ampt-charm .          # downloads + verifies dependencies, builds LHAPDF and AMPT (~5 min)
docker run --rm ampt-charm            # runs the three test suites (~15-25 min)
```

### Option B: native Linux
Requirements: `gfortran`, `g++`, `make`, `curl`, `unzip`, `patch`, `python3` with `numpy` and `matplotlib`. `gawk`
is recommended (Sec. 3).
```bash
sudo apt install build-essential gfortran g++ make gawk curl unzip patch python3-numpy python3-matplotlib   # Ubuntu/Debian
git clone https://github.com/<owner>/ampt-charm.git && cd ampt-charm
make deps      # download + SHA-256-check AMPT, LHAPDF 6.5.1, CTEQ6.1M, EPS09s; build LHAPDF into deps/ (~3 min)
make build     # patch the official AMPT and compile bin/ampt
make test      # test-results + test-repro + test-quick (Sec. 3)
```
Fedora: `sudo dnf install gcc-gfortran gcc-c++ make gawk curl unzip patch python3-numpy python3-matplotlib`.
`scripts/check_prereqs.sh` lists anything missing. macOS with Homebrew `gcc` should work but is untested; use Docker
there.

Tested: Ubuntu 24.04 (gfortran 13.3) natively and in Docker. The paper production ran on Alma/RHEL 9 (gfortran
11.5); both give byte-identical events (Sec. 3).

## 2. Running your own production

One job = one input card + one set of switches + one pair of seeds, exactly as on the batch farm:
```bash
# scripts/run_job.sh CONFIG ARM SYS BMIN BMAX NEV SEEDBASE SEEDPMUL SEEDPOFF IDX [OUTDIR]
scripts/run_job.sh reference_auau myrun auau 0 14.5 50 20270001 17 3 1
# many jobs in parallel (jobs 1..40, all cores):
scripts/run_arm.sh reference_auau myrun auau 0 14.5 300 20270001 17 3 1 40
scripts/run_arm.sh reference_pp   mypp  pp   0 2.0 30000 20270002 17 5 1 4
```
- HIJING seed of job *i* = SEEDBASE + 7919·*i*; ZPC seed = SEEDPMUL·*i* + SEEDPOFF.
- Output goes to `prod/ARM/`:
  - `red_*.dat.gz`: the reduced event file (every charm hadron, all K± and p(p̄), a weighted pion sample, and the
    flow vectors);
  - `charm_*.dat.gz`: charm records at every stage;
  - `hq_*.log.gz`.

Configurations ([`configs/`](configs/)):

| config | what | paper arms |
|---|---|---|
| `reference_auau`, `reference_pp` | the framework (paper reference) | dm24_e_a/b + dm25_e_c..h = **dm25_main**; dm24_epp + dm25_epp = **dm25_ppmain** |
| `c1_auau`, `c1_pp` | systematic variant: Jacobi baryon distance | dm25_c1, dm25_ppc1 |
| `public_auau`, `public_pp` | public AMPT physics + charm trigger only | dm25_pub, dm25_pppub |
| `stock` | all switches off = stock AMPT | — |

The seeds of every production arm are in `docs/CHANGES.pdf` (Sec. 4), so any job of the paper can be regenerated.
**Cost:** about 7 s per minimum-bias Au+Au event and 0.04 s per p+p event per core (AMD EPYC):
- the 2.4M-event reference Au+Au sample took about 4,700 core-hours;
- all arms of the paper took about 6,500;
- 100k Au+Au events (enough for R<sub>AA</sub> and Λc/D⁰, not for v1/v2) take about 200 core-hours.

Analysis:
```bash
python3 analysis/analyze_arm.py prod/myrun myrun.json --sys auau -j 8     # Au+Au
python3 analysis/analyze_arm.py prod/mypp  mypp.json  --sys pp   -j 8     # p+p reference
```
The output JSON has the same keys as the published ones in `observation/results/json/`. To compare with STAR, copy
your files over `dm25_main.json` / `dm25_ppmain.json` (keep a backup) and run `python3 analysis/audit_claims.py
out.md` and `python3 analysis/paper_figs.py figs/`.

## 3. Verification: what `make test` checks

| target | time | what it does | expected |
|---|---|---|---|
| `make test-results` | 30 s | regenerates every table and figure of `observation/results/` from the shipped analysis outputs | `PASS` (tables identical) |
| `make test-repro` | 2 min | regenerates the first events of Au+Au job 1 (arm dm24_e_a) and p+p job 1 (dm24_epp) of the paper production and compares them with the files written on the farm (`tests/reference/`) | `PASS` for `charm.dat`; `PASS` for the event file with gawk, `PASS*` without gawk |
| `make test-quick` | 10–20 min | a new small production (4×20 Au+Au, 2×2000 p+p events, new seeds) → production analysis → an **independent re-implementation** of the observables on the same events (must agree exactly) → comparison of charm yields and ⟨pT⟩ with the published full-statistics results (pulls < 3.5) | `PASS` |

`PASS*` means the event headers, flow vectors and every charm hadron are identical, but the 25 % pion sample in the
reduced file differs. That sample is drawn with awk's `rand()`, which only gawk reproduces. Charm results are
unaffected; light-hadron observables agree statistically.

## 4. Results in brief (details: `observation/results/tables/AUDIT_TABLES.md`)

All entities are compared on the same footing: the same estimator and kinematic window as STAR, and the same χ²
(published curves without model error; transport output with its statistical noise subtracted, χ²∞).

| observable (points) | this work | public AMPT | best published model shown |
|---|---|---|---|
| D⁰ v2 10–40 % (8) | 20.5 ± 6.3 | 64.0 ± 21.6 | — |
| D⁰ v2 0–80 % (8) | 14.5 ± 3.5 | 15.8 ± 6.4 | 3D hydro 3.6, TAMU 8.0, PHSD 8.0–8.5 |
| D⁰ R<sub>AA</sub> 0–10 / 10–40 / 40–80 % (8 each) | 9.7 / 10.4 / 30.8 | 15.8 / 17.9 / 29.0 | improved AMPT (Fan et al.) 5.5 / 2.8 / 0.9 |
| Λc/D⁰ vs pT, vs centrality (3, 3) | 3.1 / 2.5 | 15.2 / 15.6 | Tsinghua 0.6–1.0 / Catania 0.7 |
| Ds/D⁰ (18) | 196 ± 5 | 95 ± 10 | — (all models miss the peripheral trend) |
| D⁰ dv1/dy (STAR −0.080 ± 0.023) | −0.031 ± 0.005 | +0.006 ± 0.010 | Catania (Oliva et al.) −0.055 |

- **Better than public AMPT** on D⁰ v2, R<sub>AA</sub> (0–10, 10–40 %), Λc/D⁰ and D⁰ v1.
- **Worse than public AMPT** on Ds/D⁰. Summed over all 64 points the two are statistically equal (287 ± 9 vs 269 ± 25).
- **Main physics result:** the charm directed flow follows the charm-to-bulk tilt ratio,
  v1(c) ≈ −0.035 (1 − ε1(c)/ε1(u)). It holds across nPDF and PDF variants and in public AMPT
  (`observation/results/figures/fig_tilt_law.pdf`).
- **Scope for improvement:** the R<sub>AA</sub> bump at 1.5–2.5 GeV (coalescence), the flat Λc/D⁰ centrality
  dependence, Ds/D⁰ (fragmentation strangeness) and the bulk directed flow (`docs/CHANGES.pdf`, paper).

## 5. Repository layout
```
patches/            the change to the seven stock AMPT files (unified diff against the official release)
src_new/            new files: charm_mod.f (charm physics, switches), coal_new.f (He-Lin coalescence),
                    hqctl.inc (common block), eps09s_f.cpp (bindings to the official EPS09s code), Makefile
inputs/             input cards: auau200_lf_s6.ampt (reference), auau200_sm.ampt (stock, public AMPT)
configs/            switch sets (*.env) for each production arm
scripts/            setup_deps.sh, build.sh, run_job.sh, run_arm.sh, reduce_dmeson_v2.awk, check_prereqs.sh
analysis/           analyze_arm.py + dmlib.py (production analysis), indep_check.py (independent re-analysis),
                    audit_claims.py (all chi2 tables), paper_figs.py, fig_tilt_law.py, compare_bulk.py, ...
data/hepdata/       STAR/PHENIX HEPData records used (CSV)
data/models/        published model curves digitised from the papers' vector graphics; each JSON records its
                    source and its calibration against the experiment's own markers (extract_*.py)
observation/results json/ (analysis outputs of the paper production), tables/, figures/, bulk/
tests/              the three test suites and the farm reference files
docs/               CHANGES.pdf (+ .tex and the generator of its code tables)
```

## 6. Provenance of external inputs (downloaded by `make deps`, never modified)

| input | source | SHA-256 |
|---|---|---|
| AMPT v1.26t9b/v2.26t9b | myweb.ecu.edu/linz/ampt | `9441b5f7…5293ed55` |
| LHAPDF 6.5.1 | lhapdf.hepforge.org | `1256419e…a6a66751` |
| CTEQ6.1M (`cteq61`) | lhapdfsets.web.cern.ch | `5a6f1393…a620b92` |
| EPS09s | users.jyu.fi/~kaeskola/EPS09s | `294bce5e…ffdf2c` |

## 7. Citing
If you use this code, please cite the AMPT papers (Lin et al., PRC 72 064901 (2005); Lin & Zheng, Nucl. Sci. Tech.
32, 113 (2021)), the papers listed in `docs/CHANGES.pdf` for the components you use, and our paper (reference to be
added).
