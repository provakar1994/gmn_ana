# Stage 7: Results (R_n/p table → GMn → plots)

[Overview](../ANALYSIS_GUIDE.md) · previous stage: [06 Fits](./06_fits.md)

## Purpose of this stage

This stage turns the `fit_dx.cpp` outputs into physics results. For each production
setting, the fitted data/MC scale factor R_fit and the MC Born ratio R_MC give the measured
n/p Born cross-section ratio **R = σ_n/σ_p = R_MC · R_fit**. R and its statistical and
systematic errors are collected in one table, `scripts/lta/results.csv`. GMn/(μ_n G_D) is
then extracted from it for several proton/GEn form-factor models, with full error
propagation (`scripts/gmn_extraction/`). Finally the GMn points can be used for a u/d flavor
decomposition (`scripts/flavordecomp/`). Only the LD2 data of the production settings enter
GMn: SBS-4 30p, SBS-4 50p, SBS-9 70p, SBS-8 70p, SBS-14 70p, SBS-7 85p, SBS-11 100p. SBS-8
is formally an nTPE point and is carried as a check.

## Typical workflow

1. **Collect the fit outputs.** For each production setting (in Q² order, as above) gather
   the `fit_dx.cpp` CSVs from `pdout/fits/sbs<conf>sbs<mag>p/` ([06](./06_fits.md)): the
   production run plus the variants "OTm3" (RC prescription), `_noefficorr` and
   `_nonblkHCALcut`. `scripts/lta/grab_rnpsf.ipynb` puts them side by side.
2. **Build `scripts/lta/results.csv`** (one row per setting). The central value is
   R = R_fit × R_MC, the stat error is δR_fit × R_MC, and the systematics are absolute
   errors on R (conventions below). This step is currently done in a local notebook
   (`get_results.ipynb`) and partly by hand; see the `results.csv` card. The notebook and
   the result CSVs are intentionally not tracked in git yet: they will change, and the
   table should eventually be produced directly by this repository.
3. **Extract GMn** for all points and models:
   `extractGMn_csv.cpp("../lta/results.csv", "gmn_results.csv")` (see the
   [gmn_extraction README](../../scripts/gmn_extraction/README.md)). Use `extractGMn.cpp`
   for a single point or a quick what-if.
4. **Plots.** Use the GMn table for the final figures. `get_ff.cpp` writes model curves to
   CSV for plotting, and `customize_plots.cpp` restyles saved ROOT canvases. Run
   `flavordecomp.cpp` for the flavor-decomposition figures (the GMn values are typed into
   the macro, see its card).

### The chain in one picture

```
fit_dx.cpp CSV:  Rk ± Rkerr (k = bg model), RpMC  (per setting, per variant)
        │   R = R3·RpMC (pol2 bg = central),  R_stat = R3err·RpMC
        │   systematics: absolute δR per source (bg-model spread, NDE, RC, ...)
        ▼
scripts/lta/results.csv   (conf, field, etheta, q2, eps, R, R_stat, <sys...>)
        ▼
extractGMn_csv.cpp  →  gmn_results.csv  (GMn/(μn GD) per point × model, all errors)
        ▼
plots, flavordecomp.cpp
```

### Systematic-error conventions

These are as implemented in `extractGMn_csv.cpp` and documented in the
[gmn_extraction README](../../scripts/gmn_extraction/README.md):
- All systematic columns of `results.csv` are **absolute errors on R** (the same units as
  `R_stat`). They are **not** relative and **not** on R_fit.
- A column whose name ends in **`_sc`** is a **scale** (normalization, correlated across
  points) error. All others are **point-to-point (p2p)**. A source can have both parts
  (e.g. `NDE2` + `NDE2_sc`, `RC` + `RC_sc`). In the GMn output the p2p part is then renamed
  `X_p2p`.
- All sources are treated as uncorrelated with each other. Totals are sums in quadrature:
  `TotSys_p2p` (p2p sources + the FF errors pCS and GEn), `TotSys_sc` (scale sources),
  `StatP2P` = stat ⊕ TotSys_p2p and `TotErr` = everything. The `_sc` errors are never
  folded into the p2p total.
- An empty cell is read as 0 with a warning. `Cut V` is currently empty (not yet computed).
- In the GMn output every error is on GMn/(μ_n G_D).

---

### `scripts/lta/results.csv` (untracked, central table)

**What it does:** Holds the measured n/p Born ratio R for every production setting, with
its stat and systematic errors. It is the only input to the GMn extraction. The file is
intentionally **not tracked in git** yet: it will still change and should eventually be
produced directly by this repository. Keep a copy with every result you quote.

**Run:** not a script. It is produced as described in workflow step 2.

**Inputs:** the `fit_dx.cpp` CSVs (production + variants), and the external NDE/RC
percentages (local file `data_external_dRsys_in_pct_MCp4_Dp3.csv`, see below). `NDE1`,
`NblkC` and `FSI` currently come from an external source; the author plans to bring them
into this repository later.

**Outputs:** the table itself. **Columns** (rows in Q² order SBS-4 30p, 4 50p, 9 70p, 8 70p, 14 70p, 7 85p, 11 100p):

| Column | Meaning | Units / origin |
|---|---|---|
| `conf`, `field` | SBS config and SBS magnet % | labels |
| `etheta` | electron angle | rad; information only, recomputed from (q2, eps) in the extraction |
| `q2`, `eps` | acceptance-averaged Q² and ε_p, which define the kinematic point | GeV², – |
| `R`, `R_stat` | σ_n/σ_p (Born) and its stat error | = R3·RpMC and R3err·RpMC (fit 3, pol2 bg = production bg shape) |
| `Inel` | background-model systematic | p2p; spread of R1–R7 × RpMC (fits 1–7) |
| `NDE1` | HCAL detection-efficiency systematic | p2p; currently from an external source (to be brought into this repo later) |
| `NDE2`, `NDE2_sc` | HCAL NDE, external input: 1.41% p2p + 0.71% scale of R | p2p / scale |
| `RC`, `RC_sc` | radiative corrections, external: 0.75% p2p + 1.75% scale of R | p2p / scale |
| `NblkC` | `nblkHCAL>1` cut systematic | p2p; currently from an external source (to be brought into this repo later) |
| `Cut V` | cut-variation systematic | p2p; **empty, not yet filled** |
| `FSI` | final-state-interaction systematic | p2p; currently from an external source (to be brought into this repo later) |

**Feeds:** `extractGMn_csv.cpp`.

**Notes:**
- **R_MC fix pending.** The R values in `results.csv` were made with the pre-fix
  `fit_dx.cpp` R_MC. Until the fits are rerun, the extraction uses the local
  `results_RMCfix.csv`. It has the same columns plus `RpMC_old`, `RpMC_new` and
  `RMC_scale`, with R and all errors rescaled by R_MC(new)/R_MC(old). See "Open items" in
  the [gmn_extraction README](../../scripts/gmn_extraction/README.md).
- q2/eps are acceptance-averaged vertex quantities written by `fit_dx.cpp` (`vQ2Mean`,
  `vepsMean` in the fit CSV). In the current file they are the same for SBS-4 30p and 50p.
- Other local, untracked working files in `scripts/lta/` that this stage relies on:
  `get_results.ipynb` (computes R = Rk·RpMC for all bg models from one fit CSV),
  `bg_map.csv` (index → bg model name: R3 pol2, R5 pol3, R4 gaus, R1 anti-dy, R6 anti-dt,
  R2 inel-MC, R7 inel-MC+rndm) and `data_external_dRsys_in_pct_MCp4_Dp3.csv` (the external
  NDE2/RC percentages for p, n and ratio, p2p and scaled). Like `results.csv` (and
  `results_RMCfix.csv`), they are intentionally not tracked yet: they will change and should
  eventually be produced directly by this repository.

---

### `scripts/lta/grab_rnpsf.ipynb`

**What it does:** A Python/pandas notebook. It reads the `fit_dx.cpp` CSVs of all production
settings for four run modes, "Prod" (`OT0nm3`), "OTm3", "NOeffic" (`_noefficorr`) and
"NOnblkHCALcut" (`_nonblkHCALcut`). It builds a table of R_fit (`R3, R5, R4, R1, R6, R2, R7`)
and errors per setting × bg model × mode. This is the raw material for the central value
and the bg-model systematic, and a cross-check of the NDE and nblk variants (`NDE1` and
`NblkC` in `results.csv` currently come from an external source).

**Run:** open in Jupyter. Edit the `file_groups` dict (keys `430, 450, 970, 870, 1470, 785, 11100`
= setting + field, values = list of (CSV path, mode)) to point at your `pdout/fits/`
directory, then run all cells.

**Inputs:** the `fit_dx.cpp` CSVs.

**Outputs:** `df_final` (Prod value ± err and the three variant values). Writing
`out_grab_rnpsf.csv` is commented out.

**Feeds:** building `results.csv`.

**Notes:** the paths are absolute paths on the analysis farm and must be changed. The
values are R_fit, **not** R (multiply by RpMC).

---

### `scripts/lta/mcp4_OT0nm3_fit_dx_qelas_pass3_simc_sbs4_sbs50p_model2.csv`

**What it does:** A tracked example of a `fit_dx.cpp` QE output CSV (SBS-4 50p, production
mode, pass 3). It documents the column layout (see [06](./06_fits.md)) and is the reference
used to check the old R_MC reconstruction.

**Run / Inputs / Outputs:** none (data file). **Feeds:** reference only.

**Notes:** its `RpMC` (0.3957) is from the old R_MC code.

---

### `scripts/gmn_extraction/extractGMn_csv.cpp`

**What it does:** Extracts GMn/(μ_n G_D) for every row of a results CSV and every EMFF model
(8 models "A-B": proton FFs from A, GEn from B; primary = **Ye noTPE-Ye**). It propagates
stat, every p2p and scale systematic, and the proton-cross-section (pCS) and GEn
form-factor errors. Full physics,
models, caveats and formulas: **[gmn_extraction README](../../scripts/gmn_extraction/README.md)**
(not repeated here).

**Run** (from `scripts/gmn_extraction/`, gmn_ana installed):
```
root -l -b -q <install>/runenv/rootlogon.C \
  'extractGMn_csv.cpp("../lta/results.csv", "gmn_results.csv")'
```
arg 1: input results CSV (use `../lta/results_RMCfix.csv` until the fits are redone);
arg 2: output CSV.

**Inputs:** the results CSV (columns are found by name; `q2, eps, R, R_stat` are required,
every column after `R_stat` is a systematic, and the columns before it are labels). It also
reads the FF fit tables in `DB/FF_fits/` (path compiled in, override with env
`GMNANA_DB_DIR`).

**Outputs:** the output CSV in the current directory, one row per point × model:
`<labels>,Q2,eps,etheta_eff,R,R_stat,Model,TPE,Val,StatErr,<p2p>,pCS,GEn,<scale>,TotSys_p2p,TotSys_sc,StatP2P,TotErr`.

**Feeds:** the final GMn plots and tables, and the values typed into `flavordecomp.cpp`.

**Notes:** Ye "noTPE" is still a placeholder (it is effectively TPE corrected), and the
library prints a warning. pCS is likely underestimated. Both are described in the README.

---

### `scripts/gmn_extraction/extractGMn.cpp`

**What it does:** The same physics for **one** kinematic point, starting from R_fit (not R).
It prints R_MC, R, GMn for all models, and the error budget of the primary model.

**Run** (from `scripts/gmn_extraction/`):
```
root -l -b -q <install>/runenv/rootlogon.C \
  'extractGMn.cpp(4.5, 0.6833, 1.02, 0.01, {0.004,0.003,0.002}, {"Inel","NDE1","RC_sc"})'
```
The arguments are `Q2` (GeV²), `eps`, `R_fit`, `R_fit_stat`, a vector of systematic errors
**on R_fit** (absolute), and optionally their names (`_sc` = scale).

**Inputs:** command-line numbers only (plus the FF tables in `DB/`).

**Outputs:** terminal printout, including a CSV block.

**Feeds:** quick checks; not used in the batch chain.

---

### `scripts/lta/get_ff.cpp`

**What it does:** Tabulates one nucleon form factor vs Q² (1 to 15 GeV², 0.1 steps) for
several parametrizations, for plotting reference curves. GEp/GMp: Kelly 2004, Ye 2017/18,
Christy 2022 (no TPE). GEn: Kelly, Ye, Riordan/"Seamus", Galster. GMn: Kelly, Ye.

**Run** (from `scripts/lta/`):
```
root -l -b -q <install>/runenv/rootlogon.C 'get_ff.cpp(1.0, "GMn")'
```
arg 1 `Q2` is unused by the current code (kept for the commented-out printouts); arg 2 is one
of `GEp`, `GMp`, `GEn`, `GMn`.

**Inputs:** FF fits from the library / `DB/FF_fits/`.

**Outputs:** `<G>_parametrization.csv` in the current directory (absolute FF values, not
divided by the dipole).

**Feeds:** model curves on the GMn and FF plots.

---

### `scripts/lta/customize_plots.cpp`

**What it does:** Restyles a canvas saved in a ROOT file for talks and papers. It sets
fonts, axis title sizes, margins, grid and ticks and removes stat boxes, on every pad.
The input file and canvas name (`c0`) are hard-coded in `customize_plots()`.

**Run** (from `scripts/lta/`): edit `inFile` / the canvas name, then
`root -l <install>/runenv/rootlogon.C customize_plots.cpp`.

**Inputs:** a ROOT file with saved canvases, typically a `fit_dx.cpp` output under
`../pdout/fits/...`.

**Outputs:** none written; the canvas is drawn on screen, and you save it yourself.

**Feeds:** figures.

**Notes:** the hard-coded default file is an old pass-2 SBS-4 0p fit.

---

### Other files in `scripts/lta/`

`fEvtHdr.cpp` (prints run/event number from a replayed tree) and `rdf.cpp` (an RDataFrame
scratch example) are small utilities and are not part of the results chain.

---

### `scripts/flavordecomp/flavordecomp.cpp`

**What it does:** Flavor decomposition of the nucleon form factors, assuming charge
symmetry and no strange quarks. It builds the Dirac/Pauli FFs
F1 = (GE + τGM)/(1+τ) and F2 = (GM − GE)/(1+τ) for p and n, then
F^u = 2F^p + F^n, F^d = 2F^n + F^p, and the Sachs combinations GM^u, GE^u, GM^d, GE^d. It
uses Ye 2018 for GEp, GMp and GEn. For GMn it uses (a) Ye, a smooth curve; (b) **our SBS
GMn points**; and (c) the Hague GMn table. It compares with Cates et al. 2011. The plots
show Q⁴F1^u, Q⁴F2^u/κ_u, 2.5·Q⁴F1^d and 0.75·Q⁴F2^d/κ_d (κ_u = 1.67, κ_d = −2.03), plus
GM^q/(μ G_D) and μ·GE^q/GM^q.

**Run** (from `scripts/flavordecomp/`):
```
root -l -b -q <install>/runenv/rootlogon.C flavordecomp.cpp
```
no arguments.

**Inputs:**
- SBS GMn points **hard-coded** in the macro (`Q2v`, `SBSGMnovMuGD`, `SBSGMnErrovMuGD`;
  currently the "pass3" set at Q² = 2.989, 4.488, 7.464, 9.834, 13.465 GeV²). They are an
  old set; the macro will be revisited, and should then use the final GMn results
  (`gmn_results.csv`).
- `GMn_lookup.csv` (Hague GMn/(μ_n G_D) table) and `f1f2_cates.csv` (Cates 2011 data), read
  from the current directory. **Neither is in git**; supply them locally.

**Outputs:** `c1_pass3.pdf` … `c10_pass3.pdf` in the current directory.

**Feeds:** final flavor-decomposition figures.

**Notes:** the proton FFs come from the TPE-corrected Ye fit (`GetFFwErr`), not from the
"noTPE" set used in the GMn extraction. (verify). Errors are propagated in quadrature without
correlations. There are five SBS points, one per Q²: SBS-4 (one value for 30p/50p), SBS-9, SBS-14,
SBS-7 and SBS-11. SBS-8 is not included.
