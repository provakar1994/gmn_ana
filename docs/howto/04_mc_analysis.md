# Stage 4 — MC analysis (SIMC / g4sbs)

Back to the overview: [../ANALYSIS_GUIDE.md](../ANALYSIS_GUIDE.md) ·
Previous: [./03_data_analysis.md](./03_data_analysis.md) · Next: [./05_checks.md](./05_checks.md)

## Purpose of this stage

This stage analyzes the replayed Monte Carlo files with the same reconstruction and cuts as the data
([./03_data_analysis.md](./03_data_analysis.md)). The result is a `Tout` tree per setting and sample that
can be compared with data event by event. Each MC event also gets a weight, so that the sum of weights
is an absolute yield. The fits of stage 6 use the MC dx distributions of protons and neutrons as
templates. The MC ratio of their normalizations gives the model n/p cross-section ratio R_MC, the
reference from which GMn is extracted. MC **generation, digitization and replay are not done in this
repository**. They are done in the companion MC production / consistency-audit repository ("the MC repo").
This repository starts from the replayed MC files and the job-summary CSVs that come with them.

## MC samples

| Sample (`generator_process`) | Target | Physics | Used by |
|---|---|---|---|
| `simc_heep` | LH2 | SIMC elastic H(e,e'p), with radiative effects | `elas_ana_simu.cpp` (W² offset, LH2 data/MC, proton efficiency) |
| `simc_neen` | — | SIMC elastic scattering on a free neutron, run like the LH2 sample but with radiation off (no radiative tail) | `elas_ana_simu.cpp` (HCAL NDE systematics) |
| `simc_deeN` = `deep` + `deen` | LD2 | SIMC quasi-elastic d(e,e'p) and d(e,e'n), two separate sub-samples | `qelas_ana_simu.cpp` (signal templates, R_MC) |
| `simc_deep`, `simc_deen` | LD2 | One of the two QE sub-samples alone | `qelas_ana_simu.cpp` |
| `g4sbs_inel` | LD2 | g4sbs inelastic generator (pion production) | `qelas_ana_simu.cpp` (inelastic background shape) |

SBS-11 data were taken partly with and partly without the Al shield beside the scattering chamber. The
SBS-11 MC is therefore a **charge-weighted mix** of an Al-shield sample and a no-shield sample. For LH2 the
weights are **0.607 (Al shield) / 0.393 (no shield)**. Each sample is analyzed separately. The mix is
made downstream (`fit_W2`, `w2_offset_check` and the dx fits via their config weights). In the SBS-11
MC configs, the default blocks (`simc_heep_model1`, `simc_deeN_model2`) are the Al-shield samples. The
no-shield samples are the `..._noshield` blocks (variant `noshield`; their prefixes start with `noshield_`).

## Typical workflow

1. In the MC repo: generate, digitize and replay the samples of each setting. Copy (or point
   `rootfile_dir` to) the replayed files **and** the job-summary CSVs.
2. For each setting, in the order SBS-4 (0p, 30p, 50p), SBS-9 (0p, 70p), SBS-8 (0p, 50p, 70p, 100p),
   SBS-14 (0p, 70p), SBS-7 (0p, 85p), SBS-11 (0p, 100p), run `elas_ana_simu.cpp` on `simc heep`
   (model 1). Use the same cuts, `bbfidu_cutR` and HCAL `dx_p_cut`/`dy_p_cut` as the data (the MC peak
   positions and widths differ, so fit them in the MC output and update the MC config).
3. Check / determine `W2_offset` with `w2_offset_check.cpp` on the LH2 data and MC outputs
   ([./05_checks.md](./05_checks.md)). Put the value in **all** MC blocks of that setting (LH2 and LD2)
   and rerun step 2.
4. For the production settings (SBS-4 30p/50p, SBS-9 70p, SBS-8 70p, SBS-14 70p, SBS-7 85p,
   SBS-11 100p) run `qelas_ana_simu.cpp` on `simc deeN` and on `g4sbs inel` (model 2).
5. For SBS-11, steps 2 and 4 with the default blocks give the Al-shield sample. Repeat them with
   variant `noshield` (e.g. `./submit-run-ana-simu.sh qelas 11 100 simc deeN noshield`) for the
   no-shield sample. The two outputs are then mixed with the charge weights in the downstream configs.
6. On the ifarm, steps 2 and 4 are normally submitted with `submit-run-ana-simu.sh`.
7. Pass the `siout/*.root` files to the fits ([./06_fits.md](./06_fits.md)).

Configs are in `scripts/simulation/configs/sbs<conf>-sbs<field>p-simu/`. Existing MC configs, in Q² order:

| Setting | `conf_elas_ana_simu.json` | `conf_qelas_ana_simu.json` |
|---|---|---|
| SBS-4 | 0p, 30p, 50p | 0p, 30p, 50p |
| SBS-9 | 0p, 70p | 70p |
| SBS-8 | 0p, 50p, 70p, 100p | 0p, 50p, 70p, 100p |
| SBS-14 | 0p, 70p | 70p |
| SBS-7 | 0p, 85p | 85p |
| SBS-11 | 0p, 100p | 100p |

### Config blocks and keys

Block name = `<generator>_<process>_model<model>[_<variant>]` (e.g. `simc_heep_model1`,
`simc_deeN_model2`, `g4sbs_inel_model2`, `simc_deeN_model2_noshield`). The script builds this name
from its arguments; the optional last argument `variant` (default `""`) adds the suffix. The block
without a suffix is the **production** sample. Other samples are suffixed blocks, e.g.
`simc_heep_model1_noshield` / `simc_deeN_model2_noshield` (SBS-11 no shield),
`simc_heep_model1_sbs_OFF` (SBS-9 and SBS-14 70p field-off LH2 samples, `SBS_scale_field` = 0),
`..._thesis`, `..._radc`. Still check the `outfile_prefix`/`prefix_to_filebase` of a block before a
production run.

| Key | Meaning |
|---|---|
| `Ntype` | `p` for heep/neen, `np` for LD2. |
| `SBS_config`, `SBS_magnet_percent` | Setting. |
| `SBS_scale_field` | SBS field scale factor used in the g4sbs simulation (`/g4sbs/scalesbsfield`, e.g. 0.97 ↔ the `0p97sf` tag in file names). It is used to compute the expected proton deflection `p_def`. |
| `rootfile_dir` | Directory with the replayed MC files. For SIMC, the job summaries are in its `simcout/` subdirectory. |
| `prefix_to_filebase` | File prefix(es). For `deeN`: `[deep prefix, deen prefix]`. |
| `Njobs_to_ana` | Number of jobs to use (−1 = all). |
| `global_cut` | Raw-tree cut. **Must start with `bb.tr.n>0`** (see notes). |
| `W2_offset` | Constant added to the reconstructed W² (GeV²). `W2` in `Tout` = `W2_recon + W2_offset`. `W` is **not** shifted. |
| `dy_offset` | Constant added to the stored `dy` (not to the dy used for `pCut`). |
| `W2_cutR`, `bbfidu_cutR`, `dx_p_cut`, `dy_p_cut`, `dx_n_cut`, `dy_n_cut` | As in the data configs, but with the MC peak positions and widths. |
| `hcal_voffset`, `hcal_hoffset`, `hcal_zoffset` | HCAL position offsets for the MC geometry. |
| `h_*_lims`, `h_*_fitR` | Diagnostic histograms and fits only. |
| `outfile_prefix`, `verbose`, `verbose_function`, `max_N_tracks` | As for data. |

### How the replayed MC files are found (`util_pd::ReadSimuJobSummary`)

With `filebase = <prefix>_sbs<conf>_sbs<field>p_<generator>`:

- SIMC job summary: `<rootfile_dir>/simcout/<filebase>_<process>_summary.csv` (for `deeN`, two files:
  `..._deep_summary.csv` with the first prefix and `..._deen_summary.csv` with the second).
- g4sbs job summary: `<rootfile_dir>/<filebase>_<process>_summary.csv`.
- Each summary row is one job. The replayed file is `<rootfile_dir>/replayed_<filebase>_<process>_job_<jobid>.root`.
- Summary columns (SIMC): `jobid, Ngen_requested, Ntried, genvol, lumi, Ebeam, charge(mC), ..., using_RS, max_wt_RS`.
  For g4sbs, column 7 is the beam current (µA) instead of the charge.
- Jobs whose replayed file is missing, unreadable or has an empty tree are skipped (with a message).
  Jobs with far fewer entries than requested are kept but flagged. A missing summary file is a hard error.

### MC normalization (per-event weight)

For each job, the script takes `genvol` (ΔΩ: MeV·sr² for SIMC), `lumi` (integrated luminosity, µb⁻¹ for
SIMC), `charge` and `Ebeam` from the summary. The total number of tries N_tried is summed over all
**usable** jobs (separately for `deep` and `deen` in a `deeN` run). Then

```
weight      = σ_event · genvol · lumi / N_tried(total)      (σ_event = MC.simc_Weight for SIMC, MC.mc_sigma for g4sbs)
weight_norm = weight / charge_job                          (yield per coulomb; compare with data weight_norm = 1/(charge·livetime))
```

If the SIMC sample was generated with rejection sampling (`using_RS` = 1), σ_event is replaced by the
constant `max_wt_RS`. Weight histograms with `weight` (shapes, absolute MC yield) or `weight_norm`.

`weight_effic` / `weight_norm_effic` = the same weights multiplied by the data proton detection efficiency
at the true expected HCAL position, relative to its average (named constant `kHCALAvgEffi` = 0.9517, the mean of the map's constant fits vs x and vs y; update it together with the map). The efficiency comes from the map
`h2_effi_map` of `pDE_map.cpp` (SBS-8 LH2, all fields, pass 3). Its file is set by the named constant
`kHCALEffiMapFile` at the top of both scripts:
`../pDEmap/0p77zoff_pDE_data_sbs8_sbs-1p_model1_pass3.root`, relative to `scripts/simulation/`,
i.e. `scripts/pDEmap/`, where the map is tracked in git (how it was made: see
[./05_checks.md](./05_checks.md)). The macro stops if the file is missing.

### MC-specific `Tout` branches (in addition to the data branches of stage 3)

`weight`, `weight_norm`, `weight_effic`, `weight_norm_effic`. `mc_fnucl` (struck nucleon: 1 = p,
0 = n). `pN_t`, `thN_t`, `phN_t` (true nucleon momentum and angles). `xHCAL_exp_t`, `yHCAL_exp_t` (true
expected HCAL position). SIMC only: `veE`, `vetheta`, `vQ2`, `vepsilon` (vertex values). `phN_exp`.
There are no run/scaler/coincidence-time branches and no offline HCAL cluster choice: MC uses the
default (highest-energy) cluster.

---

### `scripts/simulation/elas_ana_simu.cpp`
**What it does:** MC counterpart of `elas_ana_data.cpp` for LH2 (heep, also neen). It reconstructs the
same kinematics and HCAL quantities, adds `W2_offset` to W², computes the per-event normalization
weights and the efficiency-weighted versions, and stores MC truth. It also makes weighted diagnostic
histograms and the dx/dy/W² fits.
**Run:** from `scripts/simulation/` (library loaded via the rootlogon, see [./01_setup.md](./01_setup.md)):
```
analyzer -b -q <rootlogon.C> 'elas_ana_simu.cpp("configs/sbs8-sbs70p-simu/conf_elas_ana_simu.json",1,"simc","heep")'
```
Arguments: `configfilename`, `model` (default 1), `generator` (`simc`/`g4sbs`), `process` (`heep`, `neen`, ...;
the default in the code is `neen`, so always pass it), optional `variant` (default `""`; e.g.
`"noshield"`, `"sbs_OFF"`).
**Inputs:** config block `<generator>_<process>_model<model>[_<variant>]`. Job-summary CSV(s) and
replayed files in `rootfile_dir`. The data efficiency map `kHCALEffiMapFile` (see above).
**Outputs:** `siout/<outfile_prefix>_elas_ana_<generator>_sbs<conf>_sbs<field>p_model<model>.root`
(+ `.pdf`). `inel_` is added to the prefix for `process = inel`.
**Feeds:** `w2_offset_check.cpp`, `pDE_simu.cpp`, `pDE_map_simu.cpp`, `eloss_check_simu.cpp` context
([./05_checks.md](./05_checks.md)). LH2 data/MC comparisons and `fit_W2.cpp` ([./06_fits.md](./06_fits.md)).
**Notes:**
- **`global_cut` must start with `bb.tr.n>0`.** Unlike the data replay, the MC replay keeps events
  without a track. The scripts evaluate the cut with `TTreeFormula::EvalInstance(0)`, which for a
  trackless event reads the previous event's track values. A trackless event could then pass and pair
  its HCAL cluster with another event's electron.
- No target energy-loss correction is applied in MC. For SIMC the beam energy is the per-event
  `MC.simc_Ebeam`. For g4sbs it is the summary `Ebeam`. Use `eloss_check_simu.cpp` to see how the MC
  momentum reconstruction treats the outgoing loss.
- `W2_offset` is the same for LH2 and LD2 of a setting: it is measured on LH2 by `w2_offset_check` and
  copied to the LD2 blocks. `W2` carries the offset and `W` does not, so check W² with W and W² cuts in mind.
- `xHCAL_exp_p` = `xHCAL_exp − p_def` of the same event (filled after `p_def` is computed).

### `scripts/simulation/qelas_ana_simu.cpp`
**What it does:** MC counterpart of `qelas_ana_data.cpp` for LD2. It runs on SIMC `deeN` (proton and
neutron QE sub-samples, each normalized with its own total N_tried) or on g4sbs `inel`. It writes the
same `Tout` structure, plus `nCut`, the neutron spot quantities and model cross-section branches
(`sigMott`, `sigRed_p`, `sigRed_n`, `sigBorn_ratio`). `mc_fnucl` lets you split the dx distribution into
proton and neutron templates.
**Run:** from `scripts/simulation/`:
```
analyzer -b -q <rootlogon.C> 'qelas_ana_simu.cpp("configs/sbs8-sbs70p-simu/conf_qelas_ana_simu.json",2,"simc","deeN")'
analyzer -b -q <rootlogon.C> 'qelas_ana_simu.cpp("configs/sbs8-sbs70p-simu/conf_qelas_ana_simu.json",2,"g4sbs","inel")'
```
For the SBS-11 no-shield sample add the variant: `...,2,"simc","deeN","noshield")`.
Arguments: `configfilename`, `model` (default 2), `generator`, `process` (`deeN`, `deep`, `deen`, `inel`),
optional `variant` (default `""`).
**Inputs:** config block `<generator>_<process>_model<model>[_<variant>]` (adds `dx_n_cut`, `dy_n_cut`,
`h_dxHCAL_n_fitR`). Summaries + replayed files. The efficiency map `kHCALEffiMapFile`.
**Outputs:** `siout/<outfile_prefix>_qelas_ana_<generator>_sbs<conf>_sbs<field>p_model<model>.root` (+ `.pdf`).
`inel_` is added to the prefix for g4sbs inelastic.
**Feeds:** dx fits (p/n templates, inelastic background) and R_MC → GMn ([./06_fits.md](./06_fits.md), [./07_results.md](./07_results.md)).
**Notes:** The same `bb.tr.n>0` rule and the same `W2_offset` rule apply as for `elas_ana_simu.cpp`.
For `deeN`, the deep and deen jobs must both have summaries in `simcout/`. `xHCAL_exp_p`
(= `xHCAL_exp − p_def` of the same event) is filled here too.

### `scripts/simulation/run-ana-simu.sh`
**What it does:** Batch wrapper for one MC job on the ifarm. It loads `analyzer`, goes to
`$GMN_ANA/scripts/simulation` and runs `elas_ana_simu.cpp` (model 1) or `qelas_ana_simu.cpp` (model 2).
**Run:** `sbatch run-ana-simu.sh <elas|qelas> <conf> <field> <simc|g4sbs> <heep|deeN|inel|...> [<variant>]`,
e.g. `sbatch run-ana-simu.sh qelas 4 30 simc deeN` or `sbatch run-ana-simu.sh qelas 11 100 simc deeN noshield`.
**Inputs:** the setting's MC config `configs/sbs<conf>-sbs<field>p-simu/conf_<process>_ana_simu.json`.
The environment variables `GMN_ANA` and `ROOTLOGON_ENV`: `export GMN_ANA=... ROOTLOGON_ENV=...` before
submitting (`sbatch` passes the environment). If unset, the original author's ifarm paths are used.
**Outputs:** as for the macro (`siout/...`).
**Feeds:** —
**Notes:** The script exits with an error if `$GMN_ANA/scripts` or `$ROOTLOGON_ENV` is not found.
The model is fixed (1 for elas, 2 for qelas). An empty `variant` selects the default (production) block.

### `scripts/simulation/submit-run-ana-simu.sh`
**What it does:** Submits `run-ana-simu.sh` with job name `simu_<process>_sbs<conf>_f<field>p_<gen>_<reaction>[_<variant>]`, 2 GB/CPU.
**Run:** `./submit-run-ana-simu.sh elas 4 30 simc heep`, `./submit-run-ana-simu.sh qelas 4 30 simc deeN`,
`./submit-run-ana-simu.sh qelas 4 30 g4sbs inel`, `./submit-run-ana-simu.sh qelas 11 100 simc deeN noshield`
(optional 6th argument = variant).
**Inputs:** as `run-ana-simu.sh` (set `GMN_ANA` and `ROOTLOGON_ENV` first).
**Outputs:** Slurm logs in `/farm_out/$USER/gmn_ana/slurm_out/` (override with `GMN_SLURM_LOGDIR`; the script creates it).
**Feeds:** —

### `scripts/simulation/g4mac/` (`get_SBSconf_info.cpp`, `sbsfield_scan_genmac.py`, example macros)
**What it does:** Helpers for setting up g4sbs simulations (the production itself is done in the MC repo).
`get_SBSconf_info.cpp` prints the nominal geometry of a setting from the library's `SBSconfig`
(beam energy, BigBite/SBS angles and distances, HCAL distance), so that you can copy the numbers into a
g4sbs macro. `sbsfield_scan_genmac.py` copies a g4sbs macro once per SBS field scale factor (replacing
the `/g4sbs/scalesbsfield` value) to scan the field. The `.mac` / `.inp` files are old examples (SIMC-input
g4sbs run for SBS-4, QE runs, nucleon guns for HCAL efficiency, SIMC input deck).
**Run:** `root -l -b -q 'get_SBSconf_info.cpp(8,70)'` (args: config, field %) from `scripts/simulation/g4mac/`.
`python3 sbsfield_scan_genmac.py`: the input macro, output base, old value and list of new values are
set in `main()` (edit before running).
**Inputs:** none / a g4sbs macro in the current directory.
**Outputs:** printout / `<scale>sf_<outfilebase>.mac` files.
**Feeds:** MC production in the MC repo.
**Notes:** The example macros contain hard-coded paths of the original account and are not kept up to date.

### `scripts/simulation/hcal/gen_xypos_DB.cpp`
**What it does:** Prints the HCAL block-centre positions (`sbs.hcal.xpos` / `sbs.hcal.ypos`, 24 rows × 12
columns) for a given vertical/horizontal HCAL offset, in the replay's coordinate convention (+x down,
+y towards the beamline). This makes the HCAL position database entries for the MC replay.
**Run:** `root -l -b -q 'gen_xypos_DB.cpp(0.75,0.0)'` from `scripts/simulation/hcal/` (args: vertical, horizontal offset in m).
**Inputs:** none.
**Outputs:** printout to paste into the replay DB (in the MC repo).
**Feeds:** MC replay → `hcal_voffset`/`hcal_hoffset` used here.

### `scripts/pcalib/` and `scripts/simulation/pcalib/` (BigBite momentum / optics calibration)
**What it does:** BigBite momentum calibration using elastic kinematics (from A. Puckett's method). The
first-order momentum coefficients are fitted so that the reconstructed p matches the elastic
momentum at the measured angle, and a new optics/momentum coefficient file is written for the replay.
- `scripts/pcalib/MakeMomentumCalibrationTree.C` — real data (LH2). Setup example: `setup_Pcalib_SBS4.txt`
  (file list, cut, then keys such as `runsheet_dir`, `rootfile_dir`, `target`, `ebeam`, angles, `dEdx`,
  `useAlShield`, `oldcoeffs`, W and `dpel` windows).
- `scripts/simulation/pcalib/momentum_calib_simu.cpp` — the same for replayed MC (pseudo-data). Setups:
  `setup_Pcalib_SBS4.txt`, `setup_Pcalib_SBS14.txt`.
- `scripts/pcalib/Optics_GEN.C` — fits the MC BigBite optics (angles, vertex, momentum) on g4sbs electron-gun
  events against the true vertex values. Setups: `setup_OnPcalib_sbs14_sbs0p_g4sbs_egun.txt`,
  `simulation/pcalib/setup_OpticsCalib_MC_SBS14.txt`. g4sbs macro: `g4mac_OnPcalib_sbs14_sbs0p_g4sbs_egun.mac`.
  It has arguments `(inputfilename, outputfilename, NMAX=1000000, BBscale=0.97)`.
**Run:** e.g. `root -l -b -q 'momentum_calib_simu.cpp("setup_Pcalib_SBS14.txt","pcsout/NewMomentumFit.root")'`
from `scripts/simulation/pcalib/`. Data: `MakeMomentumCalibrationTree.C("setup_Pcalib_SBS4.txt","pcdout/NewMomentumFit.root")`
from `scripts/pcalib/` (verify the include path for `CodaRun.h`/`Utilities.h`).
**Inputs:** setup text file (file list until `endlist`, cut until `endcut`, then key/value settings).
`oldcoeffs` = existing coefficient file whose angle/vertex terms are kept.
**Outputs:** `<output>.root` (trees, histograms), `<output>.dat` (new coefficients, `bb.preconflag = 1` format),
`<output>.pdf`. `Optics_GEN.C` writes `<output>.txt` and `foptics_<output>.txt`. The output directories
(`pcdout/`, `pcsout/`) must exist.
**Feeds:** replay databases (data replay / MC replay in the MC repo). The MC optics fit is the
reference for the energy-loss question studied by `eloss_check_simu.cpp`.
**Notes:** `Optics_GEN.C`, `setup_OnPcalib_*.txt` and `g4mac_OnPcalib_*.mac` are **not tracked in git**
at the time of writing. `Optics_GEN.C` includes `gen_tree.C`, which is not in this repository (verify).
The setup files contain old absolute paths.
