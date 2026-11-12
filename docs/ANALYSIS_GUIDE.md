# GMn analysis guide

Entry point for the `gmn_ana` repository (SBS GMn / nTPE analysis, Hall A, experiment E12-09-019).
Read this page first, then the stage file of the step you want to run.

## What is measured

The neutron magnetic form factor GMn is extracted from the ratio of quasi-elastic (QE)
e-n to e-p yields on the deuterium (LD2) target. The electron is detected in BigBite (BB);
the nucleon is detected in the hadron calorimeter HCAL behind the SBS dipole, which deflects
protons vertically, so protons and neutrons separate in HCAL `dx` (measured minus expected
vertical position). The data `dx` distribution is fit with the MC proton and neutron `dx`
shapes plus a background; the fitted n/p scale-factor ratio R_fit, times the MC Born
cross-section ratio R_MC, gives the measured n/p cross-section ratio R_n/p, from which GMn is
extracted using the proton form factors and GEn. LH2 elastic data (e-p only) are used to
calibrate and check the analysis: the MC `W2_offset`, HCAL proton/neutron detection
efficiency (NDE), HCAL position offsets, cut and fiducial-region definitions.

## Analysis chain

Steps marked `[outside]` are done in other repositories/software; only pointers are given.

```
 RAW DATA (EVIO, tape)                               MC GENERATION  [outside, "the MC repo"]
   |                                                   SIMC: QE d(e,e'p), d(e,e'n); elastic H(e,e'p)
   | [outside] detector calibration + replay           g4sbs: inelastic, particle guns, optics e-gun
   |   (SBS-offline / SBS-replay) -> pass-N             -> digitization -> replay (SBS-replay)
   |   ROOT files (data pass 3 is current)              -> MC ROOT files + job-summary CSVs
   |                                                     (MC pass 4 is current); data/MC
   | [outside, partly here] BB optics/momentum            consistency decisions are made in
   |   calibration (scripts/pcalib, simulation/pcalib)    the MC repo's consistency audit
   v                                                   v
 02 run info: good run lists, beam energy,          04 MC analysis: elas_ana_simu (LH2),
    charge, live time (DB/, scripts/epics)             qelas_ana_simu (LD2)  -> siout/
   v                                                   |
 03 data analysis: elas_ana_data (LH2),                |
    qelas_ana_data (LD2)  -> pdout/                    |
   |                                                   |
   +-------------------------+-------------------------+
                             v
 05 checks: MC W2_offset (w2_offset_check), e- energy loss (eloss_check_simu),
    HCAL proton detection efficiency (pDE_*)
                             v
 06 fits: fit_W2 (W2), fit_dx (HCAL dx -> R_fit and systematics per setting), fit_dy
                             v
 07 results: lta/results.csv (R_n/p + errors per setting)
             -> gmn_extraction/ (GMn for several FF models, full error propagation)
             -> flavordecomp/ (u/d flavor separation)
```

## Kinematic settings

Always ordered by Q², field settings ascending within a setting. "SBS field" = SBS magnet
current in % of 2100 A. E_beam and angles are the nominal values in `src/ExpConstants.cpp`
(E_beam = LD2 run-list average; the analysis itself uses the per-run value from the run
list). Q² is the acceptance-averaged value from `scripts/lta/results.csv`.

| Setting | E_beam (GeV) | BB angle (deg) | SBS angle (deg) | Q² (GeV²) | SBS fields with data | Production (GMn) |
|---|---|---|---|---|---|---|
| SBS-4  | 3.7393 | 36.0 | 31.9 | 2.99  | 0p, 30p, 50p       | **30p, 50p** |
| SBS-9  | 4.0268 | 49.0 | 22.5 | 4.49  | 70p                | **70p** |
| SBS-8  | 5.9826 | 26.5 | 29.9 | 4.53  | 0p, 50p, 70p, 100p | **70p** (check, see below) |
| SBS-14 | 5.9828 | 46.5 | 17.3 | 7.46  | 0p (LH2 only), 70p | **70p** |
| SBS-7  | 7.9308 | 40.0 | 16.1 | 9.83  | 85p                | **85p** |
| SBS-11 | 9.889  | 42.0 | 13.3 | 13.47 | 0p, 100p           | **100p** |

"Fields with data" = settings present in `DB/good_runList_GMn_nTPE_{LD2,LH2}_pass_3.csv`.
MC configs exist in addition for 0p of SBS-9 and SBS-7.

- Only the **LD2** data of the production settings enter the GMn result.
- **SBS-8** is formally an nTPE point, not GMn; its LD2 data are analyzed mainly as a check.
  Its **LH2** data (four field settings) are very important, in particular for the HCAL
  neutron detection efficiency (NDE) analysis.
- **SBS-11**: a shield beside the scattering chamber was in place for part of the running.
  The SBS-11 MC is a charge-weighted mix of Al-shield and no-shield samples (e.g.
  `simu_mix_prefixes` / `simu_mix_weights` = 0.607 / 0.393 in
  `scripts/configs/sbs11-sbs100p/conf_fit_W2.json`).
- Config `4363` (E_beam 6.373 GeV) in `ExpConstants.cpp` and `scripts/configs/sbs4363-sbs100p/`
  is a GEn setting, added at the very beginning only to verify the analysis code. It is not
  used again but is kept (it does no harm); listed in [legacy](howto/legacy.md).

## Conventions

### Directories and how scripts are run
- `include/`, `src/`: the `gmnana` shared library. `scripts/`: data analysis and data/MC
  comparison macros. `scripts/simulation/`: MC macros. `tests/`: unit tests. `DB/`: run lists
  and FF tables. `dflay/`: helper code incl. the JSON reader. Details: [01 setup](howto/01_setup.md).
- Scripts are ROOT macros run with `root` (local) or the Hall A `analyzer` (ifarm), with the
  library loaded by the installed `rootlogon.C`. Data scripts run from `scripts/`, MC scripts
  from `scripts/simulation/`; all paths inside the macros (`pdout/`, `siout/`,
  `simulation/siout/`, `../DB`) are relative to that working directory.

### Config files
- Data: `scripts/configs/sbs<K>-sbs<M>p/conf_<script>.json`
  (e.g. `scripts/configs/sbs4-sbs30p/conf_elas_ana_data.json`).
- MC: `scripts/simulation/configs/sbs<K>-sbs<M>p-simu/conf_<script>.json`.
- Cross-setting configs sit directly in `scripts/configs/` (e.g. `conf_w2_offset_check.json`,
  `conf_pDE_map.json`) or `scripts/simulation/configs/` (`conf_pDE_map_simu*.json`).
- Each file holds one JSON object whose top-level keys are **blocks**. A script reads one block,
  chosen from its arguments:
  - `elas_ana_data`, `qelas_ana_data`: `pass<P>_model<m>` (e.g. `pass3_model1`);
  - `elas_ana_simu`, `qelas_ana_simu`: `<generator>_<process>_model<m>[_<variant>]` (e.g.
    `simc_deeN_model2`, `g4sbs_inel_model2`, `simc_heep_model1`). The block without a variant
    is the production sample. A non-default sample of the same process is selected with the
    optional last macro argument `variant` (e.g. `noshield` -> `simc_deeN_model2_noshield`
    for SBS-11; `sbs_OFF` -> the SBS-9/SBS-14 field-off LH2 sample);
  - `fit_W2`, `fit_dx`, `fit_dy`: `elas` or `qelas`;
  - `pDE_*`: `prod`; `w2_offset_check`: `common`, `settings`, one block per setting.
  Other blocks (`*_thesis`, `*_aps25`, `*_mcp3`, `*_prod`, ...) are records of earlier runs.
  The MC macros can read any of them through `variant` (the suffix after the standard key).
  The other scripts (e.g. `fit_*`, `pDE_*`) have no such argument: to rerun an old block,
  rename it to the key the script expects.
- Keys such as `NOTE`, `NOTEs`, `README`, `COMMENTS` are free-text notes, ignored by the code.
  Free text **after** the closing `}` is also allowed: `JSONManager` reads only the first JSON
  object in the file.

### Output file names
Built from a config prefix (`outfile_prefix`, `output_filebase`) plus fixed parts:

| Script | Output (relative to its working directory) |
|---|---|
| `elas_ana_data.cpp`  | `pdout/<prefix>_elas_ana_data_sbs<K>_sbs<M>p_model<m>_pass<P>.root` (+ same-name `.pdf`) |
| `qelas_ana_data.cpp` | `pdout/<prefix>_qelas_ana_data_sbs<K>_sbs<M>p_model<m>_pass<P>.root` (+ `.pdf`) |
| `elas_ana_simu.cpp`  | `siout/<prefix>_elas_ana_<generator>_sbs<K>_sbs<M>p_model<m>.root` |
| `qelas_ana_simu.cpp` | `siout/<prefix>_qelas_ana_<generator>_sbs<K>_sbs<M>p_model<m>.root` |
| `fit_dx.cpp`         | `pdout/fits/sbs<K>sbs<M>p/<output_filebase>_fit_dx_<elas/qelas>_pass<P>_<gen>_sbs<K>_sbs<M>p_model<m>[_normweight][_noefficorr][_nonblkHCALcut].{root,pdf,csv}` (`.../sysstdy/` for cut-variation runs) |

The prefix records the analysis variant, e.g. `bsf_0p65zoff` (HCAL z offset 0.65 m) for data
or `bsf_mcp4_OT0_0p97sf` (MC pass 4, SBS field scale 0.97) for MC; `noshield_` marks SBS-11
no-shield MC. The empty-prefix case drops the leading `_`. `sbs-1p` in a data file name means
`SBS_magnet_percent = -1`, i.e. all field settings of that config together; the file is made with `hadd` of the
per-field outputs.

### Passes and models
- **Replay pass**: data pass 3 is current (run lists for passes 1-3 in `DB/`); MC pass 4 is
  current (`mcp4` in prefixes; `*_mcp3` blocks are the older MC pass).
- **Model** (how the expected nucleon position in HCAL is computed): 0 = from the reconstructed
  e- momentum, 1 = from the reconstructed e- angles (LH2 elastics), 2 = 4-vector calculation
  (LD2 QE).

### Where the analysis runs
- **ifarm / batch (JLab)**: the event-level steps (03, 04) read the large replayed files and
  run as Slurm jobs via `run-ana-data.sh` / `submit-run-ana-data.sh` and
  `simulation/run-ana-simu.sh` / `simulation/submit-run-ana-simu.sh`, which load the
  `analyzer` module and call the macro (`export GMN_ANA=... ROOTLOGON_ENV=...` before
  submitting, see [01 setup](howto/01_setup.md)).
- **Local**: steps 05-07 read only the `pdout/` / `siout/` output trees and run on a laptop with
  ROOT once those outputs are copied to the local data disk.

### Outputs and plots
- `scripts/pdout` (data outputs) and `scripts/simulation/siout` (MC outputs) are normally
  symlinks to a data disk; they are not in git.
- Every study writes into `pdout/` or `siout/` (sub-directories allowed, e.g. `pdout/fits/`,
  `pdout/pDE/`, `siout/w2_offset_check/`) and puts its plots into **one multi-page PDF** per
  run, not many single images.

## Where records live
- `CHANGELOG.md` (repo root): every change, dated, with files and reason.
- The MC repo's consistency audit (task specs, open-items list, audit report): data/MC
  consistency decisions (e.g. W2_offset method, SBS-11 shield mix). Decisions recorded there
  are treated as fixed in this repo.
- `scripts/lta/results.csv`: final per-setting R_n/p with stat and systematic errors, the
  input to the GMn extraction. It is intentionally **not tracked** in git yet (nor are the
  other `scripts/lta/` result files): it will still change, and it should eventually be
  produced directly by this repository.
- `scripts/gmn_extraction/README.md`: GMn extraction, FF models, error propagation.
- `DB/FF_fits/Ye2018/README.md`: provenance and caveats of the Ye FF tables.
- `NOTE`/`README` keys in config blocks: why particular cuts/offsets were chosen.

## How the stage files are organized

| File | Content |
|---|---|
| [howto/01_setup.md](howto/01_setup.md) | build, running macros, tests, library modules, `DB/`, `dflay/`, output locations |
| [howto/02_run_info.md](howto/02_run_info.md) | beam energy, charge, live time, run lists |
| [howto/03_data_analysis.md](howto/03_data_analysis.md) | LH2 elastic and LD2 QE analysis of data |
| [howto/04_mc_analysis.md](howto/04_mc_analysis.md) | MC analysis, g4sbs macros, momentum calibration |
| [howto/05_checks.md](howto/05_checks.md) | W2 offset, energy loss, detection-efficiency checks |
| [howto/06_fits.md](howto/06_fits.md) | W2, dx, dy fits |
| [howto/07_results.md](howto/07_results.md) | results table, GMn extraction, flavor decomposition |
| [howto/legacy.md](howto/legacy.md) | obsolete scripts (not maintained) |

Each script is documented by one card:

```markdown
### `path/to/script.cpp`
**What it does:** 2-4 lines in physics terms.
**Run:** example command (from scripts/ unless otherwise), with arguments explained.
**Inputs:** config file + the important keys; input files.
**Outputs:** files and where they go.
**Feeds:** which later step uses the output.
**Notes:** caveats, gotchas (optional).
```

## Script index

Paths relative to `scripts/` unless they start with `DB/`.

| Script | Purpose | Stage |
|---|---|---|
| `epics/get_avgHALLA_p_prun.cpp` | dump EPICS beam energy `HALLA_p` per EPICS event, per run | [02](howto/02_run_info.md) |
| `epics/get_avgHALLA_p_prun.ipynb` | per-run mean/std beam energy from that dump | [02](howto/02_run_info.md) |
| `epics/sbs11_ld2_ebeam_ana.ipynb` | first exploratory `HALLA_p` study (SBS-11 LD2) | [02](howto/02_run_info.md) |
| `epics/get_scalerdata_prun.cpp` | per-run BCM charge and DAQ live time from the scaler tree | [02](howto/02_run_info.md) |
| `epics/run_eORs_ana.sh` | runs the scaler macro for one or all configs | [02](howto/02_run_info.md) |
| `epics/get_charge_from_Bobs_script.sh` | per-run charge from the Hall A online charge tool | [02](howto/02_run_info.md) |
| `epics/bcmCheckCharge.cpp` | BCM current / charge / live-time check (dflay BCM tools) | [02](howto/02_run_info.md) |
| `epics/grab_epicsPV_from_evio.py`, `epics/jcache_0th_seg.sh` | read EPICS PVs (e.g. magnet currents) from raw EVIO | [02](howto/02_run_info.md) |
| `DB/add_scaler_files.ipynb`, `DB/join_n_merge_rLists.ipynb`, `DB/get_run_list.ipynb` | merge charge / beam energy into the run lists; print run numbers | [02](howto/02_run_info.md) |
| `elas_ana_data.cpp` | LH2 elastic analysis of data -> `Tout` tree + plots | [03](howto/03_data_analysis.md) |
| `qelas_ana_data.cpp` | LD2 (or dummy) QE analysis of data -> `Tout` tree + plots | [03](howto/03_data_analysis.md) |
| `run-ana-data.sh`, `submit-run-ana-data.sh` | Slurm wrappers for the two data macros | [03](howto/03_data_analysis.md) |
| `plot_gcuts_data.cpp` | plots of the variables used in the global cuts | [03](howto/03_data_analysis.md) |
| `simulation/elas_ana_simu.cpp` | LH2 elastic analysis of MC (SIMC/g4sbs), cross-section weighted | [04](howto/04_mc_analysis.md) |
| `simulation/qelas_ana_simu.cpp` | LD2 QE / inelastic analysis of MC, cross-section weighted | [04](howto/04_mc_analysis.md) |
| `simulation/run-ana-simu.sh`, `simulation/submit-run-ana-simu.sh` | Slurm wrappers for the MC macros | [04](howto/04_mc_analysis.md) |
| `simulation/g4mac/*` | g4sbs macros, SIMC input example, field-scan macro generator, `get_SBSconf_info.cpp` | [04](howto/04_mc_analysis.md) |
| `simulation/hcal/gen_xypos_DB.cpp` | HCAL block positions for the MC database | [04](howto/04_mc_analysis.md) |
| `pcalib/` (`MakeMomentumCalibrationTree.C`, setup file) | BB momentum calibration (real data) | [04](howto/04_mc_analysis.md) |
| `simulation/pcalib/` (`momentum_calib_simu.cpp`, setup files) | BB momentum/optics calibration on MC pseudo-data | [04](howto/04_mc_analysis.md) |
| `w2_offset_check.cpp` | LH2 data/MC check of the MC `W2_offset`, all settings | [05](howto/05_checks.md) |
| `simulation/eloss_check_simu.cpp` | how the MC momentum reconstruction treats the outgoing e- energy loss | [05](howto/05_checks.md) |
| `pDE_data.cpp`, `simulation/pDE_simu.cpp` | HCAL proton detection efficiency vs expected position (LH2), data / MC | [05](howto/05_checks.md) |
| `pDE_map.cpp`, `simulation/pDE_map_simu.cpp` | 2D HCAL proton detection efficiency maps, data / MC | [05](howto/05_checks.md) |
| `fit_W2.cpp` | W2 fit: MC signal (optionally a weighted mix) + background to data | [06](howto/06_fits.md) |
| `fit_dx.cpp` | HCAL dx fit: MC p and n shapes + background -> R_fit, systematic variations | [06](howto/06_fits.md) |
| `fit_dy.cpp` | HCAL dy data/MC fit (occasionally useful) | [06](howto/06_fits.md) |
| `lta/results.csv` | final per-setting R_n/p and systematics | [07](howto/07_results.md) |
| `lta/get_results.ipynb` | builds the results table from the `fit_dx` CSV outputs | [07](howto/07_results.md) |
| `lta/grab_rnpsf.ipynb` | R_fit tables for the fit variants across settings | [07](howto/07_results.md) |
| `lta/get_ff.cpp` | prints one FF at a Q² for several parametrizations | [07](howto/07_results.md) |
| `lta/customize_plots.cpp` | publication/thesis versions of plots | [07](howto/07_results.md) |
| `gmn_extraction/extractGMn.cpp`, `gmn_extraction/extractGMn_csv.cpp` | GMn from R for all FF models, full errors | [07](howto/07_results.md) |
| `flavordecomp/flavordecomp.cpp` | u/d quark flavor decomposition of the nucleon FFs | [07](howto/07_results.md) |

Obsolete scripts (e.g. the old HCAL efficiency studies `hcal_det_effi_pblk.cpp`,
`simulation/hcal_det_effi_simu.cpp`, `simulation/hcal_det_effi_gun.cpp`):
[legacy](howto/legacy.md). `thesis/`: thesis figure material, kept for
reference.
