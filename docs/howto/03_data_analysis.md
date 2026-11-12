# Stage 3 — Data analysis (LH2 elastics, LD2 quasi-elastics)

Back to the overview: [../ANALYSIS_GUIDE.md](../ANALYSIS_GUIDE.md) ·
Previous: [./02_run_info.md](./02_run_info.md) · Next: [./04_mc_analysis.md](./04_mc_analysis.md)

## Purpose of this stage

This stage turns the replayed real-data ROOT files (tree `T`) into compact per-setting analysis
trees (`Tout`). For every event that passes a loose global cut, the scripts rebuild the electron
kinematics from the BigBite track (Q², W², ε, with a target energy-loss correction applied to the beam
energy and to the track momentum). From that they predict where the struck nucleon should land on
HCAL and store the HCAL residuals `dx`, `dy`, the angles `thpq_p`/`thpq_n`, the coincidence time, the
fiducial-cut flags and the per-run normalization (charge, live time). Two scripts do this:
`elas_ana_data.cpp` for LH2 H(e,e'p) elastics (calibration, W² offset, proton detection efficiency)
and `qelas_ana_data.cpp` for LD2 D(e,e'N) quasi-elastics (the GMn yields). Every later stage reads
these `Tout` trees. It does not go back to the replayed files.

## Typical workflow

1. Make sure the run lists and per-run info exist in `DB/`
   (`good_runList_GMn_nTPE_<LH2|LD2|Dummy>_pass_<N>.csv`, see [./02_run_info.md](./02_run_info.md)).
   If `get_scaler_info` = 1, also make sure the scaler file `scripts/epics/epout/scalerdata_prun_SBS<conf>_<target>.root` exists.
   Note: LD2 run 11597 (SBS-4 50p) was added to the pass-3 list later and has charge/live time `-1`
   (not computed yet). This is fine because the charge-normalized yield is not used.
2. Optionally run `plot_gcuts_data.cpp` on a few runs to look at the distributions used in `global_cut`
   (vertex, PS energy, E/p, track momentum, W²) before you choose the cut values.
3. For each setting, in the order SBS-4 (0p, 30p, 50p), SBS-9 70p, SBS-8 (0p, 50p, 70p, 100p),
   SBS-14 (0p, 70p), SBS-7 85p, SBS-11 (0p, 100p), run `elas_ana_data.cpp` on LH2 (`model` 1).
   Fit the proton `dx`/`dy` peaks (the fit results are printed and shown on the summary page) and copy
   the means and sigmas into `dx_p_cut`/`dy_p_cut` and the other cut keys of the config. Then run it again.
   The LH2 outputs feed the W² offset check, the HCAL position offsets and the proton efficiency maps
   ([./05_checks.md](./05_checks.md)).
4. For the production settings (SBS-4 30p/50p, SBS-9 70p, SBS-8 70p, SBS-14 70p, SBS-7 85p,
   SBS-11 100p) run `qelas_ana_data.cpp` on LD2 (`model` 2). Run it on `Dummy` too if you need the
   target-wall (end-cap) background.
5. On the JLab ifarm, steps 3–4 are normally submitted as batch jobs with `submit-run-ana-data.sh`.
6. Pass the `pdout/*.root` files to the fits ([./06_fits.md](./06_fits.md)).

Configs live in `scripts/configs/sbs<conf>-sbs<field>p/` (e.g. `configs/sbs8-sbs70p/conf_elas_ana_data.json`).
Existing data configs, in Q² order:

| Setting | `conf_elas_ana_data.json` (LH2) | `conf_qelas_ana_data.json` (LD2) |
|---|---|---|
| SBS-4 | 0p, 30p, 50p | 0p, 30p, 50p |
| SBS-9 | 70p | 70p |
| SBS-8 | 0p, 50p, 70p, 100p | 0p, 50p, 70p, 100p |
| SBS-14 | 0p, 70p | 70p |
| SBS-7 | 85p | 85p |
| SBS-11 | 0p, 100p | 100p |

(`configs/sbs4363-sbs100p/` also exists, with only a qelas config. It is a GEn setting
(E_beam 6.373 GeV), added at the very beginning only to verify the code. It is not used again but
is kept.)

### Config file layout (common to `elas_ana_data` / `qelas_ana_data`)

Each JSON file holds several blocks named `pass<P>_model<M>` (e.g. `pass3_model1`). The script picks
the block from its `pass` and `model` arguments. The current production replay is pass 3. Older blocks
(pass 2, `_old`) are kept for reference. Main keys:

| Key | Meaning |
|---|---|
| `Ntype` | Struck-nucleon mass hypothesis: `p` (LH2), `np` = average n/p mass (LD2). |
| `SBS_config`, `SBS_magnet_percent` | Setting, e.g. 8 and 70. Field in % of 2100 A. |
| `runsheet_dir` | Directory of the run lists (`../DB`). The script reads `good_runList_GMn_nTPE_<target>_pass_<pass>.csv` and keeps the runs whose config and SBS magnet current match. |
| `Nruns_to_ana` | Number of runs to analyze (−1 = all). |
| `rootfile_dir` | Directory of the replayed data files. All files matching `*<run>*` are chained for each run. |
| `get_scaler_info` | 1 = attach the nearest scaler event (BCM counts/current) to each event. |
| `global_cut` | Cut on raw-tree branches (`bb.tr.vz[0]`, `bb.tr.p[0]`, `bb.ps.e`, `e.kine.W2`, ...) applied first. Only events that pass are written to `Tout`. |
| `W2_cutR` | [min, max] of the W² window that sets the `W2Cut` flag (GeV²). |
| `bbfidu_cutR` | [centre, half-width] of the BigBite fiducial cut `|fpX − 0.9·fpTh − c| ≤ w` (sets `bbfiduCut`). |
| `dx_p_cut`, `dy_p_cut` | Proton spot [mean, sigma, N sigma] in `dx`/`dy` (m). They define `pCut` and the HCAL safety margin. |
| `dx_n_cut`, `dy_n_cut` | Neutron spot, same format (qelas only, defines `nCut`). |
| `hcal_voffset`, `hcal_hoffset`, `hcal_zoffset` | HCAL position offsets (vertical, horizontal, along the HCAL axis, m) added to the nominal HCAL position when the expected hit position is computed. |
| `hcal_acl_ON` | 1 = offline choice of the "best" HCAL cluster (see `Tout` notes below). |
| `hcal_samp_frac_cutR` | Minimum HCAL sampling fraction (cluster E / ν), used only for the smallest-θpq cluster choice. |
| `max_N_tracks`, `max_N_HCAL_clusters` | Array sizes. |
| `coinT_ADC_cutR` | [mean, min, max] of the HCAL−SH ADC coincidence time (ns). The mean centres `coinT_ADC_c`. [min, max] define "in time" for the cluster choice. |
| `coinT_ADC_TOF_offset` | Offset that centres the TOF-corrected coincidence time. |
| `h_dxHCAL_lims`, `h_dyHCAL_lims`, `h_dxHCAL_p_fitR`, `h_dxHCAL_n_fitR`, `h_dyHCAL_fitR`, `h_coinT_ADC_fitR` | Histogram ranges and fit windows for the diagnostic plots only. |
| `outfile_prefix` | Prefix of the output file name (e.g. `bsf_0p77zoff`). |
| `verbose`, `verbose_function` | Verbosity. 0/0 = production. |

The meaning of every key is also written at the bottom of each `.cpp` file.

### The `Tout` tree (what the outputs hold)

One entry per event that passes `global_cut`. The real-data replay writes only events with a good
BigBite track, so every data event has a track and the `[0]` (best) track is always valid. Main branches:

| Group | Branches | Physics meaning |
|---|---|---|
| Run info | `rnum`, `segnum`, `gevnum` | Run, file segment and event number. |
| Beam / normalization | `ebeam`, `ebeam_corr`, `ebeam_std` | Run-averaged beam energy from the run list. `ebeam_corr` = after the energy loss in the target before the vertex. |
| | `dnewcharge`, `weight_norm` | Run charge (C). `weight_norm` = 1/(charge·DAQ live time), the per-event weight that gives a yield per coulomb. Not used by the production analysis. Run 11597 (SBS-4 50p LD2) has charge/live time `-1` in the pass-3 list (not computed yet, see [./02_run_info.md](./02_run_info.md)). |
| | `dnewcnt`, `dnewcurr` | Nearest scaler BCM count and current (only if `get_scaler_info`). |
| Electron kinematics | `Q2`, `W2`, `W`, `nu`, `epsilon`, `epsilon_p` (+`epsilon_n` in qelas) | Reconstructed with the chosen `model`: 0 = from the momentum, 1 = from the electron angle (elastic formula), 2 = full 4-vector. In data, `W2` has no offset. |
| | `etheta`, `ephi`, `pelas`, `dpel` | Electron angles, the elastic momentum expected at that angle, and `dpel = p/p_elastic − 1`. |
| | `pN_exp`, `thN_exp` | Expected recoil nucleon momentum and polar angle (direction of **q**). |
| Track | `vz`, `trP`, `trP_corr`, `trX/Y/Th/Ph`, `tgX/Y/Th/Ph`, `fpX/Y/Th/Ph`, `ethbend` | Vertex z, momentum (raw and energy-loss corrected), track at detector / target / focal plane, bend angle. |
| | `nhitsGEM`, `ngoodhitsGEM`, `trchi2ndf` | GEM track quality. |
| BigBite calorimeter | `ePS`, `eSH`, `EovP`, `atimeSH`, cluster block info | Preshower/shower energies. `EovP` = (PS+SH)/p. |
| HCAL (chosen cluster) | `eHCAL`, `xHCAL`, `yHCAL`, `nblkHCAL`, `rblkHCAL`, `cblkHCAL`, `idblkHCAL`, `atimeHCAL`, `tdcHCAL` | HCAL cluster used for the event (−99 if HCAL has no cluster). |
| | `xHCAL_exp`, `yHCAL_exp` | Where **q** (straight line from the vertex) hits the HCAL face. |
| | `p_def`, `xHCAL_exp_p` | Expected proton deflection in the SBS dipole (from ∫B·dl and the nucleon momentum). `xHCAL_exp_p = xHCAL_exp − p_def` (`elas_ana_data` only; `qelas_ana_data` has no `xHCAL_exp_p` branch, use `xHCAL_exp − p_def`). |
| | `dx`, `dy` | `xHCAL − xHCAL_exp`, `yHCAL − yHCAL_exp` (m). Neutrons peak near dx ≈ 0, protons near dx ≈ −p_def. |
| | `thpq_p`, `thpq_n` | Angle between **q** and the vertex→HCAL-cluster direction. `thpq_n` assumes no deflection. `thpq_p` first moves the cluster by the proton deflection. |
| Coincidence time | `coinT_ADC_c` | HCAL − SH ADC time minus the config mean (ns). Centred on 0 for true coincidences. |
| | `coinT_ADC_TOF_c`, `pTOF`, `TOF_corr`, `beta_p` (`nTOF`, `beta_n` in qelas) | Coincidence time after the nucleon time-of-flight correction, and its inputs. |
| | `coinT_trig`, `bbT_trig` | Trigger-TDC coincidence time (not used as a cut: it loses elastics). |
| Cut flags | `W2Cut` | `W2` inside `W2_cutR`. |
| | `bbfiduCut` | BigBite fiducial cut. |
| | `pCut` (`nCut`) | Inside the proton (neutron) dx–dy ellipse. `pdx_nS`, `ndx_nS`, `dy_nS` = distance from the peak in sigmas. |
| | `ARCut`, `SMCut`, `fiduCut` | HCAL active area (one block in from each edge), safety margin on the expected position (`fiduCut = ARCut && SMCut`). `SMx_nS_p`, `SMy_nS` = distance of the expected position from the active-area edge in sigmas, for a flexible margin. |
| All HCAL clusters (if `hcal_acl_ON`) | `*_acl[]`, `idclHCAL_htote`, `idclHCAL_intime`, `idclHCAL_sthpq_p` (`idclHCAL_sthpq` in qelas) | Per-cluster arrays and the index chosen by each algorithm. |
| GRINCH (SBS-8 and higher configs) | `clsizeGRINCH`, `cltmeanGRINCH`, ... | Cherenkov cluster info. |
| Cross sections (qelas only) | `sigMott`, `sigRed_p`, `sigRed_n`, `sigBorn_ratio` | Model cross sections at the event kinematics (Kelly / Seamus / Ye fits). |

HCAL cluster choice with `hcal_acl_ON` = 1: the clusters are sorted by total energy. The chosen cluster
is the highest-energy one that is **in time** (HCAL−SH ADC time inside `coinT_ADC_cutR`). If none is in
time, the highest-energy cluster is used. `idclHCAL_sthpq_p` also stores the in-time cluster with the
smallest `thpq_p` for alternative studies.

---

### `scripts/elas_ana_data.cpp`
**What it does:** Elastic H(e,e'p) analysis of LH2 data for one setting. It rebuilds the electron
kinematics with a target energy-loss correction, predicts the proton position on HCAL, picks the best
HCAL cluster and writes the `Tout` tree described above. It also fills diagnostic histograms (W, dx, dy,
coincidence time, HCAL hit maps) and fits the proton dx/dy peaks and the W² peak.
**Run:** from `scripts/`, with the gmn_ana library loaded by the rootlogon (see [./01_setup.md](./01_setup.md)):
```
analyzer -b -q <rootlogon.C> 'elas_ana_data.cpp("configs/sbs8-sbs70p/conf_elas_ana_data.json",3,1)'
```
Arguments: `configfilename`, `pass` (replay pass, default 3), `model` (0/1/2, default 1. LH2 production uses 1).
The config block `pass<pass>_model<model>` is used. `root -l -b -q` with the library loaded works the same way.
**Inputs:** config `configs/sbs<conf>-sbs<field>p/conf_elas_ana_data.json` (keys above). Run list
`DB/good_runList_GMn_nTPE_LH2_pass_<pass>.csv`. Replayed files in `rootfile_dir`. If `get_scaler_info`:
`epics/epout/scalerdata_prun_SBS<conf>_LH2.root`.
**Outputs:** `pdout/<outfile_prefix>_elas_ana_data_sbs<conf>_sbs<field>p_model<model>_pass<pass>.root`
(`Tout` + histograms + canvases) and a `.pdf` with the same base name (hit maps, dx/dy fits,
coincidence time, W², summary page with cuts, charge and fit results). The fit numbers are also printed.
**Feeds:** `dx_p_cut`/`dy_p_cut` and the HCAL offsets in the configs (iterate). `w2_offset_check.cpp`,
`pDE_data.cpp`, `pDE_map.cpp` ([./05_checks.md](./05_checks.md)). `fit_W2.cpp` and the LH2 fits
([./06_fits.md](./06_fits.md)).
**Notes:** `pdout/` must exist (on some machines it is a link to a data disk). The cut flags are only
stored, not applied: `Tout` keeps every global-cut event, so downstream scripts choose their own cuts.
The data files for SBS-8 0p/50p/100p that the checks use are one combined file named `sbs8_sbs-1p`
(see [./05_checks.md](./05_checks.md)). How it is produced is not visible in the current code (verify).

### `scripts/qelas_ana_data.cpp`
**What it does:** Quasi-elastic D(e,e'N) analysis of LD2 (or Dummy) data. It is the same machinery as
`elas_ana_data.cpp`, but with the n/p-averaged nucleon mass, both proton and neutron spot cuts (`pCut`,
`nCut`), a fiducial safety margin built from the average n/p dx width, neutron TOF, and model
cross-section branches. These trees are the input for the n/p yield ratio, from which GMn is extracted.
**Run:** from `scripts/`:
```
analyzer -b -q <rootlogon.C> 'qelas_ana_data.cpp("configs/sbs8-sbs70p/conf_qelas_ana_data.json","LD2",3,2)'
```
Arguments: `configfilename`, `target` (`"LD2"` or `"Dummy"`), `pass` (default 3), `model` (default 2.
LD2 production uses 2).
**Inputs:** `configs/sbs<conf>-sbs<field>p/conf_qelas_ana_data.json` (as above, plus `dx_n_cut`,
`dy_n_cut`, `h_dxHCAL_n_fitR`). `DB/good_runList_GMn_nTPE_<LD2|Dummy>_pass_<pass>.csv`. Replayed files.
For `Dummy`, the string `LD2` in `rootfile_dir` is replaced by `Dummy`.
**Outputs:** `pdout/<outfile_prefix>_qelas_ana_data_sbs<conf>_sbs<field>p_model<model>_pass<pass>.root` + `.pdf`.
For Dummy, `dummy_` is added to the prefix.
**Feeds:** the dx fits / n:p ratio ([./06_fits.md](./06_fits.md)) and then the GMn extraction ([./07_results.md](./07_results.md)).
**Notes:** For Dummy runs the beam energy is the nominal setting value, not the run-list value. Only the
LD2 data of the production settings enter the GMn results. SBS-8 LD2 is mainly a cross-check.

### `scripts/run-ana-data.sh`
**What it does:** Batch wrapper for one data job on the JLab ifarm. It loads the Hall A modules
(`analyzer`), goes to `$GMN_ANA/scripts` and runs `elas_ana_data.cpp` or `qelas_ana_data.cpp`.
**Run:** `sbatch run-ana-data.sh <elas|qelas> <conf> <field> <pass> <model> [<LD2|Dummy>]`,
e.g. `sbatch run-ana-data.sh qelas 8 70 3 2 LD2`. Normally called through `submit-run-ana-data.sh`.
**Inputs:** the config `configs/sbs<conf>-sbs<field>p/conf_<process>_ana_data.json` of that setting.
The environment variables `GMN_ANA` (repository path) and `ROOTLOGON_ENV` (installed `rootlogon.C`):
`export GMN_ANA=... ROOTLOGON_ENV=...` before submitting (`sbatch` passes the environment to the job).
If unset, the original author's ifarm paths are used.
**Outputs:** as for the macro it runs (`pdout/...`).
**Feeds:** —
**Notes:** The script exits with an error if `$GMN_ANA/scripts` or `$ROOTLOGON_ENV` is not found.

### `scripts/submit-run-ana-data.sh`
**What it does:** Submits `run-ana-data.sh` with `sbatch`, using a readable job name
(`data_<process>_sbs<conf>_f<field>p_<target>_pass<pass>_model<model>`) and 2 GB per CPU.
**Run:** from `scripts/`: `./submit-run-ana-data.sh elas 4 30 3 1` or `./submit-run-ana-data.sh qelas 4 30 3 2 LD2`.
For `elas` the target is set to LH2 automatically.
**Inputs:** same arguments as `run-ana-data.sh`.
**Outputs:** Slurm logs `slurm_out/<jobname>_<jobid>.out/.err` (the script creates `slurm_out/`).
Set `GMN_ANA` and `ROOTLOGON_ENV` first (see `run-ana-data.sh`).
**Feeds:** —

### `scripts/plot_gcuts_data.cpp`
**What it does:** Quick look at the variables used in `global_cut` (W, W², vertex z, track momentum
and angles, PS/SH energies, E/p, HCAL energy, trigger coincidence time, GEM hits). It copies them to a
flat tree with **no cuts**, so you can choose the cut values.
**Run:** from `scripts/`:
```
analyzer -b -q <rootlogon.C> 'plot_gcuts_data.cpp("configs/sbs8-sbs70p/conf_plot_gcuts_data.json","pdout/plot_gcuts_data")'
```
Arguments: `configfilename`, `filebase` (output path + base name, default `pdout/test_plot_gcuts_data`).
**Inputs:** `configs/sbs<conf>-sbs<field>p/conf_plot_gcuts_data.json` with **top-level** keys (no
pass/model blocks): `SBS_config`, `SBS_magnet_percent`, `target` (LH2/LD2), `replay_pass`,
`Nruns_to_ana`, `rootfile_dir_lh2`/`runnums_lh2`, `rootfile_dir_ld2`/`runnums_ld2` (explicit run lists,
not the DB).
**Outputs:** `<filebase>_sbs<conf>_sbs<field>p_pass<pass>_<target>.root` with a tree `Tout`. Draw the
distributions yourself in ROOT.
**Feeds:** your choice of `global_cut` in the elas/qelas configs.
**Notes:** The existing configs point to pass-1 directories and old run lists. Update them before you use them.
