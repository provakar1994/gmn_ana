# Stage 5 — Checks: W² offset, energy loss, HCAL detection efficiency

Back to the overview: [../ANALYSIS_GUIDE.md](../ANALYSIS_GUIDE.md) ·
Previous: [./04_mc_analysis.md](./04_mc_analysis.md) · Next: [./06_fits.md](./06_fits.md)

## Purpose of this stage

These scripts check that data and MC agree, before the MC is used as a template in the fits. They
read the `Tout` trees of stages 3 and 4 (or MC truth), not the raw replays, except where noted. They
answer three questions:
(1) **W² offset:** does the MC elastic W² peak sit where the data peak sits, i.e. is the
`W2_offset` in the MC configs right?
(2) **Energy loss:** how does the MC BigBite momentum reconstruction treat the scattered electron's
energy loss in the target, chamber and shield?
(3) **HCAL detection efficiency:** how efficient is HCAL for protons (measured with LH2 elastics) and
neutrons (MC), and how does this vary across the HCAL face? The proton efficiency map of the data is
also used to reweight the MC (`weight_effic`). SBS-8 LH2 data are the key input for this, which is why
SBS-8 matters even though it is not a GMn point.

## Typical workflow

1. After stage 3 (LH2 data) and stage 4 (LH2 `simc heep` MC) exist for all settings
   (SBS-4 0p/30p/50p, SBS-9 70p, SBS-8 0p/50p/70p/100p, SBS-14 70p, SBS-7 85p, SBS-11 100p), run
   `w2_offset_check.cpp`. For every setting marked `FLAG` (or `check`), put the new value into
   `W2_offset` of **all** MC blocks of that setting (LH2 and LD2) and rerun stage 4.
2. If the MC momentum scale is in question (e.g. after a new MC optics fit or a shield change), run
   `eloss_check_simu.cpp` on the affected samples.
3. Proton detection efficiency, data: run `pDE_map.cpp` on the LH2 `elas_ana_data` output (SBS-8, all
   fields combined, is the reference map). `pDE_data.cpp` is the simpler 1D version per setting.
4. Proton detection efficiency, MC: run `pDE_map_simu.cpp` (and/or `pDE_simu.cpp`) on the matching
   `elas_ana_simu` output and compare it with step 3.

The older HCAL efficiency studies (`hcal_det_effi_pblk.cpp`, `simulation/hcal_det_effi_simu.cpp`,
`simulation/hcal_det_effi_gun.cpp`) are obsolete; see [./legacy.md](./legacy.md).

---

## 5.1 W² offset

### `scripts/w2_offset_check.cpp`
**What it does:** Checks the MC `W2_offset` of every setting with LH2 strict elastics. It reads the
offset already in each MC file back (mean of `W2 − W·W`) and removes it, giving the zero-offset MC.
It applies identical cuts to data and MC: `thpq_p < thpq_p_cut`, `eHCAL > eHCAL_min`, `bbfiduCut`,
`|vz| < vz_cut`, and for data also a window around the coincidence-time peak. Then it measures the shift
needed to put the MC peak on the data peak, in two ways: an iterative core-Gaussian peak difference,
and a χ² scan of the MC template (data = N·MC(W² − δ) + polynomial, with an optional MC width stretch).
Best value = mean of the core result and the stretched template. Systematic = cut variations ⊕ fit-window
spread ⊕ method spread. Verdict per setting: `ok` / `check` / `FLAG`.
**Run:** from `scripts/`:
```
root -l -b -q 'w2_offset_check.cpp("configs/conf_w2_offset_check.json")'
root -l -b -q 'w2_offset_check.cpp("configs/conf_w2_offset_check.json","sbs8_70p")'   # one setting
```
Arguments: `configfilename`, optional setting name (as listed in `settings`).
**Inputs:** `scripts/configs/conf_w2_offset_check.json`:
- `common`: `data_dir` (`pdout`), `simu_dir` (`simulation/siout`), `out_dir`, `out_file_base`,
  `thpq_p_cut` (rad), `vz_cut` (m), `coinT_window` (ns).
- `settings`: the list of setting blocks, in Q² order (`sbs4_0p` … `sbs11_100p`).
- One block per setting: `data_file` (an `elas_ana_data` output), `data_run_cut` (a run selection,
  e.g. `rnum>=13459&&rnum<=13466`, for files that combine several field settings), `simu_files` +
  `simu_weights` (`elas_ana_simu` outputs and their mix weights), `eHCAL_min` (GeV).
- SBS-11: `simu_files` = Al-shield and no-shield samples with `simu_weights` `[0.607, 0.393]`. A single
  offset is determined for the mix.
- SBS-8 0p/50p/100p: the data come from the combined `..._sbs8_sbs-1p_...` file and are selected with
  `data_run_cut` (runs from `DB/good_runList_GMn_nTPE_LH2_pass_3.csv`).
**Outputs:** in `out_dir` (default `simulation/siout/w2_offset_check/`): `<out_file_base>.pdf` (3 pages
per setting + summary), `<out_file_base>.root` (the W² histograms) and `<out_file_base>.csv` (one row per
setting: baked offset, measured offset, errors, verdict). A single-setting run appends `_<setting>` to the names.
**Feeds:** `W2_offset` in all MC configs (`elas_ana_simu` and `qelas_ana_simu` blocks of the setting)
→ rerun stage 4 → fits.
**Notes:** The MC files are weighted with `weight × simu_weight`. The file names in the config must match
the current `outfile_prefix` of the stage 3/4 runs. The method replaces the earlier by-eye matching of the
LH2 W² peaks. LD2 uses the LH2 value of the same setting.

---

## 5.2 Energy loss

### `scripts/simulation/eloss_check_simu.cpp`
**What it does:** Checks how the MC BigBite momentum reconstruction handles the scattered electron's
energy loss on the way out. It histograms `dp = p_recon − p_true` (MeV) at the vertex, overall and versus
vz and electron angle, without cross-section weights (detector response only). It compares the core peak
with the per-event loss in the target + chamber (dT) and in the Al shield (dS), both computed with the
data energy-loss model (`util_pd::GetElossInTgt`). The printout gives the peak expected in two scenarios:
A (the MC optics already absorbs dT) and B (the optics absorbs nothing), for shield and no-shield samples.
**Run:** from `scripts/simulation/`:
```
analyzer -b -q <rootlogon.C> 'eloss_check_simu.cpp("configs/sbs14-sbs70p-simu/conf_qelas_ana_simu.json",2,"simc","deeN")'
```
Arguments: `configfilename`, `model` (default 2), `generator`, `process` (`heep` → LH2, otherwise LD2).
**Inputs:** the same config blocks as `elas_ana_simu` / `qelas_ana_simu` (only `SBS_config`,
`SBS_magnet_percent`, `rootfile_dir`, `prefix_to_filebase`, `Njobs_to_ana`, `global_cut`,
`max_N_tracks`, `outfile_prefix`, `verbose_function`). It reads the replayed MC files directly.
**Outputs:** `siout/<outfile_prefix>_eloss_check_simu_<gen>_<process>_sbs<conf>_sbs<field>p_model<model>.root`
+ `.pdf`, and a printed summary (core peak, mean, median, dT, dS, expected peaks).
**Feeds:** the decision on whether to apply an energy-loss correction to MC (the MC scripts currently
apply none) and the energy-loss part of the data/MC consistency audit (MC repo).
**Notes:** The script does not know whether a sample has the shield. Read the matching line of the printout.

---

## 5.3 HCAL proton / neutron detection efficiency

Definition used throughout: from LH2 elastics selected with **electron-arm cuts only** (global cut +
W² window, expected position inside the HCAL fiducial region), efficiency = fraction that also have a
matching HCAL cluster (`harm_cut`, e.g. `thpq_p<0.04&&eHCAL>…` or a dx/dy spot cut). It is studied
versus the **expected** HCAL position, so the cut on the denominator does not bias the result.

### `scripts/pDE_data.cpp`
**What it does:** 1D proton detection efficiency versus `xHCAL_exp` and `yHCAL_exp` for one setting,
from LH2 data (ratio of e-arm+HCAL cut events to e-arm-only events, fitted with a constant). It also
plots the dx–dy spots.
**Run:** from `scripts/`: `root -l -b -q 'pDE_data.cpp("configs/sbs8-sbs70p/conf_pDE_data.json")'`
**Inputs:** `configs/sbs<conf>-sbs<field>p/conf_pDE_data.json`, block `prod`: `SBS_config`,
`SBS_magnet_percent`, `model`, `pass`, `data_file_prefix` (selects
`pdout/<prefix>_elas_ana_data_sbs<conf>_sbs<field>p_model<model>_pass<pass>.root`), `global_cut`,
`w2_cut`, `harm_cut` (on `Tout` branches), fiducial settings (`apply_fidu_cut`, `AR_width_x_y`,
`sbs_kick`, `SM_width_xp_xn_y`), histogram ranges `h_xexp_lim`, `h_yexp_lim`, `h_w2_lim`.
Configs exist for SBS-4 30p/50p, SBS-9 70p, SBS-8 0p/70p.
**Outputs:** `pdout/pDE/<prefix>_pDE_data_sbs<conf>_sbs<field>p_model<model>_pass<pass>_1.png`, `_2.png`.
**Feeds:** quick look only. `pDE_map.cpp` is the full version.
**Notes:** The header comment ("stability of n/p ratio") is a leftover. The fit range is taken from the
`*_lim` keys (the `*_fitR` keys are not used). `pdout/pDE/` must exist.

### `scripts/pDE_map.cpp`
**What it does:** Full proton detection efficiency study from LH2 data. It makes 1D efficiencies vs
x/y expected and the **2D efficiency map** `h2_effi_map` versus (`yHCAL_exp`, `xHCAL_exp − p_def`),
with binomial errors. It also shows the dx–dy spots, the W² distributions and the cut envelope.
**Run:** from `scripts/`: `root -l -b -q 'pDE_map.cpp("configs/conf_pDE_map_pass3.json")'`
**Inputs:** `configs/conf_pDE_map_pass3.json` (or `conf_pDE_map.json` for pass 2). The script reads only
the block named **`prod`**. The file holds `prod_sbs4`, `prod_sbs9`, `prod_sbs8`, `prod_sbs14`,
`prod` (SBS-7), `prod_sbs11`, so rename the block you want to `prod`. Keys: `SBS_config`,
`SBS_magnet_percent` (`-1` = the combined all-field file), `model`, `pass`, `data_file_prefix`,
`outfile_prefix`, `global_cut`, `w2_cut`, `SM_cut_x`, `SM_cut_y` (expected-position window = fiducial
region), `harm_cut`, histogram ranges/fit ranges `h_xexp_*`, `h_yexp_*`, `h_w2_lim`.
**Outputs:** `pdout/pDE/<outfile_prefix><data_file_prefix>_pDE_data_sbs<conf>_sbs<field>p_model<model>_pass<pass>.root`
(with `h2_effi_map`), `.pdf` and `_1.._4.png`. Note that there is no `_` between `outfile_prefix` and `data_file_prefix`.
**Feeds:** `weight_effic` in `elas_ana_simu.cpp` / `qelas_ana_simu.cpp`. They read the SBS-8 LH2
all-field map from the file named by the constant `kHCALEffiMapFile`,
`../pDEmap/0p77zoff_pDE_data_sbs8_sbs-1p_model1_pass3.root` relative to `scripts/simulation/`, i.e.
`scripts/pDEmap/`. That copy is tracked in git so the MC weights are reproducible. It was made with the
`prod_sbs8` block of `configs/conf_pDE_map_pass3.json` with `data_file_prefix` = `bsf_0p77zoff` and an
empty `outfile_prefix`, from `pdout/bsf_0p77zoff_elas_ana_data_sbs8_sbs-1p_model1_pass3.root`. To
update it, rerun `pDE_map.cpp` and copy the output (renamed to the name above) into `scripts/pDEmap/`. Data/MC efficiency comparison with `pDE_map_simu.cpp`.
**Notes:** Produce the combined SBS-8 (and SBS-4) `sbs-1p` input file first. It is made with `hadd` of the
per-field `elas_ana_data` outputs. The fiducial safety-margin
widths drawn on the map are hard-coded (`SM_w{0.22,0.22,0.27}`).

### `scripts/simulation/pDE_simu.cpp`
**What it does:** MC version of `pDE_data.cpp`: 1D proton efficiency vs expected x/y from an
`elas_ana_simu` output.
**Run:** from `scripts/simulation/`: `root -l -b -q 'pDE_simu.cpp("configs/sbs8-sbs70p-simu/conf_pDE_simu.json")'`
**Inputs:** `configs/sbs<conf>-sbs<field>p-simu/conf_pDE_simu.json`, block `prod`: as `pDE_data`, but with
`generator` and `simu_file_prefix` (selects `siout/<prefix>_elas_ana_<gen>_sbs<conf>_sbs<field>p_model<model>.root`).
Configs exist for SBS-4 30p/50p, SBS-9 70p, SBS-8 0p/70p, SBS-14 70p, SBS-7 85p, SBS-11 100p.
**Outputs:** `siout/pDE/<prefix>_pDE_simu_<gen>_sbs<conf>_sbs<field>p_model<model>_*.png`.
**Feeds:** comparison with `pDE_data.cpp`.
**Notes:** `siout/pDE/` must exist. The numerator and denominator are both filled with `weight` (pure
MC efficiency, unlike `pDE_map_simu.cpp`).

### `scripts/simulation/pDE_map_simu.cpp`
**What it does:** MC version of `pDE_map.cpp`: weighted 1D efficiencies and the 2D map versus the
expected position (shifted by the proton deflection), with weighted binomial errors.
**Run:** from `scripts/simulation/`: `root -l -b -q 'pDE_map_simu.cpp("configs/conf_pDE_map_simu_pass4.json")'`
**Inputs:** `configs/conf_pDE_map_simu_pass4.json` (pass-4 MC. `_pass3`, plain = older). Only the block
named **`prod`** is used. The file has blocks for SBS-4, 9, 8, 14, 7, 11, including `neen` (neutron),
`force_n` and `sbs_OFF` (field off) variants. Keys as in `pDE_map.cpp`, with `generator` and
`simu_file_prefix` instead of the data keys.
**Outputs:** `siout/pDE/<outfile_prefix><simu_file_prefix>_pDE_simu_<gen>_sbs<conf>_sbs<field>p_model<model>.root`
+ `.pdf`/`.png`.
**Feeds:** data/MC comparison of the proton (and, with `neen`, neutron) detection efficiency → NDE
systematics ([./06_fits.md](./06_fits.md), [./07_results.md](./07_results.md)).
**Notes:** The denominator is weighted with `weight`, the numerator with `weight_effic` (which already
contains the data efficiency map). This is intentional: the ratio is not the pure MC efficiency but
shows the efficiency dip that the data map introduces into the MC. For the pure MC efficiency use
`pDE_simu.cpp`.
