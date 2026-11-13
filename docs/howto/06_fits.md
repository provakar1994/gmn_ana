# Stage 6: Fits (W², Δx, Δy)

[Overview](../ANALYSIS_GUIDE.md) · next stage: [07 Results](./07_results.md)

## Purpose of this stage

This stage compares the analyzed data trees with the matching MC trees, histogram against
histogram. The central fit is **`fit_dx.cpp` on LD2 quasi-elastic data**. It fits the HCAL
Δx distribution (Δx = observed minus expected hadron x position on HCAL) with the MC proton
shape, the MC neutron shape and a background shape. The fitted ratio
**R_fit = R''_n/p** is the data/MC scale factor for the n/p yield ratio. Multiplied by the
MC Born ratio R_MC (σ_n/σ_p from the SIMC form factors at ⟨Q²⟩, ⟨ε⟩), it gives the measured
n/p Born cross-section ratio R used for GMn ([07](./07_results.md)). `fit_W2.cpp` is a
supporting check. It fits the LH2 elastic W² peak with MC to set or confirm the MC
`W2_offset`, and it handles the SBS-11 shield-on/off MC mix. `fit_dy.cpp` is an occasional
quick look at Δy.

## Typical workflow

1. **Prerequisites.** You need the per-setting analysis trees (`Tout`) for data in
   `scripts/pdout/` and for MC (SIMC signal and g4sbs inelastic) in
   `scripts/simulation/siout/`. Earlier stages make them; see the [overview](../ANALYSIS_GUIDE.md).
   Load the gmn_ana library in ROOT (the `rootlogon.C` of your install's `runenv/`) and run
   everything from `scripts/`.
2. **(LH2, as needed) W² check:** run `fit_W2.cpp(...,1)` on the LH2 elastic config of a
   setting. Confirm that the MC W² peak position (`W2_peakOff_for_MC_sig`) matches the data.
   For SBS-11 100p, use the charge-weighted shield/no-shield MC mix. The systematic
   all-settings W2-offset check is done by `w2_offset_check.cpp` (earlier stage, see the overview).
3. **(LD2) Fix the MC Δx peak offsets:** run `fit_dx.cpp(...,0)` once with
   `vary_pn_xOff_ranges[0] = 1`, so the p and n peak offsets float within the given ranges.
   Copy the fitted offsets into `dx_offset_MC_for_p` / `dx_offset_MC_for_n`, set the flag
   back to 0 and rerun. (This is how the offsets are normally set; verify.)
4. **(LD2) Production fit:** run `fit_dx.cpp(conf,0)` for every production setting in Q²
   order: SBS-4 30p, SBS-4 50p, SBS-9 70p, SBS-8 70p (nTPE point, a check only),
   SBS-14 70p, SBS-7 85p, SBS-11 100p. Inspect the PDF (fit quality, residuals, χ²/NDF)
   and keep the CSV.
5. **Systematic variants:** rerun the same config with the switches flipped. Each variant
   gets its own output name. They are `RC_prescription` (production `0nm3` → variant `m3`,
   giving "OTm3"), `use_effi_corr = 0` (`_noefficorr`) and `apply_nblkHCAL>1_cut = 0`
   (`_nonblkHCALcut`). For cut-variation studies set `is_vary_cut = 1`; the output then goes
   to `.../sysstdy/`.
6. Collect the CSVs into the results table ([07](./07_results.md)).

---

### `scripts/fit_dx.cpp`

**What it does:** Fits the HCAL Δx distribution of data with MC template histograms, using
linear interpolation of the histograms (`TH1::Interpolate`) as the fit function.
- **LD2 / QE** (`is_elastic=0`) has the model
  `f(Δx) = N·[ h_p^MC(Δx) + R·h_n^MC(Δx) ] + bg(Δx)`. Here **R = R''_n/p** (shown as
  `R_{n/p}^{sf}`) is the data/MC n/p scale factor and N is a common normalization. The fit is
  done 8 times with different background models (below).
- **LH2 / elastic** (`is_elastic=1`) uses a single proton template plus background. This is
  a check, not used for GMn.

**Run** (from `scripts/`):
```
root -l -b -q 'fit_dx.cpp("configs/sbs4-sbs30p/conf_fit_dx.json", 0)'
```
- arg 1: config path. Only the block named `"qelas"` (for `0`) or `"elas"` (for `1`) is
  read. Other blocks in the file (`qelas_aps26`, `qelas_thesis`, `*_old`, ...) are archived
  variants; to use one, rename it to `qelas`.
- arg 2: `is_elastic` (default `1`). Use **`0`** for the GMn (LD2) fit.

**Inputs:** `configs/<setting>/conf_fit_dx.json`, block `qelas`. The important keys:

| Key | Meaning |
|---|---|
| `SBS_config`, `SBS_magnet_percent`, `model`, `pass`, `generator` | select the input files |
| `data_file_prefix` | data file `pdout/<prefix>_qelas_ana_data_sbs<conf>_sbs<mag>p_model<model>_pass<pass>.root` |
| `simu_file_prefix` | MC signal file `simulation/siout/<prefix>_qelas_ana_<generator>_sbs<conf>_sbs<mag>p_model<model>.root`; a `%s` in it is replaced by `RC_prescription` |
| `inel_file_prefix` | inelastic MC `simulation/siout/<prefix>_inel_qelas_ana_g4sbs_sbs<conf>_sbs<mag>p_model<model>.root` |
| `RC_prescription` | radiative-correction variant of the SIMC sample, filled into the file names. The config comment says production is `0nm3` ("Prod=>1,0nm3,1" means nblk cut 1, RC 0nm3, effi corr 1), and `m3` is the "OTm3" variant. `0nm3` = SIMC d(e,e'p) generated with `OT0` and d(e,e'n) with `OTm3` (the deen sample exists only as `OTm3`). Production configs are set to `0nm3`; set `m3` only for the RC systematic run |
| `use_effi_corr` | 1: weight MC with `weight_effic` (includes the HCAL efficiency correction); 0: plain `weight` |
| `use_norm_weight` | 1: use the charge-normalized `*_norm` weights (normally 0) |
| `apply_nblkHCAL>1_cut` | appends `nblkHCAL>1` to all signal/bg cuts (production: 1) |
| `cuts_for_signal_simu` | signal cuts. **Used for both data and MC** (identical by design) |
| `extra_data_signal_cut` | cut added to data only, e.g. the coincidence-time cut `abs(coinT_ADC_TOF_c)<=3.5` |
| `cuts_for_bg_data1` | "anti-Δy" data background sample (e.g. `abs(dy)>0.5`, other cuts as signal) |
| `cuts_for_bg_data2` | "anti-coincidence-time" (accidentals) data background sample |
| `dx_offset_MC_for_p`, `dx_offset_MC_for_n` (m) | shifts added to MC Δx for p and n before filling |
| `vary_pn_xOff_ranges` | `[flag, pmin, pmax, nmin, nmax]`; flag=1 lets the p/n offsets float within these limits |
| `Order_of_poly_bg` | order of the polynomial background in fit 5 (and the elastic-mode fit); fit 3 is always pol2 |
| `Par_guess_of_gaus_bg` | start values `[amp, mean, sigma]` of the Gaussian background (fit 4) |
| `h_dx` `[nbins,min,max]`, `dx_fit_range` | histogram binning and fit range |
| `SB_reject_points` | side-band exclusion window (elastic mode only) |
| `use_custom_fiduCut`, `AR_width_x_y`, `SM_width_xp_xn_y`, `sbs_kick` | optional custom HCAL active-area / safety-margin fiducial cut. **Retired**: off (0) in all production configs; the explicit `yHCAL_exp` and `xHCAL_exp(-p_def)` window terms in the cut strings replace it |
| `is_vary_cut`, `cut_vary_style` (0–5), `param_to_vary`, `param_to_vary_for_cs5`, `cut_iter_min_width`, `h_cut_param`, `fidu_vary_*`, `apply_to_data_only` | cut-variation study (see Notes) |
| `output_filebase` | output name stem (`%s` replaced by `RC_prescription`) |

`cuts_for_signal_data` and `cuts_for_bg_simu` are **ignored** by the current code (the
config says "hardcoded"). It also reads the run list from `../DB` for the Δx-vs-run plot.

**The 8 QE fits** (index k = canvas `ck` = CSV columns `Rk`, `Bk`, ...; names follow
`scripts/lta/bg_map.csv`):

| k | Background model | B column means |
|---|---|---|
| 0 | none | – |
| 1 | data histogram, anti-Δy (`cuts_for_bg_data1`) ("anti-dy") | scale of bg histo |
| 2 | inelastic MC, p+n summed ("inel-MC") | scale of bg histo |
| 3 | polynomial, order 2 = **production** (named constant `kProdBgPolyOrder` in `fit_dx.cpp`, not a config key) ("pol2") | constant term |
| 4 | Gaussian ("gaus") | amplitude |
| 5 | polynomial, order `Order_of_poly_bg` (3 in all production configs) ("pol3") | constant term |
| 6 | data histogram, anti-coincidence time (`cuts_for_bg_data2`) ("anti-dt") | scale of bg histo |
| 7 | B1·(inel MC p + n) + B2·(anti-dt data, pre-scaled by 1/3) ("inel-MC+rndm") | `B7` = B1; B2 in the extra columns `B7b,B7berr` at the end of the row |

With `vary_pn_xOff_ranges[0]=1`, the floating-offset variants of these fits are used
(fit 7 then drops the anti-dt term). Fits with histogram backgrounds use a weighted
log-likelihood (`"RWL"`); polynomial and Gaussian fits use χ² (`"R"`).

**R_MC:** the code fills MC histograms of the true Q², ε and θ_e (`vQ2`, `vepsilon`, `vetheta`,
weight `weight`) after the signal cuts. From ⟨Q²⟩ and ⟨ε⟩ it computes the effective
electron angle (`kine::etheta_from_epsilon`) and then
`RpMC = kine::sigmaBorn_ratio_MC(θ_eff, ⟨Q²⟩)`. That is σ_n/σ_p (Born, without the Mott
factor) with exactly the form factors that SIMC uses (Kelly GEp/GMp/GMn, Riordan GEn; see
[`scripts/gmn_extraction/README.md`](../../scripts/gmn_extraction/README.md)). The code also
writes `RMCnf`, the raw MC n/p yield ratio after cuts. It is for information only and is
not used for R.

**Outputs:** go to `pdout/fits/sbs<conf>sbs<mag>p/`, or `.../sysstdy/` when `is_vary_cut=1`.
The name is
`<output_filebase>_fit_dx_qelas_pass<pass>_<gen>_sbs<conf>_sbs<mag>p_model<model>[_normweight][_noefficorr][_nonblkHCALcut]`
followed by:
- `.csv`: one row per cut iteration (one row in production). Columns:
  `conf,field,cut,min,max`, then the MC kinematics `vethMean/Min/Max, vQ2Mean/Min/Max, vepsMean/Min/Max`,
  then `RMCnf,RMCnferr,RpMC`, then `chi20,NDF0,R0,R0err` and for k=1..7
  `chi2k,NDFk,Rk,Rkerr,Bk,Bkerr,Ypk,Ypkerr,Ynk,Ynkerr,Ybgk,Ybgkerr`, and last `B7b,B7berr`
  (fit 7's second background scale; added 2026-11-12, older CSVs end at `Ybg7err`). The Y columns are the
  fitted p, n and background yields with errors. A tracked example is
  `scripts/lta/mcp4_OT0nm3_fit_dx_qelas_pass3_simc_sbs4_sbs50p_model2.csv`.
- `.pdf`: all canvases. **c0–c7** are the 8 fits. Each has the data (points), the total fit,
  the scaled p and n MC shapes and the background, with a residual pad (c1–c7). In
  cut-variation mode each iteration adds a `cCut` page showing the cut region on the data
  and MC distributions of the varied variable.
- `_0.png` … `_7.png`: PNGs of c0–c7.
- `.root`: the canvases plus the input histograms (`h_dxHCAL_data`, `h_dxHCAL_simu_p/_n`,
  the bg histograms, `h_vQ2`, `h_vepsilon`, `h_vetheta`, `h2_dxHCAL_vs_rnum`). In
  cut-variation mode it also holds `hgist_cv_1..7` (R vs iteration for each bg model) and
  the canvas `cGist`.
- `.gif`: an animation of the cut iterations (cut-variation mode only).

Elastic mode writes the PDF/ROOT with canvases c0 (no bg, coincidence-time cut), c2 (MC + poly
bg), c3 (side-band poly fit; signal counts printed to the terminal) and c4 (MC + inelastic
MC bg). Its CSV is created but stays empty.

**Feeds:** the `Rk`, `Rkerr` and `RpMC` columns of the production CSVs (and of the variant
CSVs) go into the results table and its systematics ([07](./07_results.md)). **The central
value is fit 3 (pol2 background, the production background shape)**: `R = R3 × RpMC` and
`R_stat = R3err × RpMC`. The background-model systematic (`Inel`) is the spread of
R1–R7 × RpMC.

**Notes:**
- **Data signal cut.** Since 2026-10-21 the data signal cut is
  `cuts_for_signal_simu && extra_data_signal_cut`. All production `qelas` blocks (SBS-4
  30p/50p, SBS-9 70p, SBS-8 70p, SBS-14 70p, SBS-7 85p, SBS-11 100p) have
  `extra_data_signal_cut` with the data-only cuts (GRINCH and/or coincidence time, e.g. SBS-9:
  `cltrindexGRINCH==0&&clsizeGRINCH>2&&abs(coinT_ADC_TOF_c)<=3.5`). `cuts_for_signal_data`
  is not read; it only holds a note saying so.
- The MC kinematics (⟨Q²⟩, ⟨ε⟩) use the plain `weight`, not the efficiency weight.
- Fixed 2026-11-12 (CSVs made before that are affected): `Yn4err` held the background yield
  instead of the n-yield error, and c7 was saved as `_6.png`, overwriting the c6 PNG.
- With `use_custom_fiduCut=1`, the ε histogram is booked with range 0–16 instead of 0–1.5,
  which coarsens ⟨ε⟩.
- Cut-variation styles (`cut_vary_style`): 0 = slices of fixed width; 1 = window of growing
  width around a mean; 2 = slices of equal statistics; 3 = rising threshold; 4 = vary the
  HCAL fiducial safety margins by `fidu_vary_*` percent; 5 = vary each cut in
  `param_to_vary_for_cs5` by ±10% (the cuts being varied must **not** also be in
  `cuts_for_signal_simu`).
- Before 2026-10 R_MC was computed differently. Results made with the old code need the
  R_MC fix described in [07](./07_results.md).

---

### `scripts/fit_W2.cpp`

**What it does:** Fits the BigBite W² distribution of data with MC shapes. The main use is
LH2 elastics, where it checks the MC W² peak position (the `W2_offset` applied to the MC)
and the inelastic background under the peak. The MC signal can be a **weighted mix of
several MC samples**. This is used at SBS-11, where the run had shield-on and shield-off
periods.

**Run** (from `scripts/`):
```
root -l -b -q 'fit_W2.cpp("configs/sbs11-sbs100p/conf_fit_W2.json", 1, "pdout/fit_W2")'
```
- arg 2: `is_elastic`, default `1` (LH2, block `elas`); `0` uses block `qelas`.
- arg 3: output file stem, default `"pdout/test_fit_W2"`.

**Inputs:** `configs/<setting>/conf_fit_W2.json`, block `elas` (or `qelas`):

| Key | Meaning |
|---|---|
| `SBS_config`, `SBS_magnet_percent`, `model`, `pass`, `generator` | file selection |
| `data_file_prefix` | `pdout/<prefix>_elas_ana_data_sbs<conf>_sbs<mag>p_model<model>_pass<pass>.root` |
| `simu_file_prefix` | single MC sample `simulation/siout/<prefix>_elas_ana_<gen>_sbs..._model<model>.root` |
| `simu_mix_prefixes`, `simu_mix_weights` | optional. If present, these override `simu_file_prefix`: the list of MC prefixes and their mix fractions. For SBS-11 100p they are `[0.607, 0.393]` for Al-shield / no-shield, the charge fractions. Each file's `weight` is already normalized per generated event, so the weights are plain fractions |
| `inel_file` | inelastic MC file (relative to `simulation/siout/`) used as a background shape. SBS-11 LH2 uses the LD2 inelastic MC (no LH2 inelastic MC exists) |
| `cuts_for_inclusive_W2` | data cuts for the full W² spectrum that is fitted |
| `cuts_for_signal_data`, `cuts_for_signal_simu` | signal cuts (data template for canvas 2, MC template) |
| `cuts_for_bg_data`, `cuts_for_bg_simu` | bg selections (data bg is filled but not used in a fit; MC bg = inelastic file) |
| `coinT_cut` | coincidence-time cut for canvas 0 |
| `W2_peakOff_for_MC_sig`, `W2_peakOff_for_MC_bg` | extra shift added to MC W² (signal, bg) **on top of** the `W2_offset` already applied in the MC trees |
| `vary_W2Off_ranges` | `[flag, sig_min, sig_max, bg_min, bg_max]`. flag=1: canvas 3 floats both shifts (and ignores the fixed ones) |
| `Order_of_poly_bg`, `h_W2`, `W2_fit_range` | polynomial order, binning, fit range |

**Outputs:** one PDF, `<stem>_sbs<conf>_sbs<mag>p_model<model>.pdf`, with 4 canvases.
**c0** is data (coincidence-time cut) vs MC signal with no bg. **c1** is MC signal + polynomial
bg; the total, background and signal counts are printed. **c2** uses the data signal shape +
polynomial bg. **c3** is MC signal + inelastic-MC bg; with flag=1 its fit parameters
`xOff_s`, `xOff_bg` are the fitted MC shifts. No ROOT file is written (the `.root` name is
used only to derive the PDF name).

**Feeds:** the W² shift that matches data and MC. Enter it into the MC configs
(`W2_offset`, MC stage) or as `W2_peakOff_for_MC_sig` here. The MC Δx templates used by
`fit_dx.cpp` depend on it through the W² cut.

**Notes:**
- Only the SBS-11 100p `elas` block is up to date (pass 3, `inel_file`, mix). The other
  settings' `conf_fit_W2.json` (SBS-4 30p/50p, SBS-9 70p, SBS-14 70p, SBS-7 85p) are
  pass-2-era and have no `inel_file`, so they need updating before use. There is no SBS-8
  `conf_fit_W2.json`. (verify)
- The setting-by-setting W2-offset determination is a by-eye LH2 peak match with strict
  elastic cuts, not this fit. `fit_W2.cpp` is the cross-check.

---

### `scripts/fit_dy.cpp`

**What it does:** Fits the HCAL Δy distribution (no p/n separation in y) of data with the
MC shape, with no bg (c0, with the coincidence-time cut) and with a polynomial bg (c1).
It is a quick look at the Δy shape/width and the bg level, for example when choosing the
Δy cut. It is not part of the production chain.

**Run** (from `scripts/`):
```
root -l 'fit_dy.cpp("configs/sbs4-sbs30p/conf_fit_dy.json", 0)'
```
arg 2 `is_elastic` (default 1), arg 3 an output stem (unused). Run it interactively: the
canvases are only drawn, nothing is saved.

**Inputs:** `configs/<setting>/conf_fit_dy.json` (exists only for SBS-4 30p, pass 2):
`data_file_prefix`, `simu_file_prefix`, `cuts_for_signal_data/_simu`, `cuts_for_bg`,
`coinT_cut`, `dy_peak_offset_for_MC`, `Order_of_poly_bg`, `h_dy`, `dy_fit_range`.

**Outputs:** none (screen only).

**Feeds:** nothing directly. It informs the Δy cut / anti-Δy background choice in
`conf_fit_dx.json`.

**Notes:** the config is pass-2-era (old prefixes, `WCut`/`fiduCut` tree flags), so update
the prefixes before use.
