# Change log

## 2026-10-28 11:54
- **Files:** `include/gmnana_LinkDef.h`, `src/Utilities.cpp`
- **Change:** Fix the build with ROOT 6.40. The LinkDef guard now also checks `__CLING__` (rootcling no longer defines `__CINT__`), and it gains one positive selection rule (`namespace constant`), since ROOT 6.40 rejects a LinkDef with only `link off` rules. In `Utilities.cpp`, pass `lfname.c_str()` instead of a `std::string` to `Form("%s")`, which the new compiler rejects.

## 2026-10-28 11:54
- **Files:** `src/KinematicVar.cpp`
- **Change:** `CalcGMnError` now uses the exact error propagation, dGMn = term1·σ_R,p·dR / (2τ_n|GMn|), instead of the approximation |GMn|·dR/(2R), which ignored the GEn term, and no longer prints debug lines on every call. `ExtractGMnWithError` called the 7-argument `ExtractGMn` (the `ebeam, eeprime, etheta` overload) instead of `CalcGMnError`, so it returned a wrong error; fixed.

## 2026-10-28 11:54
- **Files:** `scripts/extractGMn.cpp`
- **Change:** Rewrite the macro around one helper (`calcGMn`) that gives GMn/(μ_n G_D) with propagated stat errors and per-source exp. systematic errors for each EMFF model. It takes named systematic sources (all uncorrelated), checks that the number of names matches the number of errors, and reports the theory systematic as two sources: σ_p (Kelly/Arrington spread + Ye GEp/GMp fit errors) and GEn (Kelly/Galster spread + Ye GEn fit error). The CSV headers are built from the source names. It also fixes printing only `sys[0]`, and removes the hard-coded CSV header.

## 2026-10-28 12:10
- **Files:** `DB/FF_fits/Ye2018/*`
- **Change:** Track the Ye 2018 form-factor material with git: the paper and its supplementary material, the official lookup tables with TPE (`ye2018_wTPE_*`), and the proton lookup table without external TPE (`ye2018_NoTPE_proton_lookup.csv`, from private communication). The library will read these tables from here, not from a CSV in the working directory.

## 2026-10-28 12:40
- **Files:** `CMakeLists.txt`, `include/EMFFFits.h`, `src/EMFFFits.cpp`, `include/LookUpTableReader.h`
- **Change:** Fix the Ye no-TPE proton FFs (`Ye2017::GetFFNoTPE`, `GetFFwErrNoTPE`). They read "GEp" from the `muGE/GM` column (the FF ratio, ~50x too large at Q2 = 4.5), took the nearest table row instead of interpolating (biasing GMp by up to ~1% and GEp by up to ~11%), and re-read a `data.csv` from the current working directory on every call. Now the table is read once from `DB/FF_fits/Ye2018/ye2018_NoTPE_proton_lookup.csv`, columns are found by name (`GEp`, `dGEp`, `GMp`, `dGMp`), and values are interpolated linearly in log(Q2); Q2 outside the table gives an error. CMake now sets `GMNANA_DB_DIR` (overridable by the env. variable of the same name), read via the new `EMFFFits::GetDBPath()`. `LookUpTableReader` gains `GetLogInterpValueByColumnName()` and strips a UTF-8 BOM from the header.

## 2026-10-28 13:05
- **Files:** `src/EMFFFits.cpp`, `include/EMFFFits.h`
- **Change:** Fix FF fits found while checking them against the original papers. (1) Galster GEn was missing the factor tau in the numerator (GEn/GD was 1.913 at Q2 = 0 instead of 0, and 0.235 instead of 0.300 at Q2 = 4.5); now GEn/GD = -mu_n*tau/(1+5.6*tau), cited from Perdrisat et al., PPNP 59, 694 (2007), Eq. (16). (2) Kelly used the neutron mass in tau for GEn and GMn; the paper (Eq. 1) and SIMC use the proton mass for all four FFs, so now Kelly does too (Kelly GMn changes by ~-0.5%). (3) Seamus GEn now also uses the proton mass in tau, as in g4sbs and SIMC; its comments now say that it is unpublished and not the fit of Riordan et al., PRL 105, 262302 (2010), and should only be used to reproduce the MC. (4) Comments: correct the Arrington2007 citation (was Christy's), note Arrington's GE validity (Q2 <= 6 GeV2), rename the misnamed `GMp_ov_mun` in the Christy fit, fix two warning texts, and replace the Ye error-unit "check needed" note with the verified explanation.

## 2026-10-28 13:30
- **Files:** `include/EMFFFits.h`, `src/EMFFFits.cpp`, `src/KinematicVar.cpp`, `scripts/extractGMn.cpp`
- **Change:** Make the MC Born CS ratio R_MC match the SIMC generator (bigbite branch, `physics_proton.f`) exactly. New class `SIMCbigbite` gives the EMFFs as SIMC computes them (Kelly GEp/GMp/GMn and Seamus GEn with tau from the proton mass, and SIMC's rounded GMp(0) = 2.79, GMn(0) = -1.913). `kine::sigmaBorn_ratio_MC` now uses it, and `extractGMn.cpp` calls `kine::sigmaBorn_ratio_MC` instead of computing R_MC by hand. R_MC was 0.8-0.9% too high before (mostly from the neutron-mass tau in Kelly GMn); it now agrees with a line-by-line copy of the SIMC code to ~1e-16. Effect on the extractGMn test point (Q2 = 4.5): GMn/(mu_n GD) 0.9347 -> 0.9307; GEn model systematic 0.0167 -> 0.0141 (Galster fix).

## 2026-10-28 13:50
- **Files:** `scripts/fit_dx.cpp`
- **Change:** Compute `Rp_MC` (written to the fit output and used by `scripts/lta/get_results.ipynb` to form R = R_fit × RpMC) with `kine::sigmaBorn_ratio_MC`, i.e. exactly the SIMC generator EMFFs, instead of by hand from `Kelly2004` and `Seamus20XX`. Before the library fixes this R_MC was 0.8-0.9% too high; R values produced with the old `fit_dx.cpp` must be regenerated (or divided by the old/new R_MC ratio).

## 2026-10-28 14:15
- **Files:** `scripts/qelas_ana_data.cpp`, `scripts/simulation/qelas_ana_simu.cpp`
- **Change:** Fill the per-event `sigBorn_ratio` branch with `kine::sigmaBorn_ratio_MC`, i.e. exactly the SIMC generator EMFFs, instead of Kelly2004 + Seamus20XX. After the Kelly/Seamus tau fix that combination is still ~0.2% low w.r.t. SIMC, because Kelly2004 uses the exact magnetic moments while SIMC uses 2.79 and -1.913 (it was ~0.8% high before the tau fix). The `sigRed_p`/`sigRed_n` branches are unchanged. Output trees made before this change keep the old values.

## 2026-10-28 15:30
- **Files:** `DB/FF_fits/Ye2018/README.md`, `DB/FF_fits/Ye2018/ye2018_noExtraTPE_proton_lookup.csv` (renamed from `ye2018_NoTPE_proton_lookup.csv`), `DB/FF_fits/Ye2018/Round6_all_sumrules_leastsq_Q21000_z12_gb5_t0fix7_noExtraTPE.dat` (new), `include/EMFFFits.h`, `src/EMFFFits.cpp`
- **Change:** Correct the label of the Ye "no-TPE" proton fit. The `noExtraTPE` fit only drops the phenomenological extra TPE term (PLB 777, Eq. 4) and is still SIFF-Blunden TPE corrected: sigma_R(noExtraTPE)/sigma_R(wTPE) = 0.990-0.997 for Q2 = 3-13.5 GeV2, nearly independent of eps, and mu_p GEp/GMp = 0.66 at Q2 = 3. Renamed the CSV to say what it is, documented this in the README, and made `Ye2017::GetFFNoTPE` print a one-time placeholder warning until a truly TPE-uncorrected Ye fit (requested from J. Arrington) is available. Also track the source `.dat` file (identical to the copy in `proton_lookup_table.zip`), which holds the z-expansion coefficients and covariance matrix; README notes that neglecting the GEp-GMp correlation overestimates the sigma_R,p error by 0-8% at the GMn kinematics. No numerical change.

## 2026-10-28 16:30
- **Files:** `CMakeLists.txt`, `tests/CMakeLists.txt`, `tests/TestUtils.h`, `tests/test_lookup_table.cpp`, `tests/test_simc_rmc.cpp`, `tests/data/log_grid.csv`, `README.md`
- **Change:** Add unit tests (CTest, no external framework), built by default (option `GMNANA_BUILD_TESTS`) and run with `ctest --output-on-failure`. `TestUtils.h` gives `CHECK_REL`/`CHECK_ABS`/`CHECK_TRUE`. `test_lookup_table` checks `LookUpTableReader` (BOM stripping, exact values at nodes, log-linear interpolation, -99999 on out-of-range/unknown column) on a small fixture. `test_simc_rmc` checks that `SIMCbigbite` and `kine::sigmaBorn_ratio_MC` reproduce a line-by-line copy of the SIMC bigbite `physics_proton.f` FFs to 1e-12 (0.5-15 GeV2 and the GMn kinematics). Verified that the SIMC test fails when a coefficient is altered.
