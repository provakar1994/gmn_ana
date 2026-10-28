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
