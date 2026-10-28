# Change log

## 2026-10-28 11:54
- **Files:** `include/gmnana_LinkDef.h`, `src/Utilities.cpp`
- **Change:** Fix the build with ROOT 6.40. The LinkDef guard now also checks `__CLING__` (rootcling no longer defines `__CINT__`), and it gains one positive selection rule (`namespace constant`), since ROOT 6.40 rejects a LinkDef with only `link off` rules. In `Utilities.cpp`, pass `lfname.c_str()` instead of a `std::string` to `Form("%s")`, which the new compiler rejects.

## 2026-10-28 11:54
- **Files:** `src/KinematicVar.cpp`
- **Change:** `CalcGMnError` now uses the exact error propagation, dGMn = term1·σ_R,p·dR / (2τ_n|GMn|), instead of the approximation |GMn|·dR/(2R), which ignored the GEn term, and no longer prints debug lines on every call. `ExtractGMnWithError` called the 7-argument `ExtractGMn` (the `ebeam, eeprime, etheta` overload) instead of `CalcGMnError`, so it returned a wrong error; fixed.
