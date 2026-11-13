# 01 Setup: build, run, test

Back to the [analysis guide](../ANALYSIS_GUIDE.md).

## Requirements
- ROOT 6 (found through `$ROOTSYS`; `include/gmnana_LinkDef.h` is written to work with ROOT >= 6.40
  as well), CMake >= 3.9, a C++17-capable compiler (ROOT's standard is used; the project does
  not set one).
- On the JLab ifarm the batch scripts use the Hall A `analyzer` module (`module load
  analyzer/1.7.12`), which brings its own ROOT.
- Python 3 with pandas / numpy / matplotlib (and `hist`, `uproot`, `seaborn` for some
  notebooks) for the `.ipynb` / `.py` helpers.

## Build the library

The library `gmnana` (`libgmnana.so` / `.dylib`) is built from `src/*.cpp` + `include/*.h` by
`CMakeLists.txt`. Use a build and an install directory outside the source tree, e.g. next to
the repository:

```
mkdir build_gmn_ana install_gmn_ana
cd build_gmn_ana
cmake -DCMAKE_INSTALL_PREFIX=$PWD/../install_gmn_ana ../gmn_ana
make -j4
make install
ctest --output-on-failure        # unit tests, see below
```

What the build does:
- Builds `libgmnana` and a ROOT dictionary (`gmnanaDict`, module `gmnana`) in the build directory.
- Bakes the path of the source `DB/` directory into the library (`GMNANA_DB_DIR`, used e.g. for the
  FF lookup tables). Override at run time with the environment variable `GMNANA_DB_DIR`.
- Generates `<build>/runenv/rootlogon.C` and `<build>/runenv/gmn_ana.rootrc` from
  `rootlogon_install.C` and `gmn_ana.rootrc`; `make install` copies them to `<install>/runenv/`
  and the headers to `<install>/include/`.
- Builds the unit tests unless `-DGMNANA_BUILD_TESTS=OFF`.

`rootlogon.C` (from `rootlogon_install.C`):
- first runs a `./rootlogon.C` or `~/rootlogon.C` if present;
- adds `-I<install>` and `-I<build>/include` to the include path;
- loads `libgmnana` **from the build directory** (`${PROJECT_BINARY_DIR}`), so the build directory
  must be kept. After a library change, `make` (+ `make install` if headers changed) is enough.

`gmn_ana.rootrc` sets `Rint.Logon` to the installed `rootlogon.C`, the dynamic path to
`<install>/lib64`, and `ACLiC.IncludePaths` to `<install>/include`. It can be used as a
`.rootrc` so that plain `root` loads the library automatically (verify); the scripts in this repo
instead pass `rootlogon.C` explicitly on the command line.

Library headers listed in `CMakeLists.txt` (`allheaders`) are installed; `LookUpTableReader.h`
is header-only and not in that list (verify it is still found by scripts that need it).

## Running scripts

Analysis scripts are ROOT macros. They include the library headers by relative path
(`#include "../include/gmn_ana.h"`, `"../../include/gmn_ana.h"` from `simulation/`) and the JSON
reader by source (`#include "../dflay/src/JSONManager.cxx"`), and need `libgmnana` loaded by
`rootlogon.C`. General form:

```
cd scripts                               # data / data-MC scripts
root -l -b -q <install>/runenv/rootlogon.C 'fit_W2.cpp("configs/sbs11-sbs100p/conf_fit_W2.json",1)'

cd scripts/simulation                    # MC scripts
analyzer -b -q <install>/runenv/rootlogon.C \
  'qelas_ana_simu.cpp("configs/sbs14-sbs70p-simu/conf_qelas_ana_simu.json",2,"simc","deeN")'
```

- `-l` no splash, `-b` batch (no graphics windows), `-q` quit when done.
- The first macro argument is always the config file path, **relative to the working
  directory**; further arguments select the config block (pass/model, generator/process,
  `elas`/`qelas`); see each card in the stage files.
- The working directory matters: output and input paths inside the macros (`pdout/`, `siout/`,
  `simulation/siout/`, `../DB`, `epics/epout/`) are relative. Data and data/MC scripts run from
  `scripts/`; MC scripts run from `scripts/simulation/`; `gmn_extraction/` scripts run from their
  own directory; `epics/` scripts run from `scripts/epics/`.
- `w2_offset_check.cpp` documents a bare `root -l -b -q 'w2_offset_check.cpp(...)'`, which
  works only if the library is auto-loaded through `gmn_ana.rootrc` / a `rootlogon.C` in the
  path (verify).
- On the ifarm, use `analyzer` instead of `root` (same syntax). The batch wrappers
  `run-ana-data.sh` and `simulation/run-ana-simu.sh` need `GMN_ANA` (repository path) and
  `ROOTLOGON_ENV` (installed `rootlogon.C`). Set both in your shell before submitting;
  `sbatch` passes the environment to the job:
  ```
  export GMN_ANA=<path to gmn_ana> ROOTLOGON_ENV=<install>/runenv/rootlogon.C
  ```
  If they are not set, the original author's ifarm paths are used. The wrapper stops with an
  error if `$GMN_ANA/scripts` or `$ROOTLOGON_ENV` is not found. The wrappers pass the config
  as `configs/sbs<K>-sbs<M>p[-simu]/conf_...json`.

## Unit tests

`tests/test_*.cpp`, one executable per topic, linked to `libgmnana`; small fixtures in
`tests/data/`. Run all from the build directory with `ctest --output-on-failure` (one test:
`ctest -R test_ff_fits --output-on-failure`). Add or update a test with every library change
that affects physics results.

| Test | Checks |
|---|---|
| `test_lookup_table` | `LookUpTableReader`: BOM stripping, CRLF, node values, log-Q² interpolation, edges, error value |
| `test_simc_rmc` | `SIMCbigbite` FFs reproduce SIMC `physics_proton.f`; `sigmaBorn_ratio_MC` uses exactly them |
| `test_gmn_extraction` | `kine::` GMn extraction round trip and analytic error derivatives vs numerical ones |
| `test_ff_fits` | Ye2017 (with/without TPE), Kelly, Galster against official tables/papers |
| `test_gmn_module` | `gmnx::` model rows, TPE groups, R = R_MC*R_fit, error totals |
| `test_eloss_tgt` | `util_pd::GetElossInTgt` (data run-number overload vs explicit-shield core, SBS-11 shield periods) |
| `test_load_simu_tree` | `util_pd::LoadSimuROOTTree` drops missing/empty/corrupt MC files so the normalization stays right |

## Library modules (`include/` + `src/`)

`include/gmn_ana.h` includes all modules except `GMnExtraction.h` (include that one directly).

| Header | Namespace / type | Provides |
|---|---|---|
| `Constants.h` | `constant` | physical constants (masses, magnetic moments, charge, ...) |
| `ExpConstants.h` / `.cpp` | `expconst`, `class SBSconfig` | detector geometry (BBCAL, HCAL blocks), per-config E_beam, BB/SBS/HCAL angles and distances, target and shield materials (`ShieldType`), SBS max current; `SBSconfig` bundles them for one setting |
| `CodaRun.h` | `struct CodaRun` | one data run as read from the good run list (run, config, fields, E_beam, charge, live time) |
| `SimuJob.h` | `struct SimuJob` | one MC job as read from the job summary CSV (files, generator, process, tries, volume, luminosity, charge, ...) |
| `KinematicVar.h` / `.cpp` | `kine` | elastic/QE kinematics (Q², W², epsilon, expected nucleon direction and HCAL position), Mott/Born/reduced cross sections, n/p Born ratio (incl. the SIMC-FF version R_MC), GMn extraction and its error propagation |
| `Cut.h` / `.cpp` | `cut` | spot cuts, HCAL active-area and safety-margin (fiducial) cuts for data and MC |
| `Fit.h` / `.cpp` | `fit` | fit helpers: parameter names, values/errors, standard fit functions |
| `FitFns.h` | `class FitFn` | TF1 functors: histogram templates (MC p/n signals, backgrounds) + polynomials, used by the dx/W2 fits |
| `SetROOTVar.h` / `.cpp` | `setrootvar` | `setbranch` helpers to switch on and address tree branches |
| `Utilities.h` / `.cpp` | `util_pd` | read run lists (`ReadRunList`), chain data files (`LoadROOTTree`), read MC job summaries and chain MC files (`ReadSimuJobSummary`, `LoadSimuROOTTree`), total charge/tries, luminosity, target energy loss (`GetElossInTgt`), histogram/plot styling helpers |
| `EMFFFits.h` / `.cpp` | `class EMFFFits` + fits | nucleon FF parametrizations: Galster1971, Kelly2004, Seamus20XX (Riordan GEn), SIMCbigbite (the MC generator FFs), Ye2017 (with/without TPE, via `DB/FF_fits/`), Christy2022 (no TPE only), Arrington2007 |
| `LookUpTableReader.h` | `class LookUpTableReader` | CSV lookup tables with closest-key and log interpolation (FF tables) |
| `GMnExtraction.h` / `.cpp` | `gmnx` | GMn for a list of FF models from R, with stat, p2p, scale, pCS and GEn errors; see `scripts/gmn_extraction/README.md` |
| `PlotCustomizer.h` / `.cpp` | `class PlotCustomizer` | thesis/publication plot styling |
| `gmnana_LinkDef.h` | - | ROOT dictionary selection |

## `dflay/`

Copied on 2022-08-22 from D. Flay's `SBS-ANA` repository (see `dflay/README.md`); not part of
the `gmnana` library, included by source where needed.
- `include/JSONManager.h`, `src/JSONManager.cxx` (nlohmann `json.hpp`): the config reader used
  by almost every script (`GetValueFromKey`, `GetValueFromSubKey`, `GetVectorFromSubKey`, ...).
  It reads only the first JSON object in a file, so notes may follow the closing `}`.
- BCM/beam tools (`BCMManager`, `bcm*.cxx`, `unser*.cxx`, `beamPlot.cxx`, ...): BCM calibration
  code; only `scripts/epics/bcmCheckCharge.cpp` uses parts of it (see
  [02](02_run_info.md) and [legacy](legacy.md)).

## `DB/`

| Path | Content |
|---|---|
| `good_runList_GMn_nTPE_<target>_pass_<P>.csv` | good-run lists per target (LD2, LH2; Dummy, Optics for pass 2) and replay pass (1, 2, 3). Columns: `sbsconfig, runnum, target, sbsmagcurr(A), bbmagcurr(A), ebeam_mean(GeV), ebeam_std(GeV), totcharge_bob(C), totcharge_pd(C), daqlvtm_bbhi, daqlvtm_ts1` (pass 1 has fewer columns). Read by `util_pd::ReadRunList`; pass 0 uses the pass-1 file |
| `LD2/`, `LH2/`, `Dummy/`, `Optics/` | per-config inputs to the run lists: `avg_ebeam_prun_sbs<K>_<target>.csv` (EPICS beam energy, before the 1.003 correction), `beamCh_n_daqLvTm_prun_SBS<K>_<target>.csv` (scaler charge and live time), `all_beamCh_n_daqLvTm_prun_<target>.csv` (all configs + online charge) |
| `*.ipynb` | `get_run_list.ipynb` (print run numbers as bash arrays), `add_scaler_files.ipynb`, `join_n_merge_rLists.ipynb` (merge the above into the run lists); see [02](02_run_info.md) |
| `FF_fits/Ye2018/` | Ye et al. (PLB 777 (2018) 8) FF tables (with TPE, and the "noExtraTPE" proton fit used as a no-TPE placeholder) + paper; read `README.md` there before using them |

## Output locations

- `scripts/pdout/`: data analysis outputs (and data/MC fits, e.g. `pdout/fits/`, `pdout/pDE/`).
- `scripts/simulation/siout/`: MC outputs (and e.g. `siout/w2_offset_check/`, `siout/pDE/`).
- Both are normally **symlinks to a data disk** (on the ifarm: a volatile/work area; locally: an
  external disk with the copied outputs). They are not in git. Create them before the first run,
  e.g. `ln -s <data disk>/pdout scripts/pdout`.
- Rule: every study writes into `pdout/` or `siout/` (sub-directories allowed) and puts its plots
  into **one multi-page PDF** (ROOT `canvas->SaveAs("file.pdf[")`, pages, `"file.pdf]"`), with the
  ROOT file / CSV of the same base name next to it.
- Small run-info outputs go to `scripts/epics/epout/` (see [02](02_run_info.md)).
- Slurm logs of the batch wrappers go to `/farm_out/$USER/gmn_ana/slurm_out/` (the JLab farm only
  accepts log paths under `/farm_out/`; override with `GMN_SLURM_LOGDIR`). The `submit-*.sh` scripts
  create the directory.
