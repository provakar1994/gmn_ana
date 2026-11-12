# Legacy (obsolete) scripts

Back to the [analysis guide](../ANALYSIS_GUIDE.md).

The scripts below are kept for reference only. They are **not maintained**: they may not
compile against the current library, may read old file names or config keys, and their
results are superseded. Do not use them for results. Paths relative to `scripts/`.

## Obsolete

| Script | What it was |
|---|---|
| `fit_qelas_data.cpp` | early fit of QE data distributions (header has no description); replaced by `fit_dx.cpp` |
| `hcal_dxdy_ld2.cpp` | first (2022) LD2 macro plotting p and n spots in HCAL dx/dy to define QE cuts (`.cfg` configs); replaced by `qelas_ana_data.cpp` |
| `simulation/hcal_dxdy_simu.cpp` | MC version of `hcal_dxdy_ld2.cpp`; replaced by `simulation/qelas_ana_simu.cpp` |
| `npratio_data.cpp` | n/p ratio stability vs cuts, data; superseded by the cut-variation mode of `fit_dx.cpp` |
| `npratio_all.cpp` | n/p ratio stability plots for all settings together; superseded as above |
| `simulation/npratio_simu.cpp` | MC counterpart of `npratio_data.cpp` |
| `yield_per_run.cpp` | yields per run (from `conf_fit_dx.json` `*_prod` blocks, plus `conf_yield_per_run.json`) to check run-by-run stability |
| `simulation/simu_data_fit.cpp` | first data/MC dx fit with the MC shapes (`conf_simu_data_fit.json`); replaced by `fit_dx.cpp` |
| `simulation/old_simu_data_fit.cpp` | older single-parameter-R version of `simu_data_fit.cpp` |
| `pedAlgo_study.cpp` | 2022 comparison of old vs new FADC pedestal-subtraction algorithms |
| `lta/rdf.cpp` | scratch file for learning ROOT RDataFrame |
| `lta/fEvtHdr.cpp` | scratch macro reading the event-header (`fEvtHdr`) branch |
| `simulation/varyR.sh` | loop running `simu_data_fit.cpp` for a range of R values |
| `simulation/chi2_vs_R_plot.py` | parabola fit of chi² vs R from the `varyR.sh` scan |
| `simulation/ML/` (`nu_det_hcal_ML.ipynb`) | machine-learning study of HCAL neutron detection on `qelas_ana_simu` ML output (`gen_ML_data` option) |
| `pass3/` (`python_utils.py`) | plotting helpers for pass-3 cut-scan trend plots |
| `configs/*/hist_lib_practice.ipynb` | practice notebook for the Python `hist` library |
| `hcal_det_effi_pblk.cpp` (+ `configs/sbs8-sbs0p/conf_hcal_det_effi_pblk.json`) | 2022 per-block HCAL proton efficiency from replayed LH2 data (bad-block search); uses `Cut.h` signatures no longer in the library |
| `simulation/hcal_det_effi_simu.cpp` (+ `simulation/configs/sbs4-sbs50p-simu/conf_hcal_det_effi_simu.json`) | early MC HCAL nucleon-detection study with an old normalization; superseded by `simulation/pDE_map_simu.cpp` |
| `simulation/hcal_det_effi_gun.cpp` (+ `simulation/configs/sbs4-sbs0p-simu/conf_hcal_det_effi_ngun.json`, `conf_hcal_det_effi_pgun.json`) | 2022 HCAL n/p detection efficiency vs momentum from g4sbs particle guns (threshold E_peak/4) |

## Folder kept for reference

| Folder | What it is |
|---|---|
| `thesis/` | thesis (chapter 4) figure macros, notebooks and the produced PDFs/PNGs; kept for reference, not part of the analysis chain |

## Probably obsolete - confirm

Found in `git ls-files`; no in-use script reads them, or they belong to the obsolete scripts
above. Confirm before moving them here.

| File(s) | Why it looks obsolete |
|---|---|
| `configs/*/conf_fit_qelas_data.json` (6 settings) | config of `fit_qelas_data.cpp` |
| `configs/*/conf_npratio_data.json`, `configs/conf_npratio_all_pass{2,3}.json`, `simulation/configs/sbs9-sbs70p-simu/conf_npratio_simu.json` | configs of the npratio scripts |
| `configs/*/conf_yield_per_run.json` | config of `yield_per_run.cpp` |
| `simulation/configs/*/conf_simu_data_fit.json`, `configs/*/sbs*_hcaldxdy.cfg`, `simulation/configs/*/sbs*_hcaldxdy_simu.cfg` | configs of `simu_data_fit.cpp` / `hcal_dxdy_*.cpp` |
| `configs/*/conf_hcal_det_effi_data.json` (sbs4-0p, sbs4-50p, sbs9-70p, sbs7-85p) | config for `hcal_det_effi_data.cpp`, which no longer exists in the repo |
| `configs/sbs4-sbs50p/conf_physics_ana.json` | no script reads it |
| `configs/sbs4-sbs50p/systematic_studies.ipynb`, `varyR_vs_{W2,coinT,eHCAL,ePS}.csv`, `databg_sbs4_sbs50p_sys_study.csv`, `poly3_sbs4_sbs50p_sys_study.csv` | 2023 pass-1 systematic studies (R vs cut), superseded by `fit_dx.cpp` cut variations |
| `configs/sbs11-sbs100p/temp_conf_qelas_ana_data.json` | temporary copy of `conf_qelas_ana_data.json` |
| `configs/sbs4363-sbs100p/conf_qelas_ana_data.json` (and config `4363` in `src/ExpConstants.cpp`) | GEn setting (E_beam 6.373 GeV) added at the very beginning only to verify the code; never used again, kept on purpose (no harm) |
| `epics/bcmCheckCharge.cpp`, `epics/conf_bcmCheckCharge.json` | BCM check with a pass-0 path and a run list not in git |
| `epics/sbs11_ld2_ebeam_ana.ipynb` | exploratory study, superseded by `epics/get_avgHALLA_p_prun.ipynb` |
| `dflay/` except `JSONManager` (and the BCM classes used by `bcmCheckCharge.cpp`) | BCM calibration code copied from SBS-ANA in 2022; not used by the analysis |
| `DB/good_runList_GMn_nTPE_*_pass_1.csv` | old 8-column format, no longer readable by `CodaRun`; keep as history only |
| `DB/LH2/.....` | stray file (SBS-11 LH2 beam energies, old column names) |
| `DB/**/.ipynb_checkpoints/*`, `scripts/thesis/**/.ipynb_checkpoints/*` (23 files) | Jupyter autosaves committed by mistake (`.gitignore` now excludes them) |
| `scripts/thesis/**/root_hist` (10 files) | ROOT command-history files committed by mistake |
