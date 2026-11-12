# 02 Run information: beam energy, charge, live time, run lists

Back to the [analysis guide](../ANALYSIS_GUIDE.md).

The data analysis takes everything it knows about a run from the **good run lists**
`DB/good_runList_GMn_nTPE_<target>_pass_<P>.csv` (read by `util_pd::ReadRunList` into
`CodaRun` objects). The scripts in `scripts/epics/` and the notebooks in `DB/` produce the
per-run beam energy, beam charge and DAQ live time that were merged into those lists. This is
done once per replay pass; the current lists (pass 3) are already in `DB/`, so this stage only
needs rerunning if a run list has to be rebuilt.

## Run-list columns and how they are used

| Column | Source | Used as |
|---|---|---|
| `sbsconfig`, `runnum`, `target`, `sbsmagcurr(A)`, `bbmagcurr(A)` | run database / logbook (verify origin) | run selection: config, target, SBS field (% = current / 21 A) |
| `ebeam_mean(GeV)`, `ebeam_std(GeV)` | EPICS `HALLA_p` per run x 1.003 (see below) | per-run beam energy in the kinematics |
| `totcharge_bob(C)` | Hall A online charge tool (`get_charge_from_Bobs_script.sh`) | run charge, **if > 0** |
| `totcharge_pd(C)` | `get_scalerdata_prun.cpp` (BCM `dnew`) | run charge if the online value is <= 0 (e.g. all SBS-4 runs, where the online tool fails) |
| `daqlvtm_bbhi`, `daqlvtm_ts1` | `get_scalerdata_prun.cpp` | DAQ live time; `daqlvtm_bbhi` is the one used |

`CodaRun::SetDataRunSheet` (in `include/CodaRun.h`) requires exactly these 11 columns, so the
pass-1 lists (8 columns) can no longer be read by the current code.

## Data flow

```
 replayed ROOT files (E tree)  --get_avgHALLA_p_prun.cpp-->  epics/epout/get_HALLA_p_prun_SBS<K>_<target>.csv
                               --get_avgHALLA_p_prun.ipynb-> DB/<target>/avg_ebeam_prun_sbs<K>_<target>.csv
                               --DB/join_n_merge_rLists.ipynb (x1.003)--> ebeam columns of the run list

 replayed ROOT files (TSsbs)   --get_scalerdata_prun.cpp (via run_eORs_ana.sh)-->
        epics/epout/scalerdata_prun_SBS<K>_<target>.root   (per-event scaler tree, used by 03)
        epics/epout/beamCh_n_daqLvTm_prun_SBS<K>_<target>.csv  -> copied to DB/<target>/
 online charge tool            --get_charge_from_Bobs_script.sh--> run,charge list
        --DB/add_scaler_files.ipynb--> DB/<target>/all_beamCh_n_daqLvTm_prun_<target>.csv
        --DB/join_n_merge_rLists.ipynb (scaler part)--> charge / live-time columns of the run list
```

The pass-3 lists are the pass-2 lists plus one LD2 run: 11597 (SBS-4, 50p), added by hand
later. The LH2 lists are identical.

> **Note - run 11597:** its four charge/live-time columns are `-1` because they have not
> been computed yet. Anything that uses run charge or live time for this run (e.g.
> `weight_norm`, the charge-normalized yield) gets a meaningless value. This is accepted:
> the analysis does not use the charge-normalized yield (the dx fits use shapes and fitted
> scale factors). Compute the charge/live time first if you ever need absolute yields.

## Scripts

### `scripts/epics/get_avgHALLA_p_prun.cpp`
**What it does:** Loops over the EPICS tree `E` of every good run of one SBS config and target,
and writes one line per EPICS event: run number, global event number, `HALLA_p` (beam energy
EPICS PV, MeV). This is the raw input for the per-run beam energy.
**Run:** from `scripts/epics/`:
`analyzer -b -q <install>/runenv/rootlogon.C 'get_avgHALLA_p_prun.cpp("conf_get_avgHALLA_p_prun.json")'`
(one config/target per run; edit the config between runs).
**Inputs:** `conf_get_avgHALLA_p_prun.json` (flat, no blocks): `target`, `SBS_config`,
`replay_pass` (run list pass; 0/1 -> pass-1 list), `Nruns_to_ana` (-1 = all), `runsheet_dir`
(`../../DB`), `rootfile_dir` (replayed files, matched as `<rootfile_dir>/*<run>*`).
**Outputs:** `scripts/epics/epout/get_HALLA_p_prun_SBS<K>_<target>.csv` (12 files, all configs
x LH2/LD2, are in git).
**Feeds:** `get_avgHALLA_p_prun.ipynb`.
**Notes:** It calls `util_pd::ReadRunList` with 6 arguments; the current library has only 7-
and 8-argument versions, so the macro needs a one-line fix (add a verbosity argument) before it
compiles again (verify). With the current run-list format it also needs pass >= 2.

### `scripts/epics/get_avgHALLA_p_prun.ipynb`
**What it does:** Reads the `HALLA_p` dump, drops entries with `HALLA_p = 0` and duplicated
EPICS events (same run and event number), and computes mean and standard deviation of the beam
energy per run (GeV). Shows per-run plots for checking stability.
**Run:** Jupyter, from `scripts/epics/`, after editing `conf_get_avgHALLA_p_prun.json` to the
config/target wanted.
**Inputs:** `conf_get_avgHALLA_p_prun.json`; `epout/get_HALLA_p_prun_SBS<K>_<target>.csv`.
**Outputs:** `DB/<target>/avg_ebeam_prun_sbs<K>_<target>.csv` (`runnum, ebeam_mean(GeV),
ebeam_std(GeV)`). The writing cell is commented out on purpose ("execute only when super
confident"); uncomment to write.
**Feeds:** `DB/join_n_merge_rLists.ipynb`.
**Notes:** These files hold the raw EPICS value. The 1.003 correction (HALLA_p reads 0.3% low
compared with a 2016 Hall A energy measurement, logbook entry 3448841) is applied only when
merging into the run list.

### `scripts/epics/sbs11_ld2_ebeam_ana.ipynb`
**What it does:** First exploratory look at `HALLA_p` for SBS-11 LD2 (zero values, duplicate
EPICS events, outlier runs). The method it developed is what `get_avgHALLA_p_prun.ipynb` does.
**Run:** Jupyter, from `scripts/epics/`.
**Inputs:** `epout/get_HALLA_p_prun_SBS11_LD2.csv`.
**Outputs:** none (plots in the notebook).
**Feeds:** nothing; reference only (candidate for [legacy](legacy.md), confirm).

### `scripts/epics/get_scalerdata_prun.cpp`
**What it does:** Loops over the scaler tree `TSsbs` of every good run, integrates the BCM
`dnew` counts into charge (gain 3317.99 Hz/uA from the BCM calibration) and computes the DAQ
live time as accepted BBCal-singles triggers (read from the replay logs) / BBCal-Hi (or TS1)
scaler counts. Writes a per-event scaler tree and a per-run summary.
**Run:** from `scripts/epics/`, normally through `run_eORs_ana.sh` (below), or directly:
`analyzer -l -q 'get_scalerdata_prun.cpp("LD2",4,-1,-1,"<data dir>",0)'`
Arguments: target (`LH2`/`LD2`/...), SBS config, SBS field in % (-1 = all fields), number of runs
(-1 = all), data directory (must contain `rootfiles/` and `logs/`), verbosity (>0 = debug and
`test_` prefix on the outputs).
**Inputs:** run list `../../DB/good_runList_GMn_nTPE_<target>_pass_2.csv` (pass **hard-coded to
2** in the macro); replayed ROOT files and replay logs.
**Outputs:** in `scripts/epics/epout/`:
`scalerdata_prun_SBS<K>[_<M>p]_<target>.root` (tree `Tout`: run, segment, event indices,
`dnewcnt`, `dnewcurr`, cumulative `dnewcharge`, trigger scalers) and
`beamCh_n_daqLvTm_prun_SBS<K>[_<M>p]_<target>.csv` (`runnum, totcharge_pd(C), daqlvtm_bbhi,
daqlvtm_ts1`). Neither is in git.
**Feeds:** the CSV (copied to `DB/<target>/`) -> `DB/add_scaler_files.ipynb` -> run list. The
ROOT file is read by `elas_ana_data.cpp` / `qelas_ana_data.cpp` when their config has
`get_scaler_info = 1` (path `epics/epout/scalerdata_prun_SBS<K>_<target>.root` relative to
`scripts/`), to attach scaler information to each event.
**Notes:** Include path is `../../include/gmn_ana.h`; the library must be loaded (the call in
`run_eORs_ana.sh` does not pass `rootlogon.C`, so it relies on an auto-loaded logon (verify)).

### `scripts/epics/run_eORs_ana.sh`
**What it does:** Wrapper that runs `get_scalerdata_prun.cpp` for one SBS config or for all six
(4, 7, 11, 14, 8, 9), with the right replayed-data directory for each.
**Run:** from `scripts/epics/`: `./run_eORs_ana.sh s LD2 -1 -1 -1 0`
Arguments: mode (`s` = scaler; `e` = EPICS is listed but not implemented), target, config (-1 =
all), SBS field (-1 = no restriction), number of runs, verbosity.
**Inputs:** data directories hard-coded in the script (pass-2 locations on the JLab work disk;
edit for another pass/user).
**Outputs:** as `get_scalerdata_prun.cpp`.
**Feeds:** as `get_scalerdata_prun.cpp`.

### `scripts/epics/get_charge_from_Bobs_script.sh`
**What it does:** Calls the Hall A online charge tool (`getChargeGMn`, on the Hall A online
account) for every good LH2, LD2, Optics and Dummy run and prints `run,charge` (-99 on failure
or timeout). This is the preferred run charge (`totcharge_bob`).
**Run:** on a Hall A online machine, `bash get_charge_from_Bobs_script.sh > <file>.csv`.
**Inputs:** run numbers hard-coded in the script (from the pass-2 lists, via
`DB/get_run_list.ipynb`); path of the charge tool hard-coded.
**Outputs:** stdout; saved as `get_charge_from_Bobs_script_gmn_all_tg.csv` (name expected by
`DB/add_scaler_files.ipynb`; not in git).
**Feeds:** `DB/add_scaler_files.ipynb`.
**Notes:** Fails for SBS-4 (clock-count issue); for those runs the run list falls back to
`totcharge_pd`.

### `scripts/epics/bcmCheckCharge.cpp`
**What it does:** Uses the `dflay/` BCM tools to plot the beam current vs event for the first
run of a list and print, per run, the `dnew` charge, run time and DAQ live time. A cross-check
of `get_scalerdata_prun.cpp`.
**Run:** from `scripts/epics/`: `root -l 'bcmCheckCharge.cpp("conf_bcmCheckCharge.json")'`
(verify; it does not need `libgmnana`).
**Inputs:** `conf_bcmCheckCharge.json`: `ROOTfile-path` (replayed files; currently a pass-0
path), `run-path` (`epout/runlist.csv`, not in git).
**Outputs:** canvas + printout only.
**Feeds:** nothing directly (check). Probably obsolete - confirm.

### `scripts/epics/grab_epicsPV_from_evio.py` and `scripts/epics/jcache_0th_seg.sh`
**What it does:** Reads chosen EPICS PVs (default: SBS and BB magnet currents `MSUPERBIGBITEM`,
`MBIGBITEM`) directly from the first segment of the raw EVIO file of every good run, to check
the magnet settings in the run list. `jcache_0th_seg.sh <run>` brings segment 0 of a run from
tape to the cache disk.
**Run:** on the ifarm, from `scripts/epics/`: edit the user inputs at the top (`process` 0 =
jcache, 1 = grep; `target`; `replaypass`; `epicsPVs`), then `python3 grab_epicsPV_from_evio.py`
(run once with `process = 0`, wait for the files, then with `process = 1`).
**Inputs:** `../../DB/good_runList_GMn_nTPE_<target>_pass_<P>.csv`; raw EVIO files on tape/cache.
**Outputs:** `epout/mag_PVs_per_run_<target>.csv` (`runnum,<PV>,...`; -9999 if not found),
`epout/grab_epics_PV_from_evio_<target>.log`.
**Feeds:** check of the `sbsmagcurr(A)`/`bbmagcurr(A)` run-list columns (verify).

### `DB/add_scaler_files.ipynb`
**What it does:** Concatenates the per-config `DB/<target>/beamCh_n_daqLvTm_prun_SBS<K>_<target>.csv`
files (scaler charge + live time) and joins them with the online-charge list.
**Run:** Jupyter, from `DB/`; set `targetch` (and `replay_pass`) in the second step.
**Inputs:** `get_charge_from_Bobs_script_gmn_all_tg.csv`;
`<target>/beamCh_n_daqLvTm_prun_SBS<K>_<target>.csv` for K = 4, 7, 11, 14, 8, 9.
**Outputs:** `DB/<target>/all_beamCh_n_daqLvTm_prun_<target>.csv` (`runnum, totcharge_bob(C),
totcharge_pd(C), daqlvtm_bbhi, daqlvtm_ts1`).
**Feeds:** `DB/join_n_merge_rLists.ipynb`.

### `DB/join_n_merge_rLists.ipynb`
**What it does:** Two independent parts. (1) Beam energy: concatenates
`<target>/avg_ebeam_prun_sbs<K>_<target>.csv` for all configs, multiplies mean and std by
**1.003** (HALLA_p calibration), prints the average and spread per config, and merges the
columns into the run list. (2) Scalers: merges
`<target>/all_beamCh_n_daqLvTm_prun_<target>.csv` into the run list (on `runnum`).
**Run:** Jupyter, from `DB/`; set `target` and `replay_pass` at the top of each part. **Restart
the kernel before re-running**, otherwise the 1.003 factor is applied twice. The final writing
cells of part (1) are commented out on purpose.
**Inputs:** the files above and `good_runList_GMn_nTPE_<target>_pass_<P>.csv`.
**Outputs:** updated `good_runList_GMn_nTPE_<target>_pass_<P>.csv` (overwritten in place - commit
before and diff after).
**Feeds:** every analysis script (through `util_pd::ReadRunList`).
**Notes:** The notebook's third section ("Merging Scaler Tree Values", per config, Optics) is
marked obsolete inside the notebook.

### `DB/get_run_list.ipynb`
**What it does:** Prints the run numbers of the pass-2 LH2, LD2, Optics and Dummy lists as
space-separated bash arrays (used to fill `get_charge_from_Bobs_script.sh`).
**Run:** Jupyter, from `DB/`.
**Inputs:** `good_runList_GMn_nTPE_<target>_pass_2.csv`.
**Outputs:** printout only.
**Feeds:** `get_charge_from_Bobs_script.sh`.
