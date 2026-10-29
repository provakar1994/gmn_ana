# GMn extraction

Scripts that turn the measured quasi-elastic n/p cross-section ratio into GMn, for
several EMFF models, with full error propagation. The physics and error propagation
live in the library (`gmnx::` in `include/GMnExtraction.h`, `kine::` in
`include/KinematicVar.h`); the scripts here are thin wrappers.

| Script | What it does |
|---|---|
| `extractGMn.cpp` | One kinematic point: takes R_fit (data/MC ratio fit parameter) and its errors, prints R, GMn for all models, the error budget, and a CSV block. |
| `extractGMn_csv.cpp` | All kinematic points of a results CSV (e.g. `scripts/lta/results.csv`, which holds R itself): writes a CSV table with GMn and all errors for every point and model. |

## Running

Build and install `gmn_ana` first, then from this directory:

```
root -l -b -q ../../../install_gmn_ana/runenv/rootlogon.C \
  'extractGMn.cpp(4.5, 0.6833, 1.02, 0.01, {0.004,0.003,0.002}, {"Inel","NDE1","RC_sc"})'
```

Arguments: `Q2` (GeV²) and `eps` (acceptance averaged), `R_fit`, `R_fit_stat`, the
systematic errors on `R_fit` (one per source, all uncorrelated) and, optionally, their
names. A source whose name ends in **`_sc`** is a **scale** error; all others are
**point-to-point**.

### All points from a results CSV

```
root -l -b -q ../../../install_gmn_ana/runenv/rootlogon.C \
  'extractGMn_csv.cpp("../lta/results.csv", "gmn_results.csv")'
```

**Input** (one row per kinematic point; columns are found by name in the header):

| Columns | Meaning |
|---|---|
| `q2`, `eps` | acceptance-averaged Q² (GeV²) and ε (required) |
| `R`, `R_stat` | n/p Born cross-section ratio and its stat error (required) |
| every column after `R_stat` | a systematic error on R (absolute; `*_sc` = scale, others p2p). Empty → 0, with a warning |
| every other column before `R_stat` | label (e.g. `conf`, `field`, `etheta`), copied to the output |

**Output**: one row per (kinematic point × model), all errors on GMn/(μ_n G_D):

```
<labels>,Q2,eps,etheta_eff,R,R_stat,Model,TPE,Val,StatErr,
<p2p sources>,pCS,GEn,<scale sources>,TotSys_p2p,TotSys_sc,StatP2P,TotErr
```

| Column | Meaning |
|---|---|
| `etheta_eff` | effective e⁻ angle (rad) reproducing ε at Q² (information only, see below) |
| `Model`, `TPE` | model name and TPE treatment of its proton FFs (see Models) |
| `Val`, `StatErr` | GMn/(μ_n G_D) and its stat error |
| p2p sources | one column per point-to-point systematic; `X` is written as `X_p2p` if a scale source `X_sc` also exists (e.g. `NDE2_p2p`, `RC_p2p`) |
| `pCS`, `GEn` | FF systematics (same for all models of a point) |
| scale sources | one column per `*_sc` systematic |
| `TotSys_p2p`, `TotSys_sc`, `StatP2P`, `TotErr` | totals, see Error propagation |

## Kinematics: (Q², ε)

Each point is defined by the **acceptance-averaged Q² and ε** (ε = ε_p, the virtual photon
polarization for the proton), which are taken as the true kinematics. The Born cross
sections, and so R and GMn, depend on the kinematics only through Q² and ε. The FFs depend
on Q² only.

The `kine::` functions are written for per-event use and take (θ_e, Q²), from which they
compute both ε_p and ε_n (which differ slightly at the same θ_e because of the nucleon mass
difference). `gmnx` therefore converts ε to the effective electron angle that gives
ε_p = ε at Q² (`kine::etheta_from_epsilon`),

```
tan²(θ_eff/2) = (1/ε − 1) / (2(1 + τ_p))
```

and uses θ_eff for all calculations; ε_n follows consistently from it. This is exact, not an
approximation. θ_eff is printed for information only.

## From R_fit to R

The data/MC fit returns R_fit, the factor by which the MC n/p ratio must be scaled to
match the data. Assuming the data/MC difference is entirely due to the Born cross-section
ratio,

```
R = σ_n/σ_p (Born) = R_MC · R_fit,      δR = R_MC · δR_fit   (for every error)
```

R_MC is the Born ratio of the MC generator, computed with exactly the EMFFs SIMC
(bigbite branch, `physics_proton.f`) uses: Kelly GEp, GMp, GMn and the Riordan/"Seamus"
GEn, with τ = Q²/4M_p² and GMp(0) = 2.79, GMn(0) = −1.913 (`kine::sigmaBorn_ratio_MC`,
checked by `tests/test_simc_rmc.cpp`), evaluated at (Q², θ_eff).

## Models

Model names are **"A‑B": the proton FFs (and so the reduced proton cross section σ_R,p)
come from fit A, GEn from fit B.** "noTPE" marks TPE-uncorrected (effective) proton FFs.

| # | Model | σ_R,p from | GEn from | Proton FFs |
|---|---|---|---|---|
| 1 | Ch noTPE‑Ye | Christy 2022 | Ye 2018 | TPE uncorrected |
| 2 | **Ye noTPE‑Ye** (primary) | Ye 2018, no TPE | Ye 2018 | TPE uncorrected |
| 3 | Ar noTPE‑Ye | Arrington 2007, Table IV | Ye 2018 | TPE uncorrected |
| 4 | Ye noTPE‑Ke | Ye 2018, no TPE | Kelly 2004 | TPE uncorrected |
| 5 | Ye noTPE‑Ga | Ye 2018, no TPE | Galster | TPE uncorrected |
| 6 | Ye‑Ye | Ye 2018 | Ye 2018 | TPE corrected (reference) |
| 7 | Ar‑Ye | Arrington 2007 | Ye 2018 | TPE corrected (reference) |
| 8 | Ke‑Ye | Kelly 2004 | Ye 2018 | mixed (reference) |

**Why TPE-uncorrected proton FFs.** No TPE correction is applied to the neutron (it is not
known well enough), so the measured n/p ratio is TPE-uncorrected for both nucleons. Using
TPE-uncorrected proton FFs keeps σ_R,p consistent with that: GMn then contains the neutron
TPE effect, and readers can apply their own neutron TPE correction. The ratios themselves
are model independent and are published as well.

**Why Ye as the primary.** Ye et al. use the same (older) radiative-correction prescription
as SIMC. Christy 2022 includes the GMp12 data and is the better fit, but uses an updated RC
prescription, which would be inconsistent with our MC. Christy is shown on top for that
reason.

Rows 6–8 are for a visual check of model dependence only; they are not used for any
systematic.

### Caveats on the fits

- **Ye no-TPE is a placeholder.** The file used (`DB/FF_fits/Ye2018/*noExtraTPE*`) only
  drops Ye's phenomenological "extra" TPE term; the SIFF-Blunden TPE correction is still
  applied, so it is effectively TPE corrected. A truly TPE-uncorrected Ye fit has been
  requested (J. Arrington, 2026-10-28). The library prints a warning while the placeholder
  is in use. See `DB/FF_fits/Ye2018/README.md`.
- **Christy 2022** (PRL 128, 102002, Eq. 3) is a fit to cross sections only, without hard
  TPE corrections, so it is TPE uncorrected (`Christy2022::GetFFNoTPE`). It uses an updated
  RC prescription (see above).
- **Arrington 2007** no-TPE (PRC 76, 035205, Appendix A, Table IV) is a fit to TPE-uncorrected
  cross sections only. Its GE is constrained only up to Q² ≈ 6 GeV², so at higher Q² GEp is an
  extrapolation (a small effect on σ_R,p, where τ·GMp² dominates).
- **Kelly 2004 is "mixed"** (PRC 70, 068202). Above Q² = 1 GeV² Kelly's GEp comes from
  recoil-polarization GEp/GMp data (Rosenbluth GEp is omitted), so it is close to the **Born**
  value, while GMp comes from cross-section data **without** TPE corrections ("Nor have
  corrections … for exchange of two hard photons been applied"), so it is an **effective**
  FF. σ_R,p = ε·GEp² + τ·GMp² built from them is neither the Born cross section (GMp still
  contains TPE) nor the measured one (the ε-dependence that TPE adds to the measured cross
  section was absorbed into a larger Rosenbluth GE, which Kelly replaced by the smaller
  polarization GE). This is the GE/GM discrepancy itself, with Kelly on both sides of it. At
  our Q², τ·GMp² dominates σ_R,p, so Kelly behaves mostly like an (older) uncorrected fit.
  (SIMC also uses Kelly; that is irrelevant for R_MC, which only has to match the generator.)
- **Galster** GEn = −μ_n τ/(1 + 5.6τ)·G_D with τ = Q²/4M_n² (Perdrisat et al., Prog. Part.
  Nucl. Phys. 59, 694, Eq. 16). An old parametrization, shown for comparison only.

## Error propagation

All errors are quoted on GMn/(μ_n G_D) and all sources are uncorrelated.

| Error | Source | Group |
|---|---|---|
| stat | δR_stat | — |
| each experimental systematic | δR of that source | p2p, or scale if the name ends in `_sc` |
| **pCS** | error of σ_R,p from the proton FF fit errors | p2p |
| **GEn** | GEn fit error | p2p |

```
TotSys_p2p = ⊕ (p2p experimental sources, pCS, GEn)
TotSys_sc  = ⊕ (scale sources)
StatP2P    = stat ⊕ TotSys_p2p
TotErr     = stat ⊕ TotSys_p2p ⊕ TotSys_sc          (⊕ = sum in quadrature)
```

Stat and the experimental systematics are propagated **per model**, with that model's FFs.
pCS and GEn are computed **once, at the primary model**, and the same values are used for
all models. They come from the Ye fit errors, the only fits with readily available
uncertainties; the spread between models is *not* used as a systematic (some fits are old
and not independent).

### 1. R as a function of GMn

For a nucleon N = p, n at the same electron angle θ and Q²:

```
τ_N = Q²/(4 M_N²)
ε_N = [1 + 2(1 + τ_N) tan²(θ/2)]⁻¹
σ_R,N = ε_N G_E,N² + τ_N G_M,N²                 (reduced cross section)
σ_N  ∝ σ_Mott · σ_R,N / [ε_N (1 + τ_N)]           (Rosenbluth)
```

The Mott cross section is the same for p and n, so it cancels in the ratio:

```
R = [σ_R,n / (ε_n(1+τ_n))] / [σ_R,p / (ε_p(1+τ_p))]
```

### 2. Inverting for GMn

With term1 = ε_n(1+τ_n) / [ε_p(1+τ_p)]:

```
σ_R,n = term1 · σ_R,p · R
GMn²  = (term1 · σ_R,p · R − ε_n GEn²) / τ_n          (GMn < 0: negative root)
```

This is `kine::ExtractGMn`.

### 3. First-order propagation

GMn = f(R, σ_R,p, GEn). For uncorrelated inputs x_i with errors δx_i,

```
δGMn² = Σ_i (∂f/∂x_i)² δx_i²
```

so there is one term per uncertain input: R (measurement), σ_R,p (proton FF fit) and GEn
(neutron FF fit). They are independent, so all covariance terms vanish. With F ≡ GMn² and
the chain rule, ∂GMn/∂x = (∂F/∂x)/(2 GMn):

```
∂GMn/∂R      = term1 · σ_R,p / (2 τ_n GMn)
∂GMn/∂σ_R,p  = term1 · R     / (2 τ_n GMn)
∂GMn/∂GEn    = −ε_n GEn / (τ_n GMn)
```

Each error is |∂GMn/∂x|·δx. The sign only gives the direction of the shift (a larger σ_R,p
gives a larger |GMn|, a larger |GEn| a smaller one); it drops out in the quadrature sum and
would only matter for correlated sources. So:

```
dGMn_R    = term1 · σ_R,p · δR    / (2 τ_n |GMn|)        kine::CalcGMnError
dGMn_pCS  = term1 · R · δσ_R,p    / (2 τ_n |GMn|)        kine::CalcGMnError_sigmaRp
dGMn_GEn  = ε_n · |GEn| · δGEn    / (τ_n |GMn|)          kine::CalcGMnError_GEn
```

### 4. The error on σ_R,p

GMn depends on GEp and GMp only through σ_R,p, so their errors enter via σ_R,p:

```
∂σ_R,p/∂GEp = 2 ε_p GEp,      ∂σ_R,p/∂GMp = 2 τ_p GMp
δσ_R,p² = (2ε_p GEp δGEp)² + (2τ_p GMp δGMp)² + 2ρ (2ε_p GEp δGEp)(2τ_p GMp δGMp)
```

We take **ρ = 0** (`kine::sigmaReducedError`), as suggested by J. Arrington: GE and GM are
anti-correlated in an L/T separation, so this overestimates δσ_R,p, which is the
conservative choice, and global-fit errors themselves should be treated conservatively.
Cross-check with the full covariance matrix of the Ye noExtraTPE fit: ρ(GE, GM) = −0.23 to
+0.08 at our kinematics, and ρ = 0 overestimates δσ_R,p by only 0–8%.

### 5. Relative forms

Using term1·σ_R,p·R = τ_n GMn² + ε_n GEn² and S ≡ ε_n GEn² / (τ_n GMn²) (the electric
fraction of σ_R,n, a few % at our kinematics):

```
dGMn_pCS / |GMn| = ½ (1 + S) · δσ_R,p / σ_R,p
dGMn_GEn / |GMn| = S · δGEn / |GEn|
```

GMn ∝ √σ_R,p, so the relative σ_R,p error enters GMn halved. The GEn error is suppressed by
S, which is why even a large relative GEn error gives a moderate GMn error. All of the above
is checked against numerical derivatives in `tests/test_gmn_extraction.cpp` and
`tests/test_gmn_module.cpp`.

### 6. Which FF errors are used

- **pCS**: δGEp, δGMp of the Ye no-TPE fit (the `dGEp`, `dGMp` columns of the noExtraTPE
  table). **These are fit-covariance errors only.** They do not include Ye's data-tension
  term (the Mainz vs world-data tension, extended to high Q² with an ad hoc power law),
  which for the TPE-corrected fit dominates the GMp error up to Q² ≈ 10 GeV² (≈0.7–0.95% vs
  0.2–0.6% from the covariance; rough estimate). pCS is therefore likely underestimated
  (by a factor of up to ~2–4 at low Q²). To be revisited when the TPE-uncorrected Ye fit
  arrives; the wTPE errors are deliberately not mixed in.
- **GEn**: the full tabulated Ye error dG (`Ye2017::GetFFwErr`), not Ye's smooth
  parametrization of it (dG_Par, `Ye2017::GetFFwParErr`), which deviates from dG by up to
  ~15%. For the neutron, Ye's errors are from the fit covariance; data-set tensions were
  handled by the data selection.

### Possible improvement

pCS and GEn move coherently with Q² between models. Following J. Arrington, they are
treated as point-to-point for now; if their impact matters, they could be split into a
normalization part plus a point-to-point part (as done for other data sets in global
fits). With only a few widely separated Q² points, it is not clear how correlated they are.

## Open items (as of 2026-10-29)

- [ ] **Truly TPE-uncorrected Ye fit** (requested from J. Arrington). When it arrives: convert it
      to CSV in `DB/FF_fits/Ye2018/`, point `Ye2017::GetFFNoTPE` at it, remove the placeholder
      warning, and check what its errors contain (fit covariance only, or also Ye's data-tension
      term, which is missing from pCS now; see "Which FF errors are used").
- [ ] **Regenerate R with the fixed `fit_dx.cpp`** (R_MC with the SIMC FFs at (⟨Q²⟩, ⟨ε⟩)), then
      rebuild `scripts/lta/results.csv`. Until then `scripts/lta/results_RMCfix.csv` is used: the old
      R rescaled by R_MC(new)/R_MC(old), with R_MC(old) reconstructed from the pre-fix code
      (`4a3a823`; validated against the one tracked fit_dx output). This assumes every row was made
      with that code at the row's (etheta, q2).
- [ ] **Fill the `Cut V` systematic** (empty in `results.csv`, read as 0 with a warning).
- [ ] **Possibly split pCS and GEn** into a normalization plus a point-to-point part (J. Arrington's
      suggestion; see "Possible improvement"), if their impact matters for the global fit.
