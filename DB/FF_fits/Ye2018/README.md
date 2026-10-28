## Contents

- `1-s2*` are from the publication (Z. Ye et al., PLB 777 (2018) 8).
- `proton_lookup_table.zip` is from private communication (Z. Ye). It contains five
  variants of the final proton fit (`*_t0fix7.dat` is the default fit).
- `Round6_all_sumrules_leastsq_Q21000_z12_gb5_t0fix7_noExtraTPE.dat` is extracted
  (unchanged) from `proton_lookup_table.zip`. It holds a 2000-row lookup table, the
  z-expansion coefficients of GEp and GMp, and the 16x16 covariance matrix of the
  free parameters (a1..a8 of GEp, then a1..a8 of GMp).
- `ye2018_wTPE_*` are CSV files extracted from `1-s2*` (`supp*pdf` as well).
- `ye2018_noExtraTPE_proton_lookup.csv` is the lookup table of the `noExtraTPE` `.dat`
  file above, converted to CSV.

## WARNING: "noExtraTPE" is NOT a TPE-uncorrected fit

Ye et al. apply the SIFF-Blunden TPE correction to *all* cross-section data. On top of
that, the default fit applies a phenomenological "extra" TPE term at high Q2
(PLB 777, Eq. 4): delta_2gamma -> delta_2gamma + 0.01 (eps - 1) ln(Q2)/ln(2.2).
The `noExtraTPE` fit only drops this extra term, so it is still Blunden-TPE corrected.
Checks:
- sigma_R(noExtraTPE)/sigma_R(wTPE) is 0.990-0.997 for Q2 = 3-13.5 GeV2 and nearly
  independent of eps (a TPE-uncorrected fit would differ by ~5% at low eps).
- mu_p GEp/GMp(noExtraTPE) = 0.66 at Q2 = 3 GeV2, i.e. the polarization value.

A truly TPE-uncorrected Ye fit has been requested from J. Arrington (2026-10-28).
Until it arrives, `Ye2017::GetFFNoTPE` uses `ye2018_noExtraTPE_proton_lookup.csv`
as a placeholder (and prints a warning).

## Notes on the noExtraTPE errors

- The `dGEp`/`dGMp` columns equal exactly the errors from the covariance matrix
  (fit covariance only, no data-tension term; the wTPE errors also include data
  tension and TPE uncertainty).
- The GEp-GMp correlation at Q2 = 3-13.5 GeV2 is -0.23 to +0.08. Neglecting it
  (adding GE and GM errors in quadrature) overestimates the sigma_R,p error by
  0-8% at the GMn kinematics. The GMn analysis neglects it (conservative choice).
