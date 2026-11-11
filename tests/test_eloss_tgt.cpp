/*
  Unit test: target energy-loss correction util_pd::GetElossInTgt.
  ----
  The function has two overloads:
  - core:       the shield beside the scattering chamber is an explicit input (used for MC);
  - run number: the shield is deduced from the run number (used for data) and the core is called.
  Checks, for LH2 and LD2 at SBS-11 (the only setting with shields), at a few (vz, etheta) points
  in each of the three run periods (no shield, PE shield, Al shield):
  1. The run-number overload equals the core with the matching shield (exactly).
  2. Both equal the closed-form energy loss written out below, so the split into two
     overloads provably did not change the data correction.
  3. vz outside the target is clamped to the target ends.
*/
#include "Utilities.h"
#include "ExpConstants.h"
#include "TestUtils.h"
#include <cmath>
#include <string>
#include <vector>

namespace {

  int const kSBSConf = 11;   // SBS-11: the only setting where shields were installed
  int const kVerbose = 0;

  // One run number inside each shield period
  struct RunPeriod { int rnum; expconst::ShieldType shield; std::string name; };
  std::vector<RunPeriod> const kRunPeriods = {
    {expconst::firstrun_PE_shield - 10, expconst::ShieldType::kNone, "no shield"},
    {expconst::firstrun_PE_shield + 10, expconst::ShieldType::kPE,   "PE shield"},
    {expconst::firstrun_Al_shield + 10, expconst::ShieldType::kAl,   "Al shield"}
  };

  // Closed-form energy loss (GeV), written out from the formula in Utilities.cpp:
  // before scattering: target liquid up to the vertex + Al upstream window;
  // after scattering:  half a cell diameter of liquid + the Al cell wall (both along the e- path)
  //                    + the shield, if any.
  void ExpectedEloss(std::string const& target, double vz, double etheta,
		     expconst::ShieldType shield, double& Eloss_bs, double& Eloss_as) {
    bool const is_lh2 = (target == "LH2");
    // D2 has Z/A = 0.5 instead of 1 for H2
    double const ZoverA_scale = is_lh2 ? 1.0 : 0.5;
    double const rho      = is_lh2 ? expconst::lh2_TgtRho    : expconst::ld2_TgtRho;
    double const uwin     = is_lh2 ? expconst::lh2_uWinThick : expconst::ld2_uWinThick;
    double const celldiam = is_lh2 ? expconst::lh2_CellDiam  : expconst::ld2_CellDiam;
    double const cellwall = is_lh2 ? expconst::lh2_CellThick : expconst::ld2_CellThick;

    double const dEdx_tgt_bs = expconst::GetdEdxCollH(kSBSConf, 1)*ZoverA_scale;
    double const dEdx_tgt_as = expconst::GetdEdxCollH(kSBSConf, 0)*ZoverA_scale;
    double const dEdx_Al_bs  = expconst::GetdEdxCollAl(kSBSConf, 1);
    double const dEdx_Al_as  = expconst::GetdEdxCollAl(kSBSConf, 0);
    double const dEdx_PE_as  = expconst::GetdEdxCollPE(kSBSConf, 0);

    // path length in the target before the vertex (cm), clamped to the target
    double const path_bs = std::min(std::max(vz*1E2 + 0.5*expconst::tgtlen, 0.0), expconst::tgtlen);
    double const sin_th = std::sin(etheta);

    Eloss_bs = path_bs*rho*dEdx_tgt_bs + uwin*expconst::Al_Rho*dEdx_Al_bs;
    Eloss_as = 0.5*celldiam/sin_th*rho*dEdx_tgt_as + cellwall/sin_th*expconst::Al_Rho*dEdx_Al_as;
    if (shield == expconst::ShieldType::kPE)
      Eloss_as += expconst::PE_ShieldThick*expconst::PE_Rho*dEdx_PE_as;
    if (shield == expconst::ShieldType::kAl)
      Eloss_as += expconst::Al_ShieldThick*expconst::Al_Rho*dEdx_Al_as;
  }

} //::

int main() {
  double const kTolClosedForm = 1e-14; // same formula, only round-off allowed (relative)
  std::vector<std::string> const targets = {"LH2", "LD2"};
  std::vector<double> const vzs = {-0.05, 0.0, 0.06};          // m, inside the 15 cm target
  std::vector<double> const ethetas = {0.30, 0.55, 0.85};      // rad, 17-49 deg

  // 1 & 2: run-number overload == core == closed form, in every run period
  for (auto const& tgt : targets)
    for (auto const& period : kRunPeriods)
      for (double vz : vzs)
	for (double eth : ethetas) {
	  std::vector<double> byRun, byShield;
	  util_pd::GetElossInTgt(tgt, period.rnum, kSBSConf, vz, eth, kVerbose, byRun);
	  util_pd::GetElossInTgt(tgt, kSBSConf, vz, eth, period.shield, kVerbose, byShield);
	  std::string const what = tgt + ", " + period.name + Form(", vz=%.2f m, eth=%.2f rad", vz, eth);

	  CHECK_ABS(byRun[0], byShield[0], 0.0, what + ": run# == core (before)");
	  CHECK_ABS(byRun[1], byShield[1], 0.0, what + ": run# == core (after)");

	  double exp_bs, exp_as;
	  ExpectedEloss(tgt, vz, eth, period.shield, exp_bs, exp_as);
	  CHECK_REL(byShield[0], exp_bs, kTolClosedForm, what + ": closed form (before)");
	  CHECK_REL(byShield[1], exp_as, kTolClosedForm, what + ": closed form (after)");
	}

  // The shield only adds to the loss after scattering, and Al (~1.9 MeV at SBS-11) > none
  {
    std::vector<double> none, Al;
    util_pd::GetElossInTgt("LD2", kSBSConf, 0.0, 0.55, expconst::ShieldType::kNone, kVerbose, none);
    util_pd::GetElossInTgt("LD2", kSBSConf, 0.0, 0.55, expconst::ShieldType::kAl, kVerbose, Al);
    CHECK_ABS(Al[0], none[0], 0.0, "shield does not change the loss before scattering");
    CHECK_TRUE(Al[1] - none[1] > 1.5E-3 && Al[1] - none[1] < 2.5E-3, "Al shield adds ~1.9 MeV after scattering");
  }

  // 3: vz outside the target is clamped to the target ends
  {
    double const half_len = 0.5*expconst::tgtlen*1E-2; // m
    double const kTolClamp = 1e-12;  // GeV; vz at the target end scales to 0 or tgtlen only up to round-off
    std::vector<double> inside, outside;
    util_pd::GetElossInTgt("LH2", kSBSConf, -half_len, 0.55, expconst::ShieldType::kNone, kVerbose, inside);
    util_pd::GetElossInTgt("LH2", kSBSConf, -half_len - 0.1, 0.55, expconst::ShieldType::kNone, kVerbose, outside);
    CHECK_ABS(outside[0], inside[0], kTolClamp, "vz upstream of the target is clamped");
    util_pd::GetElossInTgt("LH2", kSBSConf, half_len, 0.55, expconst::ShieldType::kNone, kVerbose, inside);
    util_pd::GetElossInTgt("LH2", kSBSConf, half_len + 0.1, 0.55, expconst::ShieldType::kNone, kVerbose, outside);
    CHECK_ABS(outside[0], inside[0], kTolClamp, "vz downstream of the target is clamped");
  }

  return test::Summary("test_eloss_tgt");
}
