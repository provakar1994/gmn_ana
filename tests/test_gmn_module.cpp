/*
  Unit test for the GMn extraction module (gmnx:: in GMnExtraction.h).
  ----
  Checks, at the GMn kinematics:
  1. Model rows: names, order, TPE groups, and that the primary uses Ye no-TPE proton FFs
     and Ye GEn; Christy is only available w/o TPE.
  2. InputFromRfit: R = R_MC*R_fit (errors scaled the same way); "_sc" names are scale errors.
  3. Every row: GMn/(mu_n*GD) and the stat/experimental errors equal direct kine:: calls
     w/ that row's FFs.
  4. pCS and GEn systematics are the same on all rows and equal the relative forms
     (README): pCS/|val| = (1+S)/2 * dsigmaR_p/sigmaR_p,  GEn/|val| = S * dGEn/|GEn|,
     evaluated at the primary model w/ S = eps_n*GEn^2/(tau_n*GMn^2).
  5. Totals: p2p = p2p sources (+) pCS (+) GEn; scale = "_sc" sources; stat (+) p2p; total.
*/
#include "Constants.h"
#include "EMFFFits.h"
#include "GMnExtraction.h"
#include "KinematicVar.h"
#include "TestUtils.h"
#include <cmath>

namespace {
  double const kTol = 1e-12;   // same formulas: only floating point round-off
  double const kRfit = 1.02;   // arbitrary but realistic data/MC ratio fit parameter
  double const kRfitStat = 0.01;
}

int main() {
  // -- 0. MakeErrorSource: scale errors are marked by the "_sc" suffix only --
  CHECK_TRUE(gmnx::MakeErrorSource("RC_sc", 0.1).isScale, "RC_sc is a scale error");
  CHECK_TRUE(gmnx::MakeErrorSource("_sc", 0.1).isScale, "_sc is a scale error");
  CHECK_TRUE(!gmnx::MakeErrorSource("RC", 0.1).isScale, "RC is point-to-point");
  CHECK_TRUE(!gmnx::MakeErrorSource("sc", 0.1).isScale, "sc is point-to-point");
  CHECK_TRUE(!gmnx::MakeErrorSource("RC_scale", 0.1).isScale, "RC_scale is point-to-point");

  // -- 1. Christy is only provided w/o TPE --
  Christy2022 christy;
  CHECK_TRUE(christy.GetFF(G_t::kGMp, 4.5) == -1000, "Christy GetFF is not provided (TPE uncorrected fit)");
  CHECK_TRUE(christy.GetFFNoTPE(G_t::kGMp, 4.5) > 0, "Christy GetFFNoTPE gives GMp");

  // GMn kinematics {etheta (rad), Q2 (GeV2)} (scripts/lta/results.csv)
  double const kKine[][2] = {{0.623, 2.989}, {0.854, 4.488}, {0.465, 4.531},
			     {0.808, 7.464}, {0.694, 9.834}, {0.727, 13.465}};
  std::vector<std::string> const kOrder = {"Ch noTPE-Ye", "Ye noTPE-Ye", "Ar noTPE-Ye", "Ye noTPE-Ke",
					   "Ye noTPE-Ga", "Ye-Ye", "Ar-Ye", "Ke-Ye"};
  std::vector<gmnx::TPEGroup> const kTPE = {
    gmnx::TPEGroup::kUncorrected, gmnx::TPEGroup::kUncorrected, gmnx::TPEGroup::kUncorrected,
    gmnx::TPEGroup::kUncorrected, gmnx::TPEGroup::kUncorrected, gmnx::TPEGroup::kCorrected,
    gmnx::TPEGroup::kCorrected, gmnx::TPEGroup::kMixed};
  // experimental systematics on R_fit: two point-to-point, one scale
  std::vector<std::string> const kSysNames = {"Inel", "RC_sc", "NDE1"};
  std::vector<double> const kSysRfit = {0.004, 0.006, 0.003};

  Ye2017 ye;
  for (auto const& k : kKine) {
    double const etheta = k[0], Q2 = k[1];
    std::string const at = " at Q2 = " + std::to_string(Q2);

    // -- 1. model rows --
    std::vector<gmnx::ModelDef> const defs = gmnx::ModelDefs(Q2);
    CHECK_TRUE(defs.size() == kOrder.size(), "number of models" + at);
    for (size_t i = 0; i < defs.size() && i < kOrder.size(); ++i) {
      CHECK_TRUE(defs[i].name == kOrder[i], "model name/order: " + kOrder[i] + at);
      CHECK_TRUE(defs[i].tpe == kTPE[i], "TPE group of " + kOrder[i] + at);
    }
    gmnx::ModelDef const& prim = defs[gmnx::kPrimaryModel];
    CHECK_TRUE(prim.name == "Ye noTPE-Ye", "primary model" + at);
    CHECK_REL(prim.GEp, ye.GetFFNoTPE(G_t::kGEp, Q2), kTol, "primary GEp = Ye no-TPE" + at);
    CHECK_REL(prim.GMp, ye.GetFFNoTPE(G_t::kGMp, Q2), kTol, "primary GMp = Ye no-TPE" + at);
    CHECK_REL(prim.GEn, ye.GetFF(G_t::kGEn, Q2), kTol, "primary GEn = Ye" + at);

    // -- 2. input from R_fit --
    gmnx::Input const in = gmnx::InputFromRfit(etheta, Q2, kRfit, kRfitStat, kSysNames, kSysRfit);
    double const R_MC = kine::sigmaBorn_ratio_MC(etheta, Q2);
    CHECK_REL(in.R, R_MC*kRfit, kTol, "R = R_MC*R_fit" + at);
    CHECK_REL(in.dR_stat, R_MC*kRfitStat, kTol, "dR_stat = R_MC*R_fit_stat" + at);
    CHECK_TRUE(in.sys.size() == kSysNames.size(), "number of sys sources" + at);
    for (size_t i = 0; i < in.sys.size(); ++i)
      CHECK_REL(in.sys[i].dR, R_MC*kSysRfit[i], kTol, "dR_sys = R_MC*R_fit_sys: " + kSysNames[i] + at);

    // -- 3.-5. results --
    std::vector<gmnx::ModelResult> const res = gmnx::ExtractAll(in);
    CHECK_TRUE(res.size() == defs.size(), "one result per model" + at);
    double const absMunGD = std::abs(constant::mun*EMFFFits::GetGDip(Q2));

    // expected pCS and GEn sys from the relative forms, at the primary model
    gmnx::ModelResult const& rp = res[gmnx::kPrimaryModel];
    gmnx::FFErrors const ffErr = gmnx::PrimaryFFErrors(Q2);
    double const tau_p = kine::tau(Q2, "p"), eps_p = kine::epsilon(etheta, Q2, "p");
    double const sigRp = eps_p*prim.GEp*prim.GEp + tau_p*prim.GMp*prim.GMp;
    double const dSigRp = 2.*std::sqrt(std::pow(eps_p*prim.GEp*ffErr.dGEp, 2) + std::pow(tau_p*prim.GMp*ffErr.dGMp, 2));
    double const S = kine::epsilon(etheta, Q2, "n")*prim.GEn*prim.GEn / (kine::tau(Q2, "n")*rp.GMn*rp.GMn);
    double const pCS = 0.5*(1. + S)*dSigRp/sigRp*std::abs(rp.val);
    double const GEnSys = S*ffErr.dGEn/std::abs(prim.GEn)*std::abs(rp.val);

    for (auto const& r : res) {
      gmnx::ModelDef const& m = r.model;
      std::string const of = " of " + m.name + at;
      // 3. value, stat and experimental sys w/ the row's own FFs
      double const GMn = kine::ExtractGMn(etheta, Q2, m.GEp, m.GMp, m.GEn, in.R);
      CHECK_REL(r.GMn, GMn, kTol, "GMn" + of);
      CHECK_REL(r.val, GMn/(constant::mun*EMFFFits::GetGDip(Q2)), kTol, "GMn/(mu_n*GD)" + of);
      CHECK_REL(r.stat, kine::CalcGMnError(etheta, Q2, m.GEp, m.GMp, m.GEn, in.R, in.dR_stat)/absMunGD, kTol, "stat" + of);
      CHECK_TRUE(r.sys.size() == in.sys.size(), "one sys value per source" + of);
      for (size_t i = 0; i < r.sys.size(); ++i)
	CHECK_REL(r.sys[i], kine::CalcGMnError(etheta, Q2, m.GEp, m.GMp, m.GEn, in.R, in.sys[i].dR)/absMunGD, kTol,
		  "sys " + in.sys[i].name + of);
      // 4. FF systematics: same on all rows, from the primary model
      CHECK_REL(r.pCS, pCS, kTol*1e2, "pCS" + of);
      CHECK_REL(r.GEnSys, GEnSys, kTol*1e2, "GEn sys" + of);
      // 5. totals (sources: Inel & NDE1 point-to-point, RC_sc scale)
      double const p2p = std::sqrt(r.sys[0]*r.sys[0] + r.sys[2]*r.sys[2] + r.pCS*r.pCS + r.GEnSys*r.GEnSys);
      CHECK_REL(r.sysP2P, p2p, kTol, "p2p sys total" + of);
      CHECK_REL(r.sysScale, r.sys[1], kTol, "scale sys total" + of);
      CHECK_REL(r.statP2P, std::hypot(r.stat, p2p), kTol, "stat (+) p2p" + of);
      CHECK_REL(r.tot, std::sqrt(r.stat*r.stat + p2p*p2p + r.sys[1]*r.sys[1]), kTol, "total error" + of);
    }
  }

  return test::Summary("test_gmn_module");
}
