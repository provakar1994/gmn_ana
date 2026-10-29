/* Extracts GMn from the data/MC n/p ratio fit parameter R_fit and propagates errors.
   Inputs : etheta (rad), Q2 (GeV2), R_fit, R_fit_stat, R_fit_sys (one entry per
            systematic source; all sources assumed uncorrelated) and, optionally,
            names of the systematic sources (same order as R_fit_sys).
   Outputs: Corrected Born CS ratio R (n/p) w/ errors, GMn & GMn/(mun*GD) for the
            primary and alternate EMFF models w/ propagated errors, sigma_p and GEn
            (model + fit-error) systematics on the primary result, and CSV summaries.
   Usage  : root -l -b -q 'extractGMn.cpp(0.6, 4.5, 1.02, 0.01, {0.003,0.002}, {"inel","rc"})'
*/
#include <iostream>
#include <algorithm>

#include "EMFFFits.h"

// GMn/(mun*GD) and propagated errors for one choice of (GEp, GMp, GEn)
struct GMnResult {
  std::string name;
  double GMn;                 // GMn
  double val;                 // GMn/(mun*GD)
  double stat;                // stat error on GMn/(mun*GD)
  std::vector<double> sys;    // exp. sys error per source on GMn/(mun*GD)
  double sysTot;              // exp. sys errors added in quadrature
};

GMnResult calcGMn(std::string const& name, double const etheta, double const Q2,
		  double const GEp, double const GMp, double const GEn,
		  double const R, double const dR_stat, std::vector<double> const& dR_sys) {
  double const munGD = EMFFFits::GetGDip(Q2)*constant::mun;
  GMnResult r;
  r.name = name;
  r.GMn = kine::ExtractGMn(etheta,Q2,GEp,GMp,GEn,R);
  r.val = r.GMn/munGD;
  r.stat = kine::CalcGMnError(etheta,Q2,GEp,GMp,GEn,R,dR_stat)/std::abs(munGD);
  double sum2 = 0;
  for (double const dR : dR_sys) {
    double const err = kine::CalcGMnError(etheta,Q2,GEp,GMp,GEn,R,dR)/std::abs(munGD);
    r.sys.push_back(err);
    sum2 += err*err;
  }
  r.sysTot = std::sqrt(sum2);
  return r;
}

double quadSum(std::vector<double> const& v) {
  double sum2 = 0;
  for (double const x : v) sum2 += x*x;
  return std::sqrt(sum2);
}

int extractGMn(double const etheta, double const Q2, double const R_fit, double const R_fit_stat,
	       std::vector<double> const& R_fit_sys, std::vector<std::string> sys_names = {}) {

  // name the systematic sources generically if names are not given
  if (sys_names.empty())
    for (size_t i = 0; i < R_fit_sys.size(); ++i) sys_names.push_back(Form("sys%zu",i+1));
  if (sys_names.size() != R_fit_sys.size()) {
    std::cerr << "ERROR: # of sys names (" << sys_names.size() << ") != # of sys errors (" << R_fit_sys.size() << ")\n";
    return 1;
  }

  // EMFF fits
  Ye2017 yefit;
  Kelly2004 kellyfit;
  Christy2022 christyfit;
  Arrington2007 arfit;
  Galster1971 galfit;

  // -- Born CS ratio (n/p) --
  // MC ratio, using exactly the EMFFs of the MC generator (SIMC, bigbite branch)
  double const R_MC = kine::sigmaBorn_ratio_MC(etheta,Q2);
  // Assuming that the data/MC discrepancy in n/p ratio is entirely attributable to the
  // Born CS ratio, the corrected ratio is R_MC*R_fit. Errors scale the same way.
  double const R = R_MC * R_fit;
  double const dR_stat = std::abs(R_MC) * R_fit_stat;
  std::vector<double> dR_sys;
  for (double const s : R_fit_sys) dR_sys.push_back(std::abs(R_MC) * s);
  double const dR_sysTot = quadSum(dR_sys);
  double const dR_tot = std::sqrt(dR_stat*dR_stat + dR_sysTot*dR_sysTot);

  // -- GMn for various EMFF models: (GEp & GMp model)-(GEn model) --
  double const GEp_ye = yefit.GetFF(G_t::kGEp,Q2), GMp_ye = yefit.GetFF(G_t::kGMp,Q2), GEn_ye = yefit.GetFF(G_t::kGEn,Q2);
  // Primary: Ye 2017
  GMnResult const prim = calcGMn("ye-ye",etheta,Q2,GEp_ye,GMp_ye,GEn_ye,R,dR_stat,dR_sys);
  // Vary sigma_p w.r.t. primary
  std::vector<GMnResult> const vary_sigp = {
    calcGMn("ke-ye",etheta,Q2,kellyfit.GetFF(G_t::kGEp,Q2),kellyfit.GetFF(G_t::kGMp,Q2),GEn_ye,R,dR_stat,dR_sys),
    calcGMn("ar-ye",etheta,Q2,arfit.GetFF(G_t::kGEp,Q2),arfit.GetFF(G_t::kGMp,Q2),GEn_ye,R,dR_stat,dR_sys)
  };
  // Vary GEn w.r.t. primary
  std::vector<GMnResult> const vary_GEn = {
    calcGMn("ye-ke",etheta,Q2,GEp_ye,GMp_ye,kellyfit.GetFF(G_t::kGEn,Q2),R,dR_stat,dR_sys),
    calcGMn("ye-ga",etheta,Q2,GEp_ye,GMp_ye,galfit.GetFF(G_t::kGEn,Q2),R,dR_stat,dR_sys)
  };
  // No TPE check: effective proton FFs (Christy 2022), not included in the sys budget
  GMnResult const noTPE = calcGMn("ch-ye",etheta,Q2,christyfit.GetFFNoTPE(G_t::kGEp,Q2),christyfit.GetFFNoTPE(G_t::kGMp,Q2),
				  GEn_ye,R,dR_stat,dR_sys);

  // -- Systematics on the primary result beyond the exp. ones --
  // Two sources: (1) reduced proton CS (sigma_p) and (2) GEn parametrization. Each gets
  // a model-dependence part (largest deviation from primary within the family) and a
  // Ye fit-error part (shift the input FF by its fit error and re-extract), added in quadrature.
  auto maxDev = [&prim](std::vector<GMnResult> const& vars) {
    double dev = 0;
    for (auto const& v : vars) dev = std::max(dev, std::abs(v.val - prim.val));
    return dev;
  };
  double const munGD = EMFFFits::GetGDip(Q2)*constant::mun;
  auto shiftDev = [&](double const GEp, double const GMp, double const GEn) {
    return std::abs(kine::ExtractGMn(etheta,Q2,GEp,GMp,GEn,R)/munGD - prim.val);
  };
  double const dGEp_ye = yefit.GetFFwErr(G_t::kGEp,Q2)[1];
  double const dGMp_ye = yefit.GetFFwErr(G_t::kGMp,Q2)[1];
  double const dGEn_ye = yefit.GetFFwErr(G_t::kGEn,Q2)[1];
  // (1) sigma_p
  double const sys_sigp_model = maxDev(vary_sigp);
  double const sys_sigp_fitGEp = shiftDev(GEp_ye+dGEp_ye,GMp_ye,GEn_ye);
  double const sys_sigp_fitGMp = shiftDev(GEp_ye,GMp_ye+dGMp_ye,GEn_ye);
  double const sys_sigp_fit = quadSum({sys_sigp_fitGEp,sys_sigp_fitGMp});
  double const sys_sigp = quadSum({sys_sigp_model,sys_sigp_fit});
  // (2) GEn
  double const sys_GEn_model = maxDev(vary_GEn);
  double const sys_GEn_fit = shiftDev(GEp_ye,GMp_ye,GEn_ye+dGEn_ye);
  double const sys_GEn = quadSum({sys_GEn_model,sys_GEn_fit});
  double const sys_theo = quadSum({sys_sigp,sys_GEn});
  double const sysTot = std::sqrt(prim.sysTot*prim.sysTot + sys_theo*sys_theo);
  double const totErr = std::sqrt(prim.stat*prim.stat + sysTot*sysTot);

  // -- Print results --
  std::cout << "\n==========\n";
  std::cout << Form("etheta: %.4f rad (%.2f deg), Q2: %.3f GeV2\n",etheta,etheta*TMath::RadToDeg(),Q2);
  std::cout << Form("R_fit: %.4f +/- %.4f (stat) +/- %.4f (sys), R_MC: %.4f\n",R_fit,R_fit_stat,quadSum(R_fit_sys),R_MC);
  std::cout << Form("n/p Born ratio R: %.4f +/- %.4f (stat) +/- %.4f (sys), tot err: %.4f\n",R,dR_stat,dR_sysTot,dR_tot);
  for (size_t i = 0; i < dR_sys.size(); ++i)
    std::cout << Form("   %-8s: %.4f\n",sys_names[i].c_str(),dR_sys[i]);

  std::cout << "\nGMn/(mun*GD) +/- stat +/- exp. sys  [model: (GEp,GMp)-(GEn)]\n";
  std::cout << "----------\n";
  auto printRes = [](GMnResult const& r, std::string const& tag) {
    std::cout << Form("%-6s GMn: %.5f, GMn/(mun*GD): %.4f +/- %.4f +/- %.4f   %s\n",
		      r.name.c_str(),r.GMn,r.val,r.stat,r.sysTot,tag.c_str());
  };
  printRes(prim,"<- primary");
  for (auto const& v : vary_sigp) printRes(v,"(vary sigma_p)");
  for (auto const& v : vary_GEn) printRes(v,"(vary GEn)");
  printRes(noTPE,"(no TPE check)");

  std::cout << "\nPrimary error budget on GMn/(mun*GD):\n";
  std::cout << "----------\n";
  std::cout << Form("   stat        : %.4f\n",prim.stat);
  for (size_t i = 0; i < prim.sys.size(); ++i)
    std::cout << Form("   %-12s: %.4f\n",sys_names[i].c_str(),prim.sys[i]);
  std::cout << Form("   exp. sys    : %.4f\n",prim.sysTot);
  std::cout << Form("   sigma_p     : %.4f (model %.4f, Ye fit %.4f [GEp %.4f, GMp %.4f])\n",
		    sys_sigp,sys_sigp_model,sys_sigp_fit,sys_sigp_fitGEp,sys_sigp_fitGMp);
  std::cout << Form("   GEn         : %.4f (model %.4f, Ye fit %.4f)\n",sys_GEn,sys_GEn_model,sys_GEn_fit);
  std::cout << Form("   theo. sys   : %.4f\n",sys_theo);
  std::cout << Form("   tot sys     : %.4f\n",sysTot);
  std::cout << Form("   tot err     : %.4f\n",totErr);
  std::cout << Form("=> GMn/(mun*GD) = %.4f +/- %.4f (stat) +/- %.4f (exp) +/- %.4f (sigma_p) +/- %.4f (GEn)\n",
		    prim.val,prim.stat,prim.sysTot,sys_sigp,sys_GEn);

  // -- CSV summaries --
  // primary
  std::string header = "Q2,R,dR_stat,dR_sys,GMn_ov_munGD,stat,sys_exp,sys_sigp,sys_sigp_model,sys_sigp_fit,"
    "sys_GEn,sys_GEn_model,sys_GEn_fit,sys_tot,tot";
  std::string row = Form("%.3f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f",
			 Q2,R,dR_stat,dR_sysTot,prim.val,prim.stat,prim.sysTot,sys_sigp,sys_sigp_model,sys_sigp_fit,
			 sys_GEn,sys_GEn_model,sys_GEn_fit,sysTot,totErr);
  for (size_t i = 0; i < prim.sys.size(); ++i) {
    header += "," + sys_names[i];
    row += Form(",%.4f",prim.sys[i]);
  }
  std::cout << "\nSummary (primary):\n" << header << "\n" << row << "\n";
  // all models
  std::cout << "\nSummary (models):\nmodel,GMn_ov_munGD,stat,sys_exp\n";
  std::vector<GMnResult> all = {prim};
  all.insert(all.end(),vary_sigp.begin(),vary_sigp.end());
  all.insert(all.end(),vary_GEn.begin(),vary_GEn.end());
  all.push_back(noTPE);
  for (auto const& r : all)
    std::cout << Form("%s,%.4f,%.4f,%.4f\n",r.name.c_str(),r.val,r.stat,r.sysTot);
  std::cout << "\n";

  return 0;
}
