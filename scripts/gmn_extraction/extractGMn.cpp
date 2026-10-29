/* Extracts GMn from the data/MC n/p ratio fit parameter R_fit at one kinematic point,
   for all EMFF models of the gmnx module, w/ full error propagation.
   Inputs : Q2 (GeV2) and eps (acceptance averaged), R_fit, R_fit_stat, R_fit_sys (one entry per
            systematic source; all sources uncorrelated) and, optionally, the names of
            the systematic sources (same order as R_fit_sys). Sources whose name ends in
            "_sc" are scale errors, all others point-to-point.
   Outputs: The n/p Born CS ratio R = R_MC*R_fit w/ errors, and GMn/(mu_n*GD) for every
            model w/ its error budget (printed, plus a CSV block).
   Physics, models and error propagation: see README.md in this directory.
   Usage  : (from this directory)
            root -l -b -q ../../../install_gmn_ana/runenv/rootlogon.C \
              'extractGMn.cpp(4.5, 0.6833, 1.02, 0.01, {0.004,0.003,0.002}, {"Inel","NDE1","RC_sc"})'
*/
#include <iostream>

#include "GMnExtraction.h"
#include "KinematicVar.h"

int extractGMn(double const Q2, double const eps, double const R_fit, double const R_fit_stat,
	       std::vector<double> const& R_fit_sys, std::vector<std::string> sys_names = {}) {

  // name the systematic sources generically if names are not given
  if (sys_names.empty())
    for (size_t i = 0; i < R_fit_sys.size(); ++i) sys_names.push_back(Form("sys%zu",i+1));
  if (sys_names.size() != R_fit_sys.size()) {
    std::cerr << "ERROR: # of sys names (" << sys_names.size() << ") != # of sys errors (" << R_fit_sys.size() << ")\n";
    return 1;
  }

  // -- n/p Born CS ratio R = R_MC*R_fit, and GMn for all models --
  gmnx::Input const in = gmnx::InputFromRfit(Q2,eps,R_fit,R_fit_stat,sys_names,R_fit_sys);
  std::vector<gmnx::ModelResult> const results = gmnx::ExtractAll(in);
  gmnx::ModelResult const& prim = results[gmnx::kPrimaryModel];

  // -- Print the input --
  std::cout << "\n==========\n";
  std::cout << Form("Q2: %.3f GeV2, eps: %.4f (effective e- angle: %.4f rad = %.2f deg)\n",
		    Q2,eps,in.etheta,in.etheta*TMath::RadToDeg());
  std::cout << Form("R_fit: %.4f +/- %.4f (stat), R_MC: %.4f\n",R_fit,R_fit_stat,kine::sigmaBorn_ratio_MC(in.etheta,Q2));
  std::cout << Form("n/p Born ratio R: %.4f +/- %.4f (stat)\n",in.R,in.dR_stat);
  for (auto const& s : in.sys)
    std::cout << Form("   %-10s: %.4f %s\n",s.name.c_str(),s.dR,s.isScale ? "(scale)" : "");

  // -- Print GMn for all models --
  std::cout << "\nGMn/(mu_n*GD) for all models  [A-B: proton FFs from A, GEn from B]\n";
  std::cout << Form("%-12s %-11s %7s %7s %7s %7s %7s %7s\n","model","p TPE","value","stat","p2p","scale","stat+p2p","total");
  for (auto const& r : results)
    std::cout << Form("%-12s %-11s %7.4f %7.4f %7.4f %7.4f %8.4f %7.4f %s\n",r.model.name.c_str(),
		      gmnx::TPEGroupName(r.model.tpe).c_str(),r.val,r.stat,r.sysP2P,r.sysScale,r.statP2P,r.tot,
		      (&r == &prim) ? " <- primary" : "");

  // -- Print the error budget of the primary model --
  std::cout << Form("\nError budget of the primary model (%s) on GMn/(mu_n*GD) = %.4f:\n",prim.model.name.c_str(),prim.val);
  std::cout << Form("   %-10s: %.4f\n","stat",prim.stat);
  for (size_t i = 0; i < in.sys.size(); ++i)
    std::cout << Form("   %-10s: %.4f %s\n",in.sys[i].name.c_str(),prim.sys[i],in.sys[i].isScale ? "(scale)" : "");
  std::cout << Form("   %-10s: %.4f\n","pCS",prim.pCS);
  std::cout << Form("   %-10s: %.4f\n","GEn",prim.GEnSys);
  std::cout << Form("   %-10s: %.4f\n","p2p sys",prim.sysP2P);
  std::cout << Form("   %-10s: %.4f\n","scale sys",prim.sysScale);
  std::cout << Form("   %-10s: %.4f\n","stat+p2p",prim.statP2P);
  std::cout << Form("   %-10s: %.4f\n","total",prim.tot);

  // -- CSV block: one row per model --
  std::string header = "Model,TPE,Q2,eps,R,Val,StatErr";
  for (auto const& s : in.sys) header += "," + s.name;
  header += ",pCS,GEn,TotSys_p2p,TotSys_sc,StatP2P,TotErr";
  std::cout << "\nCSV:\n" << header << "\n";
  for (auto const& r : results) {
    std::string row = Form("%s,%s,%.3f,%.4f,%.4f,%.4f,%.4f",r.model.name.c_str(),gmnx::TPEGroupName(r.model.tpe).c_str(),
			   Q2,eps,in.R,r.val,r.stat);
    for (double const e : r.sys) row += Form(",%.4f",e);
    row += Form(",%.4f,%.4f,%.4f,%.4f,%.4f,%.4f",r.pCS,r.GEnSys,r.sysP2P,r.sysScale,r.statP2P,r.tot);
    std::cout << row << "\n";
  }
  std::cout << "\n";

  return 0;
}
