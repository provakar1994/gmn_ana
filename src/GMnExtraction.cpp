#include "../include/GMnExtraction.h"

#include <cmath>
#include <iostream>

#include "../include/Constants.h"
#include "../include/EMFFFits.h"
#include "../include/KinematicVar.h"

namespace gmnx {

  namespace {
    // Value the EMFF fits return on failure (e.g. Q2 outside a lookup table)
    double const kFFErrorValue = -1000;
    // Name suffix that marks a scale (normalization) systematic
    std::string const kScaleSuffix = "_sc";

    // mu_n*GD: GMn is quoted as GMn/(mu_n*GD)
    double MunGD(double Q2) { return constant::mun*EMFFFits::GetGDip(Q2); }
  }

  //--------------------------------------------
  std::string TPEGroupName(TPEGroup tpe) {
    switch (tpe) {
    case TPEGroup::kUncorrected: return "uncorrected";
    case TPEGroup::kCorrected:   return "corrected";
    case TPEGroup::kMixed:       return "mixed";
    }
    return "unknown";
  }

  //--------------------------------------------
  std::vector<ModelDef> ModelDefs(double Q2) {
    /* The model combinations, "A-B" = proton FFs (sigma_R,p) from A, GEn from B.
       1. TPE-uncorrected proton FFs (consistent w/ the TPE-uncorrected n/p ratio, as no
          neutron TPE correction is applied). The primary uses Ye's fit, because it uses
          the same (older) radiative corrections as SIMC. Christy uses an updated RC.
       2. TPE-corrected (or mixed) proton FFs: for reference only.
       !! "Ye noTPE" is currently the Ye "noExtraTPE" fit, which is still SIFF-Blunden TPE
       corrected; a placeholder until a truly TPE-uncorrected Ye fit is available. */
    Ye2017 ye;
    Christy2022 christy;
    Arrington2007 arrington;
    Kelly2004 kelly;
    Galster1971 galster;

    double const GEp_yeNoTPE = ye.GetFFNoTPE(G_t::kGEp,Q2);
    double const GMp_yeNoTPE = ye.GetFFNoTPE(G_t::kGMp,Q2);
    double const GEn_ye = ye.GetFF(G_t::kGEn,Q2);

    std::vector<ModelDef> defs = {
      // TPE-uncorrected proton FFs
      {"Ch noTPE-Ye", TPEGroup::kUncorrected, christy.GetFFNoTPE(G_t::kGEp,Q2), christy.GetFFNoTPE(G_t::kGMp,Q2), GEn_ye},
      {"Ye noTPE-Ye", TPEGroup::kUncorrected, GEp_yeNoTPE, GMp_yeNoTPE, GEn_ye}, // primary
      {"Ar noTPE-Ye", TPEGroup::kUncorrected, arrington.GetFFNoTPE(G_t::kGEp,Q2), arrington.GetFFNoTPE(G_t::kGMp,Q2), GEn_ye},
      {"Ye noTPE-Ke", TPEGroup::kUncorrected, GEp_yeNoTPE, GMp_yeNoTPE, kelly.GetFF(G_t::kGEn,Q2)},
      {"Ye noTPE-Ga", TPEGroup::kUncorrected, GEp_yeNoTPE, GMp_yeNoTPE, galster.GetFF(G_t::kGEn,Q2)},
      // TPE-corrected (Ye, Arrington) or mixed (Kelly) proton FFs, for reference
      {"Ye-Ye", TPEGroup::kCorrected, ye.GetFF(G_t::kGEp,Q2), ye.GetFF(G_t::kGMp,Q2), GEn_ye},
      {"Ar-Ye", TPEGroup::kCorrected, arrington.GetFF(G_t::kGEp,Q2), arrington.GetFF(G_t::kGMp,Q2), GEn_ye},
      {"Ke-Ye", TPEGroup::kMixed, kelly.GetFF(G_t::kGEp,Q2), kelly.GetFF(G_t::kGMp,Q2), GEn_ye}
    };

    for (auto const& d : defs)
      if (d.GEp == kFFErrorValue || d.GMp == kFFErrorValue || d.GEn == kFFErrorValue)
	std::cerr << "*** ERROR***, [gmnx::ModelDefs] FF evaluation failed for " << d.name
		  << " at Q2 = " << Q2 << " GeV2" << std::endl;
    return defs;
  }

  //--------------------------------------------
  FFErrors PrimaryFFErrors(double Q2) {
    /* Errors of the primary model FFs ("Ye noTPE-Ye"):
       - GEp, GMp: Ye no-TPE fit error (currently the "noExtraTPE" fit: covariance only,
         w/o Ye's data-tension term; see DB/FF_fits/Ye2018/README.md)
       - GEn: Ye fit error (full tabulated error dG) */
    Ye2017 ye;
    FFErrors e;
    e.dGEp = ye.GetFFwErrNoTPE(G_t::kGEp,Q2)[1];
    e.dGMp = ye.GetFFwErrNoTPE(G_t::kGMp,Q2)[1];
    e.dGEn = ye.GetFFwErr(G_t::kGEn,Q2)[1];
    if (e.dGEp == kFFErrorValue || e.dGMp == kFFErrorValue || e.dGEn == kFFErrorValue)
      std::cerr << "*** ERROR***, [gmnx::PrimaryFFErrors] FF error evaluation failed at Q2 = "
		<< Q2 << " GeV2" << std::endl;
    return e;
  }

  //--------------------------------------------
  ErrorSource MakeErrorSource(std::string const& name, double dR) {
    bool const isScale = name.size() >= kScaleSuffix.size() &&
      name.compare(name.size() - kScaleSuffix.size(), kScaleSuffix.size(), kScaleSuffix) == 0;
    return {name, dR, isScale};
  }

  //--------------------------------------------
  Input InputFromRfit(double etheta, double Q2, double R_fit, double R_fit_stat,
		      std::vector<std::string> const& sysNames, std::vector<double> const& R_fit_sys) {
    /* Assuming the data/MC discrepancy in the n/p ratio is entirely due to the Born CS
       ratio, the measured ratio is R = R_MC*R_fit, and its errors scale the same way. */
    if (sysNames.size() != R_fit_sys.size())
      std::cerr << "*** ERROR***, [gmnx::InputFromRfit] # of sys names (" << sysNames.size()
		<< ") != # of sys errors (" << R_fit_sys.size() << ")" << std::endl;
    double const R_MC = kine::sigmaBorn_ratio_MC(etheta,Q2);
    Input in;
    in.etheta = etheta;
    in.Q2 = Q2;
    in.R = R_MC*R_fit;
    in.dR_stat = std::abs(R_MC*R_fit_stat);
    for (size_t i = 0; i < sysNames.size() && i < R_fit_sys.size(); ++i)
      in.sys.push_back(MakeErrorSource(sysNames[i], std::abs(R_MC*R_fit_sys[i])));
    return in;
  }

  //--------------------------------------------
  std::vector<ModelResult> ExtractAll(Input const& in) {
    double const absMunGD = std::abs(MunGD(in.Q2));
    std::vector<ModelDef> const defs = ModelDefs(in.Q2);

    // -- FF systematics, evaluated once at the primary model --
    // pCS: sigma_R,p error from the proton FF errors (GE & GM uncorrelated), propagated to GMn
    ModelDef const& prim = defs[kPrimaryModel];
    FFErrors const ffErr = PrimaryFFErrors(in.Q2);
    double const tau_p = kine::tau(in.Q2,"p");
    double const eps_p = kine::epsilon(in.etheta,in.Q2,"p");
    double const dSigmaRp = kine::sigmaReducedError(tau_p,eps_p,prim.GEp,ffErr.dGEp,prim.GMp,ffErr.dGMp);
    double const pCS = kine::CalcGMnError_sigmaRp(in.etheta,in.Q2,prim.GEp,prim.GMp,prim.GEn,in.R,dSigmaRp)/absMunGD;
    // GEn: GEn error propagated to GMn
    double const GEnSys = kine::CalcGMnError_GEn(in.etheta,in.Q2,prim.GEp,prim.GMp,prim.GEn,in.R,ffErr.dGEn)/absMunGD;

    // -- GMn & errors for every model --
    std::vector<ModelResult> results;
    for (auto const& m : defs) {
      ModelResult r;
      r.model = m;
      r.GMn = kine::ExtractGMn(in.etheta,in.Q2,m.GEp,m.GMp,m.GEn,in.R);
      r.val = r.GMn/MunGD(in.Q2);
      r.stat = kine::CalcGMnError(in.etheta,in.Q2,m.GEp,m.GMp,m.GEn,in.R,in.dR_stat)/absMunGD;
      r.pCS = pCS;
      r.GEnSys = GEnSys;
      // experimental systematics, sorted into point-to-point and scale groups
      std::vector<double> p2p = {pCS, GEnSys}, scale;
      for (auto const& s : in.sys) {
	double const err = kine::CalcGMnError(in.etheta,in.Q2,m.GEp,m.GMp,m.GEn,in.R,s.dR)/absMunGD;
	r.sys.push_back(err);
	(s.isScale ? scale : p2p).push_back(err);
      }
      r.sysP2P = QuadSum(p2p);
      r.sysScale = QuadSum(scale);
      r.statP2P = QuadSum({r.stat, r.sysP2P});
      r.tot = QuadSum({r.stat, r.sysP2P, r.sysScale});
      results.push_back(r);
    }
    return results;
  }

  //--------------------------------------------
  double QuadSum(std::vector<double> const& v) {
    double sum2 = 0;
    for (double const x : v) sum2 += x*x;
    return std::sqrt(sum2);
  }

} //::gmnx
