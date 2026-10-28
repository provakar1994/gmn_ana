/*
  Unit test: the MC Born cross-section ratio R_MC must use exactly the EMFFs of the
  SIMC generator that produced the GMn MC.
  ----
  Reference: SIMC bigbite branch, physics_proton.f
  https://github.com/MarkKJones/simc_gfortran/blob/bigbite/physics_proton.f
  The two SIMC subroutines are copied below line by line (fofa_best_fit_p: Kelly GEp/GMp;
  fofa_best_fit_n: Kelly GMn and Riordan/Seamus GEn). Note that SIMC uses the proton mass
  in tau for both nucleons and the rounded magnetic moments GMp(0) = 2.79, GMn(0) = -1.913.
  Checks:
  1. SIMCbigbite::GetFF reproduces the SIMC FFs for 0.5 <= Q2 <= 15 GeV2.
  2. kine::sigmaBorn_ratio_MC equals kine::sigmaBorn_ratio evaluated with the SIMC FFs
     at the GMn kinematics.
*/
#include "EMFFFits.h"
#include "KinematicVar.h"
#include "TestUtils.h"
#include <cmath>

namespace simc {

  double const kMproton = 0.938272081; // GeV, as in physics_proton.f
  double const kLambda2 = 0.71;        // GeV2, dipole mass parameter

  // fofa_best_fit_p: proton FFs (Kelly, PRC 70, 068202)
  void ProtonFFs(double Q2, double& GE, double& GM) {
    double const tau = Q2/(4.0*kMproton*kMproton);
    double const tau2 = tau*tau, tau3 = tau*tau*tau;
    double const GEdenom = 1. + 10.98*tau + 12.82*tau2 + 21.97*tau3;
    GE = (1. - 0.24*tau)/GEdenom;
    double const GMdenom = 1. + 10.97*tau + 18.86*tau2 + 6.55*tau3;
    GM = (2.79*(1. + 0.12*tau))/GMdenom;
  }

  // fofa_best_fit_n: neutron FFs (GMn: Kelly; GEn: Riordan, PRL 105, 262302)
  void NeutronFFs(double Q2, double& GE, double& GM) {
    double const GD = 1./((1. + Q2/kLambda2)*(1. + Q2/kLambda2));
    double const tau = Q2/(4.0*kMproton*kMproton);
    double const tau2 = tau*tau, tau3 = tau*tau*tau;
    double const GEdenom = 1. + 5.222*tau + 0.04*tau2 + 11.438*tau3;
    GE = ((1.52*tau + 2.629*tau2 + 3.055*tau3)*GD)/GEdenom;
    double const GMdenom = 1. + 14.72*tau + 24.2*tau2 + 84.1*tau3;
    GM = (-1.913*(1. + 2.33*tau))/GMdenom;
  }

} //::simc

int main() {
  double const kTol = 1e-12; // same formulas, so only floating point round-off is allowed

  // -- 1. FF values on a Q2 grid covering the GMn kinematics --
  SIMCbigbite simcfit;
  for (double Q2 = 0.5; Q2 <= 15.0; Q2 += 0.5) {
    double GEp, GMp, GEn, GMn;
    simc::ProtonFFs(Q2, GEp, GMp);
    simc::NeutronFFs(Q2, GEn, GMn);
    std::string const at = " at Q2 = " + std::to_string(Q2);
    CHECK_REL(simcfit.GetFF(G_t::kGEp, Q2), GEp, kTol, "SIMC GEp" + at);
    CHECK_REL(simcfit.GetFF(G_t::kGMp, Q2), GMp, kTol, "SIMC GMp" + at);
    CHECK_REL(simcfit.GetFF(G_t::kGEn, Q2), GEn, kTol, "SIMC GEn" + at);
    CHECK_REL(simcfit.GetFF(G_t::kGMn, Q2), GMn, kTol, "SIMC GMn" + at);
  }

  // -- 2. R_MC at the GMn kinematics {etheta (rad), Q2 (GeV2)} (scripts/lta/results.csv) --
  double const kKine[][2] = {{0.623, 2.989}, {0.854, 4.488}, {0.465, 4.531},
			     {0.808, 7.464}, {0.694, 9.834}, {0.727, 13.465}};
  for (auto const& k : kKine) {
    double const etheta = k[0], Q2 = k[1];
    double GEp, GMp, GEn, GMn;
    simc::ProtonFFs(Q2, GEp, GMp);
    simc::NeutronFFs(Q2, GEn, GMn);
    double const R_simc = kine::sigmaBorn_ratio(etheta, Q2, GEp, GMp, GEn, GMn);
    CHECK_REL(kine::sigmaBorn_ratio_MC(etheta, Q2), R_simc, kTol, "R_MC at Q2 = " + std::to_string(Q2));
  }

  return test::Summary("test_simc_rmc");
}
