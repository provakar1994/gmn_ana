/*
  Unit test for the GMn extraction and its error propagation (kine:: in KinematicVar.h).
  ----
  Inputs: the SIMC generator FFs (SIMCbigbite, no lookup tables needed) at the GMn
  kinematics, so the "true" GMn is known. FF and R errors are arbitrary but realistic.
  Checks:
  1. Round trip: R = sigmaBorn_ratio(FFs), then ExtractGMn(R) gives back the input GMn.
  2. sigmaReducedError agrees with numerical derivatives of sigmaReduced (in quadrature).
  3. CalcGMnError (dR), CalcGMnError_GEn (dGEn) and CalcGMnError_sigmaRp (dGEp, dGMp via
     the chain rule) agree with numerical central-difference derivatives of ExtractGMn.
  4. The relative forms quoted in scripts/gmn_extraction/README.md hold, with
     S = eps_n*GEn^2 / (tau_n*GMn^2):
       dGMn_pCS/|GMn| = 0.5*(1+S) * dsigmaR_p/sigmaR_p
       dGMn_GEn/|GMn| = S * dGEn/|GEn|
*/
#include "EMFFFits.h"
#include "KinematicVar.h"
#include "TestUtils.h"
#include <cmath>
#include <functional>

namespace {
  double const kRoundOffTol = 1e-12; // analytic identities: only floating point round-off
  double const kNumDerivTol = 1e-6;  // analytic vs numerical derivative
  double const kStepFrac = 1e-4;     // central-difference step, as a fraction of the variable

  // Arbitrary relative errors used for the checks
  double const kRelErrGEp = 0.04, kRelErrGMp = 0.01, kRelErrGEn = 0.20, kRelErrR = 0.01;

  // |df/dx| * dx from a central difference with step h = kStepFrac*x
  double NumError(std::function<double(double)> const& f, double x, double dx) {
    double const h = kStepFrac*std::abs(x);
    return std::abs((f(x + h) - f(x - h)) / (2.*h)) * dx;
  }
}

int main() {
  // GMn kinematics {etheta (rad), Q2 (GeV2)} (scripts/lta/results.csv)
  double const kKine[][2] = {{0.623, 2.989}, {0.854, 4.488}, {0.465, 4.531},
			     {0.808, 7.464}, {0.694, 9.834}, {0.727, 13.465}};
  SIMCbigbite fit;

  for (auto const& k : kKine) {
    double const etheta = k[0], Q2 = k[1];
    std::string const at = " at Q2 = " + std::to_string(Q2);

    double const GEp = fit.GetFF(G_t::kGEp, Q2), GMp = fit.GetFF(G_t::kGMp, Q2);
    double const GEn = fit.GetFF(G_t::kGEn, Q2), GMn = fit.GetFF(G_t::kGMn, Q2);
    double const dGEp = kRelErrGEp*GEp, dGMp = kRelErrGMp*GMp, dGEn = kRelErrGEn*GEn;

    // -- 1. round trip --
    double const R = kine::sigmaBorn_ratio(etheta, Q2, GEp, GMp, GEn, GMn);
    double const dR = kRelErrR*R;
    double const GMnExtracted = kine::ExtractGMn(etheta, Q2, GEp, GMp, GEn, R);
    CHECK_REL(GMnExtracted, GMn, kRoundOffTol, "round trip GMn" + at);

    // -- 2. sigmaReducedError vs numerical derivatives --
    double const tau_p = kine::tau(Q2, "p"), eps_p = kine::epsilon(etheta, Q2, "p");
    double const sigRp = kine::sigmaReduced(tau_p, eps_p, GEp, GMp);
    double const dSigRp = kine::sigmaReducedError(tau_p, eps_p, GEp, dGEp, GMp, dGMp);
    double const dSigRp_GE = NumError([&](double x) { return kine::sigmaReduced(tau_p, eps_p, x, GMp); }, GEp, dGEp);
    double const dSigRp_GM = NumError([&](double x) { return kine::sigmaReduced(tau_p, eps_p, GEp, x); }, GMp, dGMp);
    CHECK_REL(dSigRp, std::hypot(dSigRp_GE, dSigRp_GM), kNumDerivTol, "sigmaReducedError" + at);

    // -- 3. GMn errors vs numerical derivatives of ExtractGMn --
    // from dR
    double const num_R = NumError([&](double x) { return kine::ExtractGMn(etheta, Q2, GEp, GMp, GEn, x); }, R, dR);
    CHECK_REL(kine::CalcGMnError(etheta, Q2, GEp, GMp, GEn, R, dR), num_R, kNumDerivTol, "CalcGMnError (dR)" + at);
    // from dGEn
    double const num_GEn = NumError([&](double x) { return kine::ExtractGMn(etheta, Q2, GEp, GMp, x, R); }, GEn, dGEn);
    CHECK_REL(kine::CalcGMnError_GEn(etheta, Q2, GEp, GMp, GEn, R, dGEn), num_GEn, kNumDerivTol, "CalcGMnError_GEn" + at);
    // from dGEp and dGMp separately, each mapped onto d(sigmaR_p) = 2*eps_p*GEp*dGEp or 2*tau_p*GMp*dGMp
    double const num_GEp = NumError([&](double x) { return kine::ExtractGMn(etheta, Q2, x, GMp, GEn, R); }, GEp, dGEp);
    double const num_GMp = NumError([&](double x) { return kine::ExtractGMn(etheta, Q2, GEp, x, GEn, R); }, GMp, dGMp);
    CHECK_REL(kine::CalcGMnError_sigmaRp(etheta, Q2, GEp, GMp, GEn, R, 2.*eps_p*GEp*dGEp), num_GEp, kNumDerivTol,
	      "CalcGMnError_sigmaRp (dGEp)" + at);
    CHECK_REL(kine::CalcGMnError_sigmaRp(etheta, Q2, GEp, GMp, GEn, R, 2.*tau_p*GMp*dGMp), num_GMp, kNumDerivTol,
	      "CalcGMnError_sigmaRp (dGMp)" + at);
    // both together: the pCS error equals the quadrature sum of the two
    double const dGMn_pCS = kine::CalcGMnError_sigmaRp(etheta, Q2, GEp, GMp, GEn, R, dSigRp);
    CHECK_REL(dGMn_pCS, std::hypot(num_GEp, num_GMp), kNumDerivTol, "CalcGMnError_sigmaRp (dGEp & dGMp)" + at);

    // -- 4. relative forms --
    double const S = kine::epsilon(etheta, Q2, "n")*GEn*GEn / (kine::tau(Q2, "n")*GMn*GMn);
    CHECK_REL(dGMn_pCS/std::abs(GMn), 0.5*(1. + S)*dSigRp/sigRp, kRoundOffTol*1e2, "relative form, pCS" + at);
    CHECK_REL(kine::CalcGMnError_GEn(etheta, Q2, GEp, GMp, GEn, R, dGEn)/std::abs(GMn), S*dGEn/std::abs(GEn),
	      kRoundOffTol*1e2, "relative form, GEn" + at);
  }

  return test::Summary("test_gmn_extraction");
}
