/* GMn extraction from the n/p Born cross-section ratio R, for a set of EMFF models,
   with full error propagation.
   ----
   Model naming "A-B": the proton FFs (i.e. the reduced proton CS sigma_R,p) come from
   fit A, GEn from fit B. "noTPE" marks TPE-uncorrected (effective) proton FFs.
   Errors (all on GMn/(mu_n*GD)):
   - stat and every experimental systematic on R are propagated per model.
   - pCS (from the proton FF errors) and GEn (from the GEn error) are computed once, at
     the primary model, and the same values are used for all models.
   - Systematics whose name ends in "_sc" are scale errors; all others (incl. pCS & GEn)
     are point-to-point. Each group is added in quadrature separately.
   The error propagation formulas are derived in scripts/gmn_extraction/README.md.
*/
#ifndef GMN_EXTRACTION_H
#define GMN_EXTRACTION_H

#include <string>
#include <vector>

namespace gmnx {

  // How the proton FFs used for sigma_R,p treat two-photon exchange (TPE)
  enum class TPEGroup { kUncorrected, kCorrected, kMixed };
  std::string TPEGroupName(TPEGroup tpe);

  // One EMFF model combination (one row of the output)
  struct ModelDef {
    std::string name;   // e.g. "Ye noTPE-Ye"
    TPEGroup tpe;       // TPE treatment of the proton FFs
    double GEp, GMp;    // proton FFs (-> sigma_R,p) from fit A
    double GEn;         // GEn from fit B
  };

  // All model combinations at Q2 (GeV2), in output order. The primary is ModelDefs()[kPrimaryModel].
  std::vector<ModelDef> ModelDefs(double Q2);
  int const kPrimaryModel = 1; // "Ye noTPE-Ye"

  // Errors of the primary model FFs used for the pCS and GEn systematics (absolute)
  struct FFErrors { double dGEp, dGMp, dGEn; };
  FFErrors PrimaryFFErrors(double Q2);

  // One systematic error source on R
  struct ErrorSource {
    std::string name;
    double dR;      // absolute error on R
    bool isScale;   // scale error (name ends in "_sc"), otherwise point-to-point
  };
  // Builds an ErrorSource; isScale is set from the "_sc" name suffix
  ErrorSource MakeErrorSource(std::string const& name, double dR);

  // Measured n/p Born CS ratio and its errors at one kinematic point.
  // The kinematics are defined by the (acceptance-averaged) Q2 and eps. The Born CS
  // ratio depends on the kinematics only through them; eps_n follows from the same
  // e- angle, which is derived from eps (= eps_p) and stored in etheta.
  struct Input {
    double Q2;                      // GeV2
    double eps;                     // virtual photon polarization (proton, eps_p)
    double etheta;                  // derived: e- angle (rad) that gives eps_p = eps at Q2
    double R, dR_stat;              // n/p Born CS ratio and its stat error
    std::vector<ErrorSource> sys;   // experimental systematics on R (uncorrelated)
  };
  // Input from the measured n/p Born CS ratio R and its errors (absolute)
  Input MakeInput(double Q2, double eps, double R, double R_stat,
		  std::vector<std::string> const& sysNames, std::vector<double> const& R_sys);
  // Input from the data/MC ratio fit parameter R_fit: R = R_MC*R_fit (and the same for errors),
  // with R_MC the Born CS ratio of the MC generator (kine::sigmaBorn_ratio_MC)
  Input InputFromRfit(double Q2, double eps, double R_fit, double R_fit_stat,
		      std::vector<std::string> const& sysNames, std::vector<double> const& R_fit_sys);

  // GMn and its errors for one model. All errors are on GMn/(mu_n*GD).
  struct ModelResult {
    ModelDef model;
    double GMn;               // GMn
    double val;               // GMn/(mu_n*GD)
    double stat;              // from dR_stat
    std::vector<double> sys;  // from each experimental systematic (same order as Input::sys)
    double pCS;               // from the sigma_R,p error (primary model, same for all models)
    double GEnSys;            // from the GEn error (primary model, same for all models)
    double sysP2P;            // point-to-point sys: p2p experimental sources (+) pCS (+) GEnSys
    double sysScale;          // scale sys: "_sc" sources
    double statP2P;           // stat (+) sysP2P
    double tot;               // stat (+) sysP2P (+) sysScale
  };

  // GMn for all models of ModelDefs(), in the same order
  std::vector<ModelResult> ExtractAll(Input const& in);

  // sqrt(sum of squares)
  double QuadSum(std::vector<double> const& v);

} //::gmnx

#endif
