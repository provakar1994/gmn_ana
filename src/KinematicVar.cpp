#include "../include/KinematicVar.h"

namespace kine {

  //--------------------------------------------
  double M_N(std::string Ntype) {
    double temp = 0.;
    if (Ntype.compare("p") == 0) 
      temp = constant::Mp;
    else if (Ntype.compare("n") == 0) 
      temp = constant::Mn;
    else if (Ntype.compare("np") == 0) 
      temp = 0.5*(constant::Mn + constant::Mp);
    else
      std::cerr << "[KinematicVar::M_N] Enter a valid nucleon type! **!**" << std::endl;
    return temp;
  }
  //--------------------------------------------
  double pelas(double ebeam, double etheta, std::string Ntype) {
    /* Scattered e- momentum using elastic kinematics.
       Neglecting electron rest mass, hence E = p */
    return ebeam/(1. + (ebeam/kine::M_N(Ntype))*(1.0 - cos(etheta)));
  }
  //--------------------------------------------
  double pelas(SBSconfig sbsconf, std::string Ntype) {
    /* Scattered e- momentum using elastic kinematics.
       Neglecting electron rest mass, hence E = p */
    double ebeam = sbsconf.GetEbeam();
    double etheta = sbsconf.GetBBtheta_rad();
    return ebeam/(1. + (ebeam/kine::M_N(Ntype))*(1.0 - cos(etheta)));
  }
  //--------------------------------------------
  double thelas(double ebeam, double epprime, std::string Ntype) {
    /* e- scattering angle using elastic kinematics. */
    return acos((kine::M_N(Ntype)/ebeam)*(1.-(ebeam/epprime))+1.);
  }
  //--------------------------------------------
  double thelas(SBSconfig sbsconf, std::string Ntype) {
    /* Scattered e- momentum using elastic kinematics. */
    double ebeam = sbsconf.GetEbeam();
    double etheta = sbsconf.GetBBtheta_rad();
    double epprime = pelas(sbsconf,Ntype);
    return acos((kine::M_N(Ntype)/ebeam)*(1.-(ebeam/epprime))+1.);
  }
  //--------------------------------------------
  double etheta(TLorentzVector Peprime) {
    return acos(Peprime.Pz() / Peprime.E());
  }
  //--------------------------------------------
  double ephi(TLorentzVector Peprime) {
    return atan2(Peprime.Py(), Peprime.Px());
  }
  //--------------------------------------------
  double ethbend(double thtg, double phtg, double thfp, double phfp) {
    /* trajectory bend angle calculation */

    // track direction at target
    TVector3 enhat_tgt( thtg, phtg, 1.0 );
    enhat_tgt = enhat_tgt.Unit();

    // track direction at focal plane
    TVector3 enhat_fp( thfp, phfp, 1.0 );
    enhat_fp = enhat_fp.Unit();

    // rotating focal plane Csys by 10deg clockwise about Y to match
    // target csys
    double GEMpitch = 10.0*TMath::DegToRad();
    TVector3 GEMzaxis(-sin(GEMpitch),0,cos(GEMpitch));
    TVector3 GEMyaxis(0,1,0);
    TVector3 GEMxaxis = (GEMyaxis.Cross(GEMzaxis)).Unit();

    // calculate the transformed fp coordinates
    TVector3 enhat_fp_rot = enhat_fp.X()*GEMxaxis + enhat_fp.Y()*GEMyaxis + enhat_fp.Z()*GEMzaxis;

    // trajectory bend angle
    double thbend = acos( enhat_fp_rot.Dot( enhat_tgt ) );
    return thbend;
  }
  //--------------------------------------------
  void SetPN(std::string Ntype, TLorentzVector &PN) {
    PN.SetPxPyPzE(0., 0., 0., kine::M_N(Ntype));
  } 
  //--------------------------------------------
  double pN_expect(double nu, std::string Ntype) {
    return sqrt(pow(nu, 2.) + 2. * kine::M_N(Ntype) * nu);
  }
  //--------------------------------------------
  TVector3 qVect_unit(double Ntheta, double Nphi) {
    TVector3 pNhat(sin(Ntheta) * cos(Nphi), sin(Ntheta) * sin(Nphi), cos(Ntheta));
    return pNhat;
  }
  //--------------------------------------------
  void SetHCALaxes(double sbstheta_rad,                           // SBS angle in radian 
		   std::vector<TVector3> &HCAL_axes) {
    TVector3 HCAL_zaxis(sin(-sbstheta_rad),0,cos(-sbstheta_rad)); // Clock-wise rotation about Y axis
    TVector3 HCAL_xaxis(0,-1,0);                                  // -Y axis of Hall CoS = X axis of HCAL CoS
    TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
    HCAL_axes.push_back(HCAL_xaxis);
    HCAL_axes.push_back(HCAL_yaxis);
    HCAL_axes.push_back(HCAL_zaxis);
  }
  //--------------------------------------------
  void GetxyHCALexpect(TVector3 vertex, TVector3 pNhat, TVector3 HCAL_origin, std::vector<TVector3> HCAL_axes,
		       TVector3 &HCAL_intersect, std::vector<double> &xyHCALexpect) {
    /* This function calculates the expected vertical (x) and horizontal (y) positions
     of the recoil nucleon at the face of HCAL. 
     input:
     1. vertex       : vertex vector [in Hall CoS?? It must be but haven't confirmed], 
     2. pNhat        : projected q vector, 
     3. HCAL_origin  : HCAL origin vector [in Hall CoS], 
     4. HCAL_axes    : HCAL CoS axes [in Hall CoS]
     output:
     1. HCAL_intersect 
     2. xyHCALexpect : expected x and y positions
    */
    // Intersection of a ray with a plane
    double sintersect = (HCAL_origin - vertex).Dot(HCAL_axes[2]) / (pNhat.Dot(HCAL_axes[2]));
    // ray from Hall origin onto the face of HCAL where the nucleon hit
    HCAL_intersect = vertex + sintersect*pNhat; 

    double xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot(HCAL_axes[0]);
    double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot(HCAL_axes[1]);

    xyHCALexpect.push_back(xexpect_HCAL);
    xyHCALexpect.push_back(yexpect_HCAL);
  }
  //--------------------------------------------
  double Q2(double ebeam, double eeprime, double etheta) {
    return 2.0*ebeam*eeprime*(1.0-cos(etheta));
  }
  //--------------------------------------------
  double Q2(SBSconfig sbsconf, std::string Ntype) {
    double ebeam = sbsconf.GetEbeam();
    double etheta = sbsconf.GetBBtheta_rad();
    double eeprime = kine::pelas(ebeam,etheta,Ntype);
    return 2.0*ebeam*eeprime*(1.0-cos(etheta));
  }
  //--------------------------------------------
  double tau(double Q2, std::string Ntype) {
    return Q2 / (4.0*pow(kine::M_N(Ntype),2));
  }
  //--------------------------------------------
  double tau(SBSconfig sbsconf, std::string Ntype) {
    return kine::Q2(sbsconf,Ntype) / (4.0*pow(kine::M_N(Ntype),2));
  }
  //--------------------------------------------
  double epsilon(double etheta, double Q2, std::string Ntype) {
    double tau = kine::tau(Q2,Ntype);
    return pow(1. + 2.*(1.+tau)*pow(tan(0.5*etheta),2) , -1);
  }
  //--------------------------------------------
  double epsilon(SBSconfig sbsconf, std::string Ntype) {
    double etheta = sbsconf.GetBBtheta_rad();
    double tau = kine::tau(sbsconf,Ntype);
    return pow(1. + 2.*(1.+tau)*pow(tan(0.5*etheta),2) , -1);
  }
  //--------------------------------------------
  double etheta_from_epsilon(double epsilon, double Q2, std::string Ntype) {
    /* e- scattering angle that gives the virtual photon polarization epsilon at Q2 for
       elastic scattering off nucleon Ntype, i.e. the inverse of kine::epsilon(etheta,Q2,Ntype):
         tan^2(etheta/2) = (1/eps - 1) / (2*(1 + tau))
       Useful to get the effective angle of acceptance-averaged (Q2, eps).
       Valid for 0 < eps <= 1; returns NaN (w/ an error message) otherwise. */
    if (!(epsilon > 0. && epsilon <= 1.)) {
      std::cerr << "*** ERROR***, [kine::etheta_from_epsilon] epsilon = " << epsilon
		<< " is outside (0, 1]" << std::endl;
      return std::nan("");
    }
    double tau = kine::tau(Q2,Ntype);
    return 2.0*atan(sqrt((1./epsilon - 1.) / (2.*(1. + tau))));
  }
  //--------------------------------------------
  double epsilon_general(double etheta, double Q2, double nu) {
    /* Doesn't assume Q2=2*M_N*nu ie elastic kinematics */
    return pow(1. + 2.*(1.+nu*nu/Q2)*pow(tan(0.5*etheta),2) , -1);
  }
  //--------------------------------------------
  double W2_general(double ebeam, double eeprime, double etheta, std::string Ntype) {
    /* Calculating W2 from 4-momentum conservation (No elastic assumptions) */
    double Q2 = 2.0*ebeam*eeprime*(1.0-cos(etheta));
    return pow(kine::M_N(Ntype),2.) + 2.*kine::M_N(Ntype)*(ebeam-eeprime) - Q2;
  }
  //--------------------------------------------
  double W2(double ebeam, double eeprime, double Q2, std::string Ntype) {
    return pow(kine::M_N(Ntype),2.0) + 2.0*kine::M_N(Ntype)*(ebeam-eeprime) - Q2;
  }
  //--------------------------------------------
  double W(double ebeam, double eeprime, double Q2, std::string Ntype) {
    return std::max(0., sqrt(kine::W2(ebeam, eeprime, Q2, Ntype)));
  }


  /* ###########################################
     ## Functions to calculate cross-sections ##  
     ########################################### */
  //--------------------------------------------
  double sigmaRutherford(double ebeam, double etheta) { 
    /* Rutherford scattering in relativistic limit */
    return pow(constant::alpha,2) / (4.*pow(ebeam,2)*pow(sin(0.5*etheta),4));
  }
  //--------------------------------------------
  double sigmaMott(double ebeam, double eeprime, double etheta) {
    /* Mott cross-section */
    double rf_cs = kine::sigmaRutherford(ebeam,etheta);
    return rf_cs*(eeprime/ebeam)*pow(cos(0.5*etheta),2);
  }
  //--------------------------------------------
  double sigmaMott(SBSconfig sbsconf, std::string Ntype) {
    /* Mott cross-section */
    double ebeam = sbsconf.GetEbeam();
    double etheta = sbsconf.GetBBtheta_rad();
    double eeprime = kine::pelas(ebeam,etheta,Ntype);
    double rf_cs = kine::sigmaRutherford(ebeam,etheta);
    return rf_cs*(eeprime/ebeam)*pow(cos(0.5*etheta),2);
  }
  //--------------------------------------------
  double sigmaReduced(double tau, double epsilon, double GE, double GM) {
    /* Calculates redused cross-section */
    return  epsilon*GE*GE + tau*GM*GM;
  }
  //--------------------------------------------
  double sigmaReducedError(double tau, double epsilon, double GE, double GEErr, double GM, double GMErr) {
    /* Error on the reduced cross-section sigmaR = eps*GE^2 + tau*GM^2 from the FF errors:
         d(sigmaR)/dGE = 2*eps*GE,  d(sigmaR)/dGM = 2*tau*GM
         => d(sigmaR) = 2*sqrt( (eps*GE*dGE)^2 + (tau*GM*dGM)^2 )
       GE & GM are taken as uncorrelated. In an L/T separation they are anti-correlated, so
       this is a (conservative) overestimate. See DB/FF_fits/Ye2018/README.md for its size. */
    return 2.0*std::sqrt(pow(epsilon*GE*GEErr,2) + pow(tau*GM*GMErr,2));
  }
  //--------------------------------------------
  double sigmaBorn(double ebeam, double eeprime, double etheta, double GE, double GM, std::string Ntype) {
    /* Calculates Born cross-section using Rosenbluth formula */
    double sigmaMott = kine::sigmaMott(ebeam,eeprime,etheta);
    double Q2 = kine::Q2(ebeam,eeprime,etheta);
    double tau = kine::tau(Q2,Ntype);
    double epsilon = kine::epsilon(etheta,Q2,Ntype);
    double sigmaReduced = kine::sigmaReduced(tau,epsilon,GE,GM);
    return sigmaMott * sigmaReduced * pow(epsilon*(1.+tau),-1);
  }
  //--------------------------------------------
  double sigmaBorn(SBSconfig sbsconf, double GE, double GM, std::string Ntype) {
    /* Calculates Born cross-section using Rosenbluth formula */
    double ebeam = sbsconf.GetEbeam();
    double etheta = sbsconf.GetBBtheta_rad();
    double eeprime = kine::pelas(ebeam,etheta,Ntype);
    double sigmaMott = kine::sigmaMott(ebeam,eeprime,etheta);
    double Q2 = kine::Q2(ebeam,eeprime,etheta);
    double tau = kine::tau(Q2,Ntype);
    double epsilon = kine::epsilon(etheta,Q2,Ntype);
    double sigmaReduced = kine::sigmaReduced(tau,epsilon,GE,GM);
    return sigmaMott * sigmaReduced * pow(epsilon*(1.+tau),-1);
  }
  //--------------------------------------------
  double sigmaBorn_wo_Mott(double etheta, double Q2, double GE, double GM, std::string Ntype) {
    /* Born cross-section w/o the Mott CS term */
    double tau = kine::tau(Q2,Ntype);
    double epsilon = kine::epsilon(etheta,Q2,Ntype);
    double sigmaReduced = kine::sigmaReduced(tau,epsilon,GE,GM);
    return sigmaReduced * pow(epsilon*(1.+tau),-1);
  }
  //--------------------------------------------
  double sigmaBorn_ratio(double ebeam, double eeprime, double etheta, double GEp, double GMp, double GEn, double GMn) {
    /* Calculates Born cross-section ratio */
    double sigmBorn_p = kine::sigmaBorn(ebeam,eeprime,etheta,GEp,GMp,"p");
    double sigmBorn_n = kine::sigmaBorn(ebeam,eeprime,etheta,GEn,GMn,"n");
    return sigmBorn_n / sigmBorn_p;
  }
  //--------------------------------------------
  double sigmaBorn_ratio(SBSconfig sbsconf, double GEp, double GMp, double GEn, double GMn) {
    /* Calculates Born cross-section ratio */
    double sigmBorn_p = kine::sigmaBorn(sbsconf,GEp,GMp,"p");
    double sigmBorn_n = kine::sigmaBorn(sbsconf,GEn,GMn,"n");
    return sigmBorn_n / sigmBorn_p;
  }
  //--------------------------------------------
  double sigmaBorn_ratio(double etheta, double Q2, double GEp, double GMp, double GEn, double GMn) {
    /* Calculates Born cross-section ratio assuming the Mott CS is the same for both p and n */
    double sigmaBorn_p_wo_Mott = kine::sigmaBorn_wo_Mott(etheta,Q2,GEp,GMp,"p");
    double sigmaBorn_n_wo_Mott = kine::sigmaBorn_wo_Mott(etheta,Q2,GEn,GMn,"n");
    return sigmaBorn_n_wo_Mott / sigmaBorn_p_wo_Mott;
  }
  //--------------------------------------------
  double sigmaBorn_ratio_MC(double etheta, double Q2) {
    /* Calculates Born cross-section ratio for MC, using exactly the EMFFs of the SIMC
       generator (bigbite branch) that produced the GMn MC. See SIMCbigbite in EMFFFits.h */
    SIMCbigbite simcfit;
    double GEp_simc = simcfit.GetFF(G_t::kGEp,Q2);
    double GMp_simc = simcfit.GetFF(G_t::kGMp,Q2);
    double GEn_simc = simcfit.GetFF(G_t::kGEn,Q2);
    double GMn_simc = simcfit.GetFF(G_t::kGMn,Q2);
    // calculate ratios
    double sigmaBorn_p_wo_Mott = kine::sigmaBorn_wo_Mott(etheta,Q2,GEp_simc,GMp_simc,"p");
    double sigmaBorn_n_wo_Mott = kine::sigmaBorn_wo_Mott(etheta,Q2,GEn_simc,GMn_simc,"n");
    return sigmaBorn_n_wo_Mott / sigmaBorn_p_wo_Mott;
  }
  //--------------------------------------------
  double sigmaBorn_ratio_MC(SBSconfig sbsconf, std::string Ntype) {
    /* Calculates Born cross-section ratio for MC (Using the same parametrizations used in SIMC/g4sbs) */
    double Q2 = kine::Q2(sbsconf,Ntype);
    double etheta = sbsconf.GetBBtheta_rad();
    return kine::sigmaBorn_ratio_MC(etheta,Q2);
  }

  /* #######################################################
     ## Functions to extract GMn from Born CS ratio ##  
     ####################################################### */
  //--------------------------------------------
  double ExtractGMn(double ebeam, double eeprime, double etheta, double GEp, double GMp, double GEn, double ratio) {
    /* Calculates GMn from a given ratio */
    double Q2 = kine::Q2(ebeam,eeprime,etheta);
    double tau_p = kine::tau(Q2,"p");
    double tau_n = kine::tau(Q2,"n");
    double epsilon_p = kine::epsilon(etheta,Q2,"p");
    double epsilon_n = kine::epsilon(etheta,Q2,"n");
    double sigmaReduced_p = kine::sigmaReduced(tau_p,epsilon_p,GEp,GMp);
    // defining some terms for convenience
    double term1 = epsilon_n*(1.+tau_n) / (epsilon_p*(1.+tau_p));

    return - pow((term1*sigmaReduced_p*ratio - epsilon_n*GEn*GEn)/tau_n , 0.5);
  }
  //--------------------------------------------
  double ExtractGMn(SBSconfig sbsconf, double GEp, double GMp, double GEn, double ratio) {
    /* Calculates GMn from a given ratio */
    double ebeam = sbsconf.GetEbeam();
    double etheta = sbsconf.GetBBtheta_rad();
    double eeprime_p = kine::pelas(ebeam,etheta,"p");
    double eeprime_n = kine::pelas(ebeam,etheta,"n");
    double Q2_p = kine::Q2(ebeam,eeprime_p,etheta);
    double Q2_n = kine::Q2(ebeam,eeprime_n,etheta);
    double tau_p = kine::tau(Q2_p,"p");
    double tau_n = kine::tau(Q2_n,"n");
    double epsilon_p = kine::epsilon(etheta,Q2_p,"p");
    double epsilon_n = kine::epsilon(etheta,Q2_n,"n");
    double sigmaReduced_p = kine::sigmaReduced(tau_p,epsilon_p,GEp,GMp);
    // defining some terms for convenience
    double term1 = epsilon_n*(1.+tau_n) / (epsilon_p*(1.+tau_p));

    return - pow((term1*sigmaReduced_p*ratio - epsilon_n*GEn*GEn)/tau_n , 0.5);
  }
  //--------------------------------------------
  double ExtractGMn(double etheta, double Q2, double GEp, double GMp, double GEn, double ratio) {
    /* Calculates GMn from a given ratio */
    double tau_p = kine::tau(Q2,"p");
    double tau_n = kine::tau(Q2,"n");
    double epsilon_p = kine::epsilon(etheta,Q2,"p");
    double epsilon_n = kine::epsilon(etheta,Q2,"n");
    double sigmaReduced_p = kine::sigmaReduced(tau_p,epsilon_p,GEp,GMp);
    // defining some terms for convenience
    double term1 = epsilon_n*(1.+tau_n) / (epsilon_p*(1.+tau_p));

    return - pow((term1*sigmaReduced_p*ratio - epsilon_n*GEn*GEn)/tau_n , 0.5);
  }
  //--------------------------------------------
  /* GMn error propagation (first order, one function per uncertain input)
     ----
     GMn^2 = F(R, sigmaR_p, GEn) = (term1*sigmaR_p*R - eps_n*GEn^2)/tau_n,
     with term1 = eps_n*(1+tau_n) / (eps_p*(1+tau_p)). Since F = GMn^2, for any input x:
       dGMn/dx = (dF/dx) / (2*GMn)
     Each function below returns |dGMn/dx|*dx (the sign only gives the direction of the
     shift and drops out when uncorrelated contributions are added in quadrature).
     The derivation is documented in scripts/gmn_extraction/README.md. */
  namespace {
    // Kinematic factors and |GMn| shared by the GMn error functions at fixed (etheta, Q2)
    struct GMnTerms {
      double tau_n, epsilon_n, term1, sigmaReduced_p, absGMn;
    };
    GMnTerms CalcGMnTerms(double etheta, double Q2, double GEp, double GMp, double GEn, double ratio) {
      GMnTerms t;
      double const tau_p = kine::tau(Q2,"p");
      double const epsilon_p = kine::epsilon(etheta,Q2,"p");
      t.tau_n = kine::tau(Q2,"n");
      t.epsilon_n = kine::epsilon(etheta,Q2,"n");
      t.term1 = t.epsilon_n*(1.+t.tau_n) / (epsilon_p*(1.+tau_p));
      t.sigmaReduced_p = kine::sigmaReduced(tau_p,epsilon_p,GEp,GMp);
      t.absGMn = std::abs(kine::ExtractGMn(etheta,Q2,GEp,GMp,GEn,ratio));
      return t;
    }
  }
  //--------------------------------------------
  double CalcGMnError(double etheta, double Q2, double GEp, double GMp, double GEn, double ratio, double ratioErr) {
    /* GMn error from the error on the n/p Born CS ratio R (exact, incl. GEn term):
       dF/dR = term1*sigmaR_p/tau_n  =>  dGMn = term1*sigmaR_p*dR / (2*tau_n*|GMn|) */
    GMnTerms const t = CalcGMnTerms(etheta,Q2,GEp,GMp,GEn,ratio);
    return t.term1*t.sigmaReduced_p*std::abs(ratioErr) / (2.0*t.tau_n*t.absGMn);
  }
  //--------------------------------------------
  double CalcGMnError_sigmaRp(double etheta, double Q2, double GEp, double GMp, double GEn, double ratio, double sigmaRpErr) {
    /* GMn error from the error on the reduced proton CS sigmaR_p (e.g. from kine::sigmaReducedError):
       dF/d(sigmaR_p) = term1*R/tau_n  =>  dGMn = term1*R*d(sigmaR_p) / (2*tau_n*|GMn|)
       GMn depends on GEp & GMp only through sigmaR_p, so this covers both proton FF errors. */
    GMnTerms const t = CalcGMnTerms(etheta,Q2,GEp,GMp,GEn,ratio);
    return t.term1*std::abs(ratio)*std::abs(sigmaRpErr) / (2.0*t.tau_n*t.absGMn);
  }
  //--------------------------------------------
  double CalcGMnError_GEn(double etheta, double Q2, double GEp, double GMp, double GEn, double ratio, double GEnErr) {
    /* GMn error from the error on GEn:
       dF/dGEn = -2*eps_n*GEn/tau_n  =>  dGMn = eps_n*|GEn|*dGEn / (tau_n*|GMn|) */
    GMnTerms const t = CalcGMnTerms(etheta,Q2,GEp,GMp,GEn,ratio);
    return t.epsilon_n*std::abs(GEn)*std::abs(GEnErr) / (t.tau_n*t.absGMn);
  }
  //--------------------------------------------
  void ExtractGMnWithError(double etheta,double Q2,double GEp,double GMp,double GEn,double ratio,double ratioErr,std::vector<double> &output) {
    /* Calulates GMn and the associated error from a given ratio and error */
    double GMn = kine::ExtractGMn(etheta,Q2,GEp,GMp,GEn,ratio);
    double GMnErr = kine::CalcGMnError(etheta,Q2,GEp,GMp,GEn,ratio,ratioErr);

    output = {GMn,GMnErr};
  }
  
} //::kine
