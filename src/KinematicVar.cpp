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
  void GetxyHCALexpect(TVector3 vertex, TVector3 pNhat, TVector3 HCAL_origin, 
		       std::vector<TVector3> HCAL_axes, std::vector<double> &xyHCALexpect) {
    /* This function calculates the expected vertical (x) and horizontal (y) positions
     of the recoil nucleon at the face of HCAL. 
     input:
     1. vertex       : vertex vector [in Hall CoS?? It must be but haven't confirmed], 
     2. pNhat        : projected q vector, 
     3. HCAL_origin  : HCAL origin vector [in Hall CoS], 
     4. HCAL_axes    : HCAL CoS axes [in Hall CoS]
     output:
     1. xyHCALexpect : expected x and y positions
    */
    // Intersection of a ray with a plane
    double sintersect = (HCAL_origin - vertex).Dot(HCAL_axes[2]) / (pNhat.Dot(HCAL_axes[2]));
    // ray from Hall origin onto the face of HCAL where the nucleon hit
    TVector3 HCAL_intersect = vertex + sintersect*pNhat; 

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
    /* Calculates Born cross-section ratio for MC (Using the same parametrizations used in SIMC/g4sbs) */
    // EMFF fits
    Kelly2004 kellyfit;
    Seamus20XX seamusfit;
    // EMFF extraction using same parametrization used in MC generators
    double GEp_kelly = kellyfit.GetFF(G_t::kGEp,Q2);
    double GMp_kelly = kellyfit.GetFF(G_t::kGMp,Q2);
    double GEn_seamus = seamusfit.GetFF(G_t::kGEn,Q2);
    double GMn_kelly = kellyfit.GetFF(G_t::kGMn,Q2);
    // calculate ratios
    double sigmaBorn_p_wo_Mott = kine::sigmaBorn_wo_Mott(etheta,Q2,GEp_kelly,GMp_kelly,"p");
    double sigmaBorn_n_wo_Mott = kine::sigmaBorn_wo_Mott(etheta,Q2,GEn_seamus,GMn_kelly,"n");
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
  double CalcGMnError(double etheta, double Q2, double GEp, double GMp, double GEn, double ratio, double ratioErr) {
    /* Performs GMn error propagation from a given ratio & error */
    double tau_p = kine::tau(Q2,"p");
    double tau_n = kine::tau(Q2,"n");
    double epsilon_p = kine::epsilon(etheta,Q2,"p");
    double epsilon_n = kine::epsilon(etheta,Q2,"n");
    double sigmaReduced_p = kine::sigmaReduced(tau_p,epsilon_p,GEp,GMp);
    // defining some terms for convenience
    double term1 = epsilon_n*(1.+tau_n) / (epsilon_p*(1.+tau_p));

    double numer = (term1*sigmaReduced_p/tau_n) * ratioErr;
    double denom = 2.0*pow((term1*sigmaReduced_p*ratio - epsilon_n*GEn*GEn)/tau_n , 0.5);
    return numer/denom;
  }  
  //--------------------------------------------
  void ExtractGMnWithError(double etheta,double Q2,double GEp,double GMp,double GEn,double ratio,double ratioErr,std::vector<double> &output) {
    /* Calulates GMn and the associated error from a given ratio and error */
    double GMn = kine::ExtractGMn(etheta,Q2,GEp,GMp,GEn,ratio);
    double GMnErr = kine::ExtractGMn(etheta,Q2,GEp,GMp,GEn,ratio,ratioErr);

    output = {GMn,GMnErr};
  }
  
} //::kine
