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
  double pcentral(double ebeam, double etheta, std::string Ntype) {
    return ebeam/(1. + (ebeam/kine::M_N(Ntype))*(1.0 - cos(etheta)));
  }
  //--------------------------------------------
  double pcentral(SBSconfig sbsconf, std::string Ntype) {
    double ebeam = sbsconf.GetEbeam();
    double etheta = sbsconf.GetBBtheta_rad();
    return ebeam/(1. + (ebeam/kine::M_N(Ntype))*(1.0 - cos(etheta)));
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
    double eeprime = kine::pcentral(ebeam,etheta,Ntype);
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
  double W2(double ebeam, double eeprime, double Q2, std::string Ntype) {
    return pow(kine::M_N(Ntype),2.0) + 2.0*kine::M_N(Ntype)*(ebeam-eeprime) - Q2;
  }
  //--------------------------------------------
  double W(double ebeam, double eeprime, double Q2, std::string Ntype) {
    return std::max(0., sqrt(kine::W2(ebeam, eeprime, Q2, Ntype)));
  }
  //--------------------------------------------
  double Luminosity(double ibeam, std::string targetType) {
    double lumi = 0.;
    if (targetType.compare("LH2") == 0)
      lumi = ((ibeam/constant::qe)*expconst::tgtlen*expconst::lh2_TgtRho*(constant::N_A/constant::H2_Amass));
    else if (targetType.compare("LD2") == 0)
      lumi = ((ibeam/constant::qe)*expconst::tgtlen*expconst::ld2_TgtRho*(constant::N_A/constant::D2_Amass));
    else
      std::cerr << "[KinematicVar::Luminosity] Enter a valid target type! **!**" << std::endl;
    return lumi;
  }

} //::kine
