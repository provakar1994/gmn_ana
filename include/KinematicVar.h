/* This namespace holds functions to calculate interesting physics quantities for GMn-nTPE 
   Reaction: e + e' -> N + N' (fixed target)
   4-momentum conservation: Pe + Peprime -> PN + PNprime
   Please try to be consistent with the above naming convension.
   -----
   Created P. Datta <pdbforce@jlab.org> 05-27-2022
*/

#ifndef KINE_VAR_H
#define KINE_VAR_H

#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"

#include "EMFFFits.h"
#include "Constants.h"
#include "ExpConstants.h"

namespace kine{

  // choosing proper nucleon mass depending on its type
  double M_N(std::string Ntype);    // Options: "p", "n", & "np"

  // scattered e- p and angle calc. assuming elastic kinematics
  double pelas(double ebeam, double etheta, std::string Ntype); // GeV (Use for per event calc.) 
  double pelas(SBSconfig sbsconf, std::string Ntype);
  double thelas(double ebeam, double eeprime, std::string Ntype); // rad 
  double thelas(SBSconfig sbsconf, std::string Ntype);

  // utility functions to calc. scattering angle using tr.p* variables
  double etheta(TLorentzVector Peprime);  // Scattering angle (rad)
  double ephi(TLorentzVector Peprime);    // Angle of scattering plane (rad)

  // Constructing target nucleon 4-momentum (assuming at rest)
  void SetPN(std::string Ntype, TLorentzVector &PN);

  // expected recoil nucleon momentum
  double pN_expect(double nu,             // energy of the virtual photon (GeV)
		   std::string Ntype);    // type of nucleon in the reaction
		   
  // projected nucleon 3 vector using elastic kinematics (pNhat)
  TVector3 qVect_unit(double Ntheta,      // recoil nucleon theta (rad)
		      double Nphi);       // recoil nucleon phi (rad) 

  // Constructing HCAL co-ordinate system (CoS) in terms of Hall CoS
  void SetHCALaxes(double sbstheta_rad,           // SBS angle (rad)
		   std::vector<TVector3> &HCAL_axes);  // HCAL axes in order: X, Y, Z (Output)

  // Get the expected vertical (x) and horizontal (y) positions of the recoil 
  // nucleon at the face of HCAL.
  void GetxyHCALexpect(TVector3 vertex,                 // vertex vector [in Hall CoS]
		       TVector3 pNhat,                  // projected q vector
		       TVector3 HCAL_origin,            // HCAL origin vector [in Hall CoS]
		       std::vector<TVector3> HCAL_axes,      // HCAL CoS axes [in Hall CoS]
		       std::vector<double> &xyHCALexpect);   // expected x and y positions (Output)       

  double Q2(double ebeam, double eeprime, double etheta);                // GeV, GeV, rad (Use for per event calc.)
  double Q2(SBSconfig sbsconf, std::string Ntype);                       
  double tau(double Q2, std::string Ntype);                              // GeV2 (Use for per event calc.)
  double tau(SBSconfig sbsconf, std::string Ntype);
  double epsilon(double etheta, double Q2, std::string Ntype);           // rad, GeV2 (Use for per event calc.)
  double epsilon(SBSconfig sbsconf, std::string Ntype);
  double epsilon_general(double etheta, double Q2, double nu);
  double W2_general(double ebeam, double eeprime, double etheta, std::string Ntype); // GeV, GeV, rad
  double W2(double ebeam, double eeprime, double Q2, std::string Ntype); // GeV, GeV, GeV2
  double W(double ebeam, double eeprime, double Q2, std::string Ntype);  // GeV, GeV, GeV2

  /* #####################################
     ## Functions to get Cross-sections ##  
     ##################################### */
  double sigmaRutherford(double ebeam, double etheta);                   // GeV, rad
  double sigmaMott(double ebeam, double eeprime, double etheta);         // GeV, GeV, rad (Use for per event calc.)
  double sigmaMott(SBSconfig sbsconf, std::string Ntype);
  double sigmaReduced(double tau, double epsilon, double GE, double GM); 
  double sigmaBorn(double ebeam, double eeprime, double etheta, double GE, double GM, std::string Ntype);
  double sigmaBorn(SBSconfig sbsconf, double GE, double GM, std::string Ntype);
  double sigmaBorn_wo_Mott(double etheta, double Q2, double GE, double GM, std::string Ntype); // w/o sigmaMott term
  double sigmaBorn_ratio(double ebeam, double eeprime, double etheta, double GEp, double GMp, double GEn, double GMn);
  double sigmaBorn_ratio(SBSconfig sbsconf, double GEp, double GMp, double GEn, double GMn);
  // the following fn calculates the Born CS ratio assuming Mott CS are the same for p and n
  double sigmaBorn_ratio(double etheta, double Q2, double GEp, double GMp, double GEn, double GMn);
  double sigmaBorn_ratio_MC(double etheta, double Q2);
  double sigmaBorn_ratio_MC(SBSconfig sbsconf, std::string Ntype);
 
  /* #######################################################
     ## Functions to extract GMn from Born CS ratio ##  
     ####################################################### */
  double ExtractGMn(double ebeam, double eeprime, double etheta, double GEp, double GMp, double GEn, double ratio);
  double ExtractGMn(SBSconfig sbsconf, double GEp, double GMp, double GEn, double ratio);
  double ExtractGMn(double etheta, double Q2, double GEp, double GMp, double GEn, double ratio);
}

#endif
