#ifndef EXP_CONSTANTS_H
#define EXP_CONSTANTS_H

#include <iostream>
#include <unordered_map>

#include "TMath.h"
#include "TString.h"

namespace expconst {

  // Detector variables
  // BBCAL
  static const int shcol = 7;
  static const int shrow = 27;
  static const int pscol = 2;
  static const int psrow = 26;
  // HCAL 
  // 1. dimesions found from g4sbs: G4SBSHArmBuilder::MakeHCALV2
  static const int hcalcol = 12;
  static const int hcalrow = 24;
  static const double hcalblk_w = 0.1524;       //m, width of a HCAL block
  static const double hcalblk_h = 0.1524;       //m, height of a HCAL block
  static const double hcalblk_cTc_h = 0.15494;  //m, horizontal center-to-center dist.
  static const double hcalblk_cTc_v = 0.15875;  //m, vertical center-to-center dist.
  static const double hcalblk_gap_h = 0.00254;  //m, horiz. gap bet. two blocks
  static const double hcalblk_gap_v = 0.00635;  //m, vert. gap bet. two blocks
  // 2. x & y pos of blocks sitting at all 4 edges
  // pass 1 (assuming HCAL frame introduces ~ 45cm vertical offset)
  /* static const double xHCAL_t_DB = -2.190625;   //m, center of top row blocks (from SBS-replay/DB) */
  /* static const double xHCAL_b_DB = 1.460625;    //m, center of bottom row blocks (from SBS-replay/DB) */
  /* static const double yHCAL_r_DB = -0.85217;    //m, center of right most blocks (from SBS-replay/DB) */
  /* static const double yHCAL_l_DB = 0.85217;     //m, center of left most blocks (from SBS-replay/DB) */
  // pass 2 (assuming HCAL frame introduces 75cm vertical offset, confirmed by GEp CAD file)
  static const double xHCAL_t_DB = -2.575625;   //m, center of top row blocks (from SBS-replay/DB)
  static const double xHCAL_b_DB = 1.075625;    //m, center of bottom row blocks (from SBS-replay/DB)
  static const double yHCAL_r_DB = -0.85217;    //m, center of right most blocks (from SBS-replay/DB)
  static const double yHCAL_l_DB = 0.85217;     //m, center of left most blocks (from SBS-replay/DB)
  // --- above => data, below => simu
  // pass 1 (assuming HCAL frame introduces 45cm vertical offset)
  /* static const double xHCAL_t_DB_MC = -2.27563; //m, center of top row blocks (from SBS-replay/DB_MC) */
  /* static const double xHCAL_b_DB_MC = 1.37562;  //m, center of bottom row blocks (from SBS-replay/DB_MC) */
  /* static const double yHCAL_r_DB_MC = -0.85217; //m, center of right most blocks (from SBS-replay/DB_MC) */
  /* static const double yHCAL_l_DB_MC = 0.85217;  //m, center of left most blocks (from SBS-replay/DB_MC) */
  // pass 2 (assuming HCAL frame introduces 75cm vertical offset, confirmed by GEp CAD file)
  static const double xHCAL_t_DB_MC = -2.575625; //m, center of top row blocks (from SBS-replay/DB_MC)
  static const double xHCAL_b_DB_MC = 1.075625;  //m, center of bottom row blocks (from SBS-replay/DB_MC)
  static const double yHCAL_r_DB_MC = -0.85217; //m, center of right most blocks (from SBS-replay/DB_MC)
  static const double yHCAL_l_DB_MC = 0.85217;  //m, center of left most blocks (from SBS-replay/DB_MC)
  /* // 3. Offsets adjusted by looking at deltax and deltay distributions */
  /* static const double hcaloffset_v_data = -0.38; //m, vert. offset of HCAL origin w.r.t DB (data) */
  /* static const double hcaloffset_h_data = 0.15;  //m, horiz. offset of HCAL origin w.r.t DB (data) */
  /* static const double hcaloffset_v_simu = 0.0;   //m, vert. offset of HCAL origin w.r.t DB (simu) */
  /* static const double hcaloffset_h_simu = 0.0;   //m, horiz. offset of HCAL origin w.r.t DB (simu) */
  
  // Constant for the entire experiment
  /*
    NOTES:
    *. Target cell dimensions:
       - Cell diameter is still a guess. Using the value Andrew used.
       - Cell thickness is also a guess. Again using what Andrew used.
    *. Target density: Look at ~/OneDrive - University of Connecticut/workPD/proj/gmn-ana-resources/target
    *. Stopping power: 
       - Useful website by NIST ESTAR: https://physics.nist.gov/PhysRefData/Star/Text/ESTAR.html 
       - According to the website stopping power of hydrogen ranges from 5.74 to 5.82 MeV*cm2/g in the scattered
         e- energy of interest (2.1-3.6 GeV) for our analysis.
       - Above website don't have data on deuterium but I found a master's thesis from 1952 that shows that the 
         stopping power of deuterium is pretty close to hydrogen. https://open.library.ubc.ca/media/stream/pdf/831/1.0085416/1
  */
  // target
  static const double tgtlen = 15.0;  //cm 
  // LH2
  static const double lh2_TgtRho = 0.0725;      //g/cc, target density
  static const double lh2_CellThick = 0.02;     //cm, target cell thickness (Andrew's guess)
  static const double lh2_CellDiam  = 1.6*2.54; //cm, target cell diameter (Andrew's guess)
  static const double lh2_uWinThick = 0.0145;   //cm, upstream window thickness
  static const double lh2_dWinThick = 0.0158;   //cm, downstream window (tip) thickness
  static const double lh2_dWallThick = 0.0143;  //cm, downstream wall thickness
  //static const double lh2_dEdx = 0.005771;      //GeV*cm2/g, collisional stopping power (2.5GeV energy), NIST ESTAR 
  // LD2
  static const double ld2_TgtRho = 0.167;       //g/cc, target density
  static const double ld2_CellThick = 0.02;     //cm, target cell thickness (Andrew's guess)
  static const double ld2_CellDiam  = 1.6*2.54; //cm, target cell diameter (Andrew's guess)
  static const double ld2_uWinThick = 0.0125;   //cm, upstream window thickness
  static const double ld2_dWinThick = 0.0138;   //cm, downstream window (tip) thickness
  static const double ld2_dWallThick = 0.0136;  //cm, downstream wall thickness
  //static const double ld2_dEdx = 0.005771;      //GeV*cm2/g, collisional stopping power (2GeV energy), NIST ESTAR 
  // magnet
  static const double bbmaxcurr = 750;   //A, 100% BB magnet current
  static const double sbsmaxcurr = 2100; //A, 100% SBS magnet current
  static const double sbsdipolegap = (48.0*2.54) / 100.;  // ~1.22 m
  static const double sbsmaxfield = 3.1*atan(0.85 / (11.0-2.25-(sbsdipolegap/2.))) / (0.3*sbsdipolegap*0.7); // ~1.26 T (?)
  // Polyethylene (PE) shield near scattering chamber (Installed during SBS-11). See Utilities::GetElossInTgt for more info.
  static const double PE_Rho = 0.91;         //g/cc
  static const double PE_ShieldThick = 1.0;  //cm, 10 mm
  // Al shield near scattering chamber (Installed during SBS-11). See Utilities::GetElossInTgt for more info.
  static const double Al_Rho = 2.7;              //g/cc
  static const double Al_ShieldThick = 2.54/8.0; //cm, 1/8th inch
  //static const double Al_dEdx = 0.0021;          //GeV*cm2/g, collisional stopping power (1-4GeV energy), NIST ESTAR 

  // Following quantities vary with configuration
  double ebeam(int config);     //GeV
  double pcentral(int config);  //GeV
  double bbtheta(int config);   //deg
  double bbdist(int config);    //m
  double sbstheta(int config);  //deg
  double sbsdist(int config);   //m
  double hcaltheta(int config); //deg
  double hcaldist(int config);  //m
  double GetdEdxCollH(int const config,                 // SBS config
		      bool const is_before_scattering); // 1=>YES, 0=>NO (i.e. after scattering)
  double GetdEdxCollAl(int const config, bool const is_before_scattering);  //GeV*cm2/g
  double GetdEdxCollPE(int const config, bool const is_before_scattering);  //GeV*cm2/g
  /* double dEdx_coll_bs(int config, std::string material); //GeV*cm2/g */
  /* double dEdx_coll_as(int config, std::string material); //GeV*cm2/g */
}

// a class for SBS config
class SBSconfig {
 public:

  int    GetSBSconf()       const { return fSBSconf; }
  int    GetSBSmag()        const { return fSBSmag; }
  double GetEbeam()         const { return fEbeam; }
  double GetBBtheta()       const { return fBBtheta; }
  double GetBBtheta_rad()   const { return fBBtheta_rad; }
  double GetBBdist()        const { return fBBdist; }
  double GetSBStheta()      const { return fSBStheta; }
  double GetSBStheta_rad()  const { return fSBStheta_rad; }
  double GetSBSdist()       const { return fSBSdist; }
  double GetHCALtheta()     const { return fHCALtheta; }
  double GetHCALtheta_rad() const { return fHCALtheta_rad; }
  double GetHCALdist()      const { return fHCALdist; }

  // constructor
  SBSconfig(int conf, int sbsmag) {
    fSBSconf       = conf;
    fSBSmag        = sbsmag;
    fEbeam         = expconst::ebeam(conf);
    fBBtheta       = expconst::bbtheta(conf);
    fBBtheta_rad   = expconst::bbtheta(conf)*TMath::DegToRad();
    fBBdist        = expconst::bbdist(conf);
    fSBStheta      = expconst::sbstheta(conf);
    fSBStheta_rad  = expconst::sbstheta(conf)*TMath::DegToRad();
    fSBSdist       = expconst::sbsdist(conf);
    fHCALtheta     = expconst::hcaltheta(conf);
    fHCALtheta_rad = expconst::hcaltheta(conf)*TMath::DegToRad();
    fHCALdist      = expconst::hcaldist(conf);
  }

  // define an ostream operator to print to screen conveniently
  friend std::ostream& operator <<(std::ostream &out, const SBSconfig& sbsconf) {
    out  << " -------------------------- "                         << std::endl
	 << Form(" SBS Config: %d, "                   , sbsconf.fSBSconf)   << std::endl
	 << Form(" SBS Magnet Settings: %d (p), "      , sbsconf.fSBSmag)    << std::endl
    	 << Form(" Beam energy: %0.4f (GeV),"          , sbsconf.fEbeam)     << std::endl
    	 << Form(" BigBite angle: %0.1f (deg),"        , sbsconf.fBBtheta)   << std::endl
      	 << Form(" BigBite distance: %0.5f (m),"       , sbsconf.fBBdist)    << std::endl
    	 << Form(" Super BigBite angle: %0.1f (deg),"  , sbsconf.fSBStheta)  << std::endl
      	 << Form(" Super BigBite distance: %0.2f (m)," , sbsconf.fSBSdist)   << std::endl
	 << Form(" HCAL angle: %0.1f (deg),"           , sbsconf.fHCALtheta) << std::endl
    	 << Form(" HCAL distance: %0.1f (m)"           , sbsconf.fHCALdist)  << std::endl
	 << " -------------------------- "                        << std::endl << std::endl;
    return out;
  }

 private:
  int    fSBSconf;             // SBS configuration number
  int    fSBSmag;              // SBS magnet settings (%)
  double fEbeam;               // beam energy (better to get this from tree) (GeV)
  double fBBtheta;             // BigBite magnet angle (deg)
  double fBBtheta_rad;         // BigBite magnet angle (rad)
  double fBBdist;              // BigBite magnet distance from target (m)
  double fSBStheta;            // Super BigBite magnet angle (deg)
  double fSBStheta_rad;        // Super BigBite magnet angle (rad)
  double fSBSdist;             // Super BigBite magnet distance from target (m)
  double fHCALtheta;           // HCAL angle (deg)
  double fHCALtheta_rad;       // HCAL angle (rad)
  double fHCALdist;            // HCAL distance from target (m)
};

#endif
