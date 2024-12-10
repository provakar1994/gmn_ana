/* 
   This macro will perform QE analysis for GMn using MC data. 
   E.g. Config. File: sbs14-sbs70p-simu/conf_qelas_ana_simu.json
   * A brief description of all the config. file parameters can be
   found at the bottom of this script. 
   -----
   P. Datta  Created  11-05-2022 
*/

// TO-DO
// 1. Energy loss calculation for g4sbs generator

// ***
// NOTE: 08/18/24
// T_W2 is W2recon + W2_offset but T_W has no offset correction

#include <vector>
#include <iostream>
#include <unordered_map>

#include "TCut.h"
#include "TH1F.h"
#include "TFile.h"
#include "TLatex.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TEllipse.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"

#include "../../include/gmn_ana.h"
#include "../../dflay/src/JSONManager.cxx"

/* this script will only analyze LD2 data */
static const std::string target = "LD2";

int qelas_ana_simu (const char *configfilename, 
		    int model=2, //Analysis model
		    std::string generator="simc",
		    std::string process="deeN")
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // clock to keep macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);
  std::string key = generator + "_" + process + "_model" + std::to_string(model);

  // setting verbosity
  int verbose = jmgr->GetValueFromSubKey<int>(key,"verbose");
  int verbosefn = jmgr->GetValueFromSubKey<int>(key,"verbose_function");

  // seting up the desired SBS configuration
  int conf = jmgr->GetValueFromSubKey<int>(key,"SBS_config");
  int sbsmag = jmgr->GetValueFromSubKey<int>(key,"SBS_magnet_percent");
  double sbsscalefield = jmgr->GetValueFromSubKey<double>(key,"SBS_scale_field");
  SBSconfig sbsconf(conf, sbsmag);
  std::cout << sbsconf;

  // reading job summary and parsing ROOT trees
  std::string rfd = jmgr->GetValueFromSubKey_str(key,"rootfile_dir");
  std::string prefix = jmgr->GetValueFromSubKey_str(key,"prefix_to_filebase");
  int njobs = jmgr->GetValueFromSubKey<int>(key,"Njobs_to_ana"); // # MC jobs to analyze
  std::vector<SimuJob> sjobs; util_pd::ReadSimuJobSummary(rfd,prefix,conf,sbsmag,generator,process,njobs,verbosefn,sjobs);
  TChain *C = new TChain("T"); util_pd::LoadSimuROOTTree(sjobs,verbosefn,C);

  // Choosing the model of calculation
  // model 0 => uses reconstructed p as independent variable
  // model 1 => uses reconstructed angles as independent variable
  // model 2 => uses 4-vector calculation
  if (model == 0) std::cout << "Using model 0 [recon. p as indep. var.] for analysis.." << std::endl;
  else if (model == 1) std::cout << "Using model 1 [recon. angle as indep. var.] for analysis.." << std::endl;
  else if (model == 2) std::cout << "Using model 2 [4-vector calculation] for analysis.." << std::endl;
  else { std::cerr << "Enter a valid model number! **!**" << std::endl; throw; }

  // choosing nucleon type 
  std::string Ntype = jmgr->GetValueFromSubKey_str(key,"Ntype");

  // setting up global cuts
  std::string gcut = jmgr->GetValueFromSubKey_str(key,"global_cut");
  std::vector<std::string> gCutList; util_pd::SplitString('&',gcut,gCutList);
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut",(TCut)gcut.c_str(),C);

  // setting up ROOT tree branch addresses ---------------------------------------
  int maxNtr = jmgr->GetValueFromSubKey<int>(key,"max_N_tracks");
  C->SetBranchStatus("*",0);
  // bbsh clus var
  double eSH, xSH, ySH, rblkSH, cblkSH, idblkSH, atimeSH;
  std::vector<std::string> bbshclvar = {"e","x","y","rowblk","colblk","idblk","atimeblk"};
  std::vector<void*> bbshclvar_mem = {&eSH,&xSH,&ySH,&rblkSH,&cblkSH,&idblkSH,&atimeSH};
  setrootvar::setbranch(C, "bb.sh", bbshclvar, bbshclvar_mem);

  // bbps clus var
  double ePS, rblkPS, cblkPS, idblkPS, atimePS;
  std::vector<std::string> bbpsclvar = {"e","rowblk","colblk","idblk","atimeblk"};
  std::vector<void*> bbpsclvar_mem = {&ePS,&rblkPS,&cblkPS,&idblkPS,&atimePS};
  setrootvar::setbranch(C, "bb.ps", bbpsclvar, bbpsclvar_mem);
 
  // hcal clus var
  double eHCAL, xHCAL, yHCAL, rblkHCAL, cblkHCAL, idblkHCAL, atimeHCAL, tdcHCAL[maxNtr];
  std::vector<std::string> hcalclvar = {"e","x","y","rowblk","colblk","idblk","atimeblk","clus_blk.tdctime"};
  std::vector<void*> hcalclvar_mem = {&eHCAL,&xHCAL,&yHCAL,&rblkHCAL,&cblkHCAL,&idblkHCAL,&atimeHCAL,&tdcHCAL};
  setrootvar::setbranch(C, "sbs.hcal", hcalclvar, hcalclvar_mem);

  // track var
  double ntrack,p[maxNtr],px[maxNtr],py[maxNtr],pz[maxNtr],xTr[maxNtr],yTr[maxNtr],thTr[maxNtr],phTr[maxNtr];
  double vx[maxNtr],vy[maxNtr],vz[maxNtr];
  double xtgt[maxNtr],ytgt[maxNtr],thtgt[maxNtr],phtgt[maxNtr],xfp[maxNtr],yfp[maxNtr],thfp[maxNtr],phfp[maxNtr];
  std::vector<std::string> trvar = {"n","p","px","py","pz","x","y","th","ph","vx","vy","vz","tg_x","tg_y","tg_th","tg_ph",
				    "r_x","r_y","r_th","r_ph"};
  std::vector<void*> trvar_mem = {&ntrack,&p,&px,&py,&pz,&xTr,&yTr,&thTr,&phTr,&vx,&vy,&vz,&xtgt,&ytgt,&thtgt,&phtgt,
				  &xfp,&yfp,&thfp,&phfp};
  setrootvar::setbranch(C,"bb.tr",trvar,trvar_mem);

  // GEM variables
  double nhitsGEM[maxNtr], ngoodhitsGEM[maxNtr], trchi2ndf[maxNtr];
  std::vector<std::string> gemvar = {"nhits","ngoodhits","chi2ndf"};
  std::vector<void*> gemvar_mem = {&nhitsGEM,&ngoodhitsGEM,&trchi2ndf};
  setrootvar::setbranch(C,"bb.gem.track",gemvar,gemvar_mem);

  //hcal variables
  bool gen_ML_data = 0; //= jmgr->GetValueFromSubKey<int>(key,"gen_ML_data"); // Generate ML data?
  int aHCALndata;
  double aHCAL[maxNtr],apHCAL[maxNtr],acHCAL[maxNtr],ampHCAL[maxNtr],amppHCAL[maxNtr],elHCAL[maxNtr];
  std::vector<std::string> hcalvar = {"a","a","a_p","a_c","a_amp","a_amp_p","adcelemID"};
  std::vector<void*> hcalvar_mem = {&aHCAL,&aHCALndata,&apHCAL,&acHCAL,&ampHCAL,&amppHCAL,&elHCAL};
  if (gen_ML_data) setrootvar::setbranch(C,"sbs.hcal",hcalvar,hcalvar_mem,1);

  //MC variables
  double mc_sigma;  // mc_sigma => Cross-section weight
  double mc_ebeam, mc_veE, mc_vetheta; // MC truth info for electrons
  double mc_fnucl, mc_np, mc_npx, mc_npy, mc_npz; // MC truth info for outgoing nucleon
  double mc_vx, mc_vy, mc_vz; // true vertex info
  std::vector<std::string> mc = {"mc_sigma","mc_fnucl","mc_np"};   // Default: g4sbs gen.
  std::vector<void*> mc_mem = {&mc_sigma,&mc_fnucl,&mc_np}; 
  if (generator.compare("simc")==0) {
    mc = {"simc_Weight","simc_fnucl","simc_Ebeam","simc_veE","simc_vetheta","simc_p_n","simc_px_n","simc_py_n","simc_pz_n","simc_vx","simc_vy","simc_vz"};  
    mc_mem = {&mc_sigma,&mc_fnucl,&mc_ebeam,&mc_veE,&mc_vetheta,&mc_np,&mc_npx,&mc_npy,&mc_npz,&mc_vx,&mc_vy,&mc_vz}; 
  }
  setrootvar::setbranch(C,"MC",mc,mc_mem);

  // turning on the remaining branches we use for the globalcut
  C->SetBranchStatus("bb.gem.track.nhits", 1);
  C->SetBranchStatus("bb.etot_over_p", 1);

  // Reading in HCAL efficiency map
  double avg_effi = 0.9485; //error weighted average of x and y efficiency
  TFile *fEffi = util_pd::ReadRootFile("~/gmn_ana/scripts/pdout/pDE/0p77zoff_pDE_data_sbs8_sbs-1p_model1_pass2.root");
  TH1F *h2_effi_map = (TH1F*)fEffi->Get("h2_effi_map");  

  // defining the outputfile
  std::string filebase = jmgr->GetValueFromSubKey_str(key,"outfile_prefix");
  filebase = filebase.empty() ? "" : filebase + "_";
  if (process.compare("inel")==0) filebase = filebase + "inel_";
  TString outFile = Form("siout/%sqelas_ana_%s_sbs%d_sbs%dp_model%d.root",filebase.c_str(),generator.c_str(),conf,sbsmag,model);
  TFile *fout = new TFile(outFile.Data(),"RECREATE");

  // defining histograms
  TH1F *h_W = util_pd::TH1FhW("h_W");
  TH1F *h_W_cut = util_pd::TH1FhW("h_W_cut");
  TH1F *h_W_acut = util_pd::TH1FhW("h_W_acut");
  TH1F *h_W2_cut = new TH1F("h_W2_cut","",200,-1,4);
  TH1F *h_W2_cut_noOff = new TH1F("h_W2_cut_noOff","",200,-1,4);  
  TH1F *h_dpel = new TH1F("h_dpel",";p/p_{elastic}(#theta)-1;",100,-0.3,0.3);
  
  TH1F *h_Q2 = util_pd::TH1FhQ2("h_Q2", conf);
  std::vector<double> hdx_lim; jmgr->GetVectorFromSubKey<double>(key,"h_dxHCAL_lims", hdx_lim);
  std::vector<double> hdy_lim; jmgr->GetVectorFromSubKey<double>(key,"h_dyHCAL_lims", hdy_lim);
  TH1F *h_dxHCAL = new TH1F("h_dxHCAL","W & fiducial cuts;x_{HCAL}^{obs} - x_{HCAL}^{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);
  TH1F *h_dxHCAL_nfc = new TH1F("h_dxHCAL_nfc","W cut;x_{HCAL}^{obs} - x_{HCAL}^{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);
  TH1F *h_dyHCAL = new TH1F("h_dyHCAL","W & fiducial cuts;y_{HCAL}^{obs} - y_{HCAL}^{exp} (m);",int(hdy_lim[0]),hdy_lim[1],hdy_lim[2]);
  TH1F *h_dyHCAL_nfc = new TH1F("h_dyHCAL_nfc","W cut;y_{HCAL}^{obs} - y_{HCAL}^{exp} (m);",int(hdy_lim[0]),hdy_lim[1],hdy_lim[2]);
  TH1F *h_dxHCAL_n = new TH1F("h_dxHCAL_n","mc_fnucl = n;x_{HCAL}^{obs} - x_{HCAL}^{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);
  TH1F *h_dxHCAL_p = new TH1F("h_dxHCAL_p","mc_fnucl = p;x_{HCAL}^{obs} - x_{HCAL}^{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);

  TH1F *h_dx_w_p_def = new TH1F("h_dx_w_p_def",";dx+p_def (m);",200,-1,1);
    
  TH2F *h2_rcHCAL = util_pd::TH2FHCALface_rc("h2_rcHCAL");
  TH2F *h2_dxdyHCAL = util_pd::TH2FdxdyHCAL("h2_dxdyHCAL");

  // Defining interesting ROOT tree branches 
  TTree *Tout = new TTree("Tout", "");
  Tout->SetMaxTreeSize(4000000000LL);
  //cuts
  bool WCut;              Tout->Branch("WCut", &WCut, "WCut/O");
  bool bbfiduCut;         Tout->Branch("bbfiduCut", &bbfiduCut, "bbfiduCut/O");
  bool pCut;              Tout->Branch("pCut", &pCut, "pCut/O");
  bool nCut;              Tout->Branch("nCut", &nCut, "nCut/O");
  double pdx_nS;          Tout->Branch("pdx_nS", &pdx_nS, "pdx_nS/D"); //# sigma away from p dx peak
  double ndx_nS;          Tout->Branch("ndx_nS", &ndx_nS, "ndx_nS/D"); //# sigma away from n dx peak
  double dy_nS;           Tout->Branch("dy_nS", &dy_nS, "dy_nS/D"); //# sigma away from dy peak 
  // --
  bool SMCut;             Tout->Branch("SMCut", &SMCut, "SMCut/B");
  bool ARCut;             Tout->Branch("ARCut", &ARCut, "ARCut/B");
  bool fiduCut;           Tout->Branch("fiduCut", &fiduCut, "fiduCut/O");
  //MC related
  double weight;          Tout->Branch("weight", &weight, "weight/D");
  double weight_effic;    Tout->Branch("weight_effic", &weight_effic, "weight_effic/D");
  double weight_norm;     Tout->Branch("weight_norm", &weight_norm, "weight_norm/D");  
  double weight_norm_effic; Tout->Branch("weight_norm_effic", &weight_norm_effic, "weight_norm_effic/D");  
  int T_mc_fnucl;         Tout->Branch("mc_fnucl", &T_mc_fnucl, "mc_fnucl/I");
  //
  double T_ebeam;         Tout->Branch("ebeam", &T_ebeam, "ebeam/D");
  //kine
  // -- vertex (only for SIMC) ---
  double T_veE;           if (generator.compare("simc")==0) Tout->Branch("veE", &T_veE, "veE/D");
  double T_vetheta;       if (generator.compare("simc")==0) Tout->Branch("vetheta", &T_vetheta, "vetheta/D");
  double T_vQ2;           if (generator.compare("simc")==0) Tout->Branch("vQ2", &T_vQ2, "vQ2/D");
  double T_vepsilon;      if (generator.compare("simc")==0) Tout->Branch("vepsilon", &T_vepsilon, "vepsilon/D"); // calculated using general eqn.
  // --
  double T_nu;            Tout->Branch("nu", &T_nu, "nu/D");
  double T_Q2;            Tout->Branch("Q2", &T_Q2, "Q2/D");
  double T_W2;            Tout->Branch("W2", &T_W2, "W2/D");
  //double T_W2_c;          Tout->Branch("W2_c", &T_W2_c, "W2_c/D"); //centered properly based on LH2 data/MC comparison
  double T_W;             Tout->Branch("W", &T_W, "W/D");
  double T_dpel;          Tout->Branch("dpel", &T_dpel, "dpel/D");
  double T_ephi;          Tout->Branch("ephi", &T_ephi, "ephi/D");
  double T_etheta;        Tout->Branch("etheta", &T_etheta, "etheta/D");
  double T_pelas;         Tout->Branch("pelas", &T_pelas, "pelas/D");
  double T_ethbend;       Tout->Branch("ethbend", &T_ethbend, "ethbend/D");
  double T_pN_exp;        Tout->Branch("pN_exp", &T_pN_exp, "pN_exp/D"); //exp. nucleon momentum
  double T_thN_exp;       Tout->Branch("thN_exp", &T_thN_exp, "thN_exp/D"); //exp. nucelon theta
  double T_phN_exp;       Tout->Branch("phN_exp", &T_phN_exp, "phN_exp/D"); //exp. nucelon phi
  double T_epsilon;       Tout->Branch("epsilon", &T_epsilon, "epsilon/D"); // calculated using general eqn.
  double T_epsilon_p;     Tout->Branch("epsilon_p", &T_epsilon_p, "epsilon_p/D");
  double T_epsilon_n;     Tout->Branch("epsilon_n", &T_epsilon_n, "epsilon_n/D");
  double T_thpq_p;        Tout->Branch("thpq_p", &T_thpq_p, "thpq_p/D");
  double T_thpq_n;        Tout->Branch("thpq_n", &T_thpq_n, "thpq_n/D");
  //track
  double T_vz;            Tout->Branch("vz", &T_vz, "vz/D");
  double T_trP;           Tout->Branch("trP", &T_trP, "trP/D");
  double T_trX;           Tout->Branch("trX", &T_trX, "trX/D");
  double T_trY;           Tout->Branch("trY", &T_trY, "trY/D");
  double T_trTh;          Tout->Branch("trTh", &T_trTh, "trTh/D");
  double T_trPh;          Tout->Branch("trPh", &T_trPh, "trPh/D");
  double T_tgX;           Tout->Branch("tgX", &T_tgX, "tgX/D");
  double T_tgY;           Tout->Branch("tgY", &T_tgY, "tgY/D");
  double T_tgTh;          Tout->Branch("tgTh", &T_tgTh, "tgTh/D");
  double T_tgPh;          Tout->Branch("tgPh", &T_tgPh, "tgPh/D");
  double T_fpX;           Tout->Branch("fpX", &T_fpX, "fpX/D");
  double T_fpY;           Tout->Branch("fpY", &T_fpY, "fpY/D");
  double T_fpTh;          Tout->Branch("fpTh", &T_fpTh, "fpTh/D");
  double T_fpPh;          Tout->Branch("fpPh", &T_fpPh, "fpPh/D");
  //cross-section predicted from MC
  double T_sigMott;       Tout->Branch("sigMott", &T_sigMott, "sigMott/D");
  double T_sigRed_p;      Tout->Branch("sigRed_p", &T_sigRed_p, "sigRed_p/D");
  double T_sigRed_n;      Tout->Branch("sigRed_n", &T_sigRed_n, "sigRed_n/D");
  double T_sigBorn_ratio; Tout->Branch("sigBorn_ratio", &T_sigBorn_ratio, "sigBorn_ratio/D");
  //BBCAL
  double T_ePS;           Tout->Branch("ePS", &T_ePS, "ePS/D"); 
  double T_rblkPS;        Tout->Branch("rblkPS", &T_rblkPS, "rblkPS/D"); 
  double T_cblkPS;        Tout->Branch("cblkPS", &T_cblkPS, "cblkPS/D"); 
  double T_idblkPS;       Tout->Branch("idblkPS", &T_idblkPS, "idblkPS/D"); 
  double T_atimePS;       Tout->Branch("atimePS", &T_atimePS, "atimePS/D"); 
  double T_eSH;           Tout->Branch("eSH", &T_eSH, "eSH/D"); 
  double T_xSH;           Tout->Branch("xSH", &T_xSH, "xSH/D"); 
  double T_ySH;           Tout->Branch("ySH", &T_ySH, "ySH/D"); 
  double T_rblkSH;        Tout->Branch("rblkSH", &T_rblkSH, "rblkSH/D"); 
  double T_cblkSH;        Tout->Branch("cblkSH", &T_cblkSH, "cblkSH/D"); 
  double T_idblkSH;       Tout->Branch("idblkSH", &T_idblkSH, "idblkSH/D"); 
  double T_atimeSH;       Tout->Branch("atimeSH", &T_atimeSH, "atimeSH/D"); 
  double T_EovP;          Tout->Branch("EovP", &T_EovP, "EovP/D"); 
  //HCAL
  double T_eHCAL;         Tout->Branch("eHCAL", &T_eHCAL, "eHCAL/D"); 
  double T_xHCAL;         Tout->Branch("xHCAL", &T_xHCAL, "xHCAL/D"); 
  double T_yHCAL;         Tout->Branch("yHCAL", &T_yHCAL, "yHCAL/D"); 
  double T_idblkHCAL;     Tout->Branch("idblkHCAL", &T_idblkHCAL, "idblkHCAL/D"); 
  double T_rblkHCAL;      Tout->Branch("rblkHCAL", &T_rblkHCAL, "rblkHCAL/D"); 
  double T_cblkHCAL ;     Tout->Branch("cblkHCAL", &T_cblkHCAL, "cblkHCAL/D"); 
  double T_atimeHCAL;     Tout->Branch("atimeHCAL", &T_atimeHCAL, "atimeHCAL/D"); 
  double T_tdcHCAL;       Tout->Branch("tdcHCAL", &T_tdcHCAL, "tdcHCAL/D"); 
  double T_xHCAL_exp;     Tout->Branch("xHCAL_exp", &T_xHCAL_exp, "xHCAL_exp/D"); 
  double T_yHCAL_exp;     Tout->Branch("yHCAL_exp", &T_yHCAL_exp, "yHCAL_exp/D");
  double T_xHCAL_exp_p;   Tout->Branch("xHCAL_exp_p", &T_xHCAL_exp_p, "xHCAL_exp_p/D"); // average sbs_kick included   
  double T_dx;            Tout->Branch("dx", &T_dx, "dx/D"); 
  double T_dy;            Tout->Branch("dy", &T_dy, "dy/D");
  double T_p_def;         Tout->Branch("p_def", &T_p_def, "p_def/D"); // expected proton deflection
  double T_ToF_n;         Tout->Branch("ToF_n", &T_ToF_n, "ToF_n/D");
  //GEM
  double T_nhitsGEM;      Tout->Branch("nhitsGEM", &T_nhitsGEM, "nhitsGEM/D");
  double T_ngoodhitsGEM;  Tout->Branch("ngoodhitsGEM", &T_ngoodhitsGEM, "ngoodhitsGEM/D");
  double T_trchi2ndf;     Tout->Branch("trchi2ndf", &T_trchi2ndf, "trchi2ndf/D");
  // Variables based on MC truth info
  double T_pN_t;          Tout->Branch("pN_t", &T_pN_t, "pN_t/D"); //TRUE nucleon momentum
  double T_thN_t;         Tout->Branch("thN_t", &T_thN_t, "thN_t/D"); //TRUE nucelon theta
  double T_phN_t;         Tout->Branch("phN_t", &T_phN_t, "phN_t/D"); //TRUE nucelon theta
  double T_xHCAL_exp_t;   Tout->Branch("xHCAL_exp_t", &T_xHCAL_exp_t, "xHCAL_exp_t/D"); 
  double T_yHCAL_exp_t;   Tout->Branch("yHCAL_exp_t", &T_yHCAL_exp_t, "yHCAL_exp_t/D");

  // Do the energy loss calculation here (only for g4sbs generator)
  double ebeam = sbsconf.GetEbeam(); // gets overwritten in the event loop

  // reading W2 offset and W cut limits
  double dy_offset = jmgr->GetValueFromSubKey<double>(key,"dy_offset");
  double W2_offset = jmgr->GetValueFromSubKey<double>(key,"W2_offset");
  std::vector<double> W_cutR; jmgr->GetVectorFromSubKey<double>(key,"W_cutR",W_cutR);
  // reading BB fiducial cut limits
  std::vector<double> bbfidu_cutR; jmgr->GetVectorFromSubKey<double>(key,"bbfidu_cutR",bbfidu_cutR);

  // HCAL cut definitions
  std::vector<double> dx_p_cut; jmgr->GetVectorFromSubKey<double>(key,"dx_p_cut", dx_p_cut);
  std::vector<double> dy_p_cut; jmgr->GetVectorFromSubKey<double>(key,"dy_p_cut", dy_p_cut);
  std::vector<double> dx_n_cut; jmgr->GetVectorFromSubKey<double>(key,"dx_n_cut", dx_n_cut);
  std::vector<double> dy_n_cut; jmgr->GetVectorFromSubKey<double>(key,"dy_n_cut", dy_n_cut);
  double sbs_kick = jmgr->GetValueFromSubKey<double>(key,"sbs_kick"); // let's data guide the amount of sbs_kick (pPeak-nPeak) - 04/21/24
  double avg_dx_np_sig = (dx_n_cut[1]+dx_p_cut[1])/2.; 
  std::vector<double> hcal_active_area = cut::hcal_active_area_simu(1,1); // Exc. 1 blk from all 4 sides
  std::vector<double> hcal_safety_margin = cut::hcal_safety_margin(avg_dx_np_sig,avg_dx_np_sig,dy_p_cut[1],hcal_active_area);
  // varying safety margin width by +/- 10% and +/- 20% in vertical direction: - means shorter margin so more stats
  std::vector<double> hcal_safety_margin_p10p = cut::hcal_safety_margin(avg_dx_np_sig*1.1,avg_dx_np_sig*1.1,dy_p_cut[1],hcal_active_area);
  std::vector<double> hcal_safety_margin_m10p = cut::hcal_safety_margin(avg_dx_np_sig*.9,avg_dx_np_sig*.9,dy_p_cut[1],hcal_active_area);
  std::vector<double> hcal_safety_margin_p20p = cut::hcal_safety_margin(avg_dx_np_sig*1.2,avg_dx_np_sig*1.2,dy_p_cut[1],hcal_active_area);
  std::vector<double> hcal_safety_margin_m20p = cut::hcal_safety_margin(avg_dx_np_sig*.8,avg_dx_np_sig*.8,dy_p_cut[1],hcal_active_area);
  TH2F *h2_xyHCAL_p = util_pd::TH2FHCALface_xy_simu("h2_xyHCAL_p",sbs_kick);
  TH2F *h2_xyHCAL_n = util_pd::TH2FHCALface_xy_simu("h2_xyHCAL_n",0);

  // costruct axes of HCAL CoS in Hall CoS
  double hcal_voffset = jmgr->GetValueFromSubKey<double>(key,"hcal_voffset");
  double hcal_hoffset = jmgr->GetValueFromSubKey<double>(key,"hcal_hoffset");
  double hcal_zoffset = jmgr->GetValueFromSubKey<double>(key,"hcal_zoffset");
  //std::vector<TVector3> HCAL_axes; kine::SetHCALaxes(sbsconf.GetSBStheta_rad(), HCAL_axes);
  std::vector<TVector3> HCAL_axes; kine::SetHCALaxes(sbsconf.GetHCALtheta_rad(), HCAL_axes);
  TVector3 HCAL_origin = (sbsconf.GetHCALdist()+hcal_zoffset)*HCAL_axes[2] + hcal_voffset*HCAL_axes[0] + hcal_hoffset*HCAL_axes[1];

  // calculating MC normalizetion factors 
  double mc_omega, lumi;
  //std::vector<double> totNtriesnCh; util_pd::GetTotNtriesnCh(sjobs,totNtriesnCh);
  std::vector<long double> totNtries, totCharge;
  bool is_simcdeeN = false;
  if (generator.compare("simc")==0 && process.compare("deeN")==0) {
    is_simcdeeN = true;
    std::vector<long double> totNtriesnCh_p; util_pd::GetTotNtriesnCh(sjobs,"deep",totNtriesnCh_p);
    std::vector<long double> totNtriesnCh_n; util_pd::GetTotNtriesnCh(sjobs,"deen",totNtriesnCh_n);
    totNtries = {totNtriesnCh_p[0],totNtriesnCh_n[0]};
    totCharge = {totNtriesnCh_p[1],totNtriesnCh_n[1]};
  } else {
    std::vector<long double> totNtriesnCh; util_pd::GetTotNtriesnCh(sjobs,totNtriesnCh);
    totNtries = {totNtriesnCh[0]}; totCharge = {totNtriesnCh[1]}; 
  }
  // -- Rejection sampling (RS) flags and weight
  bool usingRS = sjobs[0].usingRS;
  double maxwtRS = sjobs[0].maxwtRS; //gets updated per job in the event loop
  double maxwtRS_n;       //needed for summary canvas
  // report
  if (usingRS) {
    std::cout << "----\n";
    std::cout << "Using RS: " << usingRS << std::endl;
    std::cout << "Max Weight RS: " << maxwtRS << std::endl;
    if (is_simcdeeN) {
      std::cout << "Total ntries (p): " << totNtries[0] << std::endl;
      std::cout << "Total ntries (n): " << totNtries[1] << std::endl;
    }
    std::cout << "----\n";
  }

  // implementing hashtable with norm info for efficiency
  std::unordered_map<std::string,SimuJob> mnorm;
  for (auto & sjob: sjobs) mnorm[sjob.rfname] = sjob;

  // ML data
  TString data_file = "siout/hcal_all_ML.csv"; ofstream ml_data;
  double ml_arr[288]; 
  if (gen_ML_data) {
    ml_data.open(data_file);
    ml_data << "Event #, True nucleon p (GeV/c), n(0) or p(1), ADC, (ped. sub.), (pC), for, all, 288, channels, ...,\n";  
    for(int i=0; i<288; i++) {ml_arr[i] = 0.0;}
  }

  // EMFF fits
  Ye2017 yefit;
  Kelly2004 kellyfit;
  Seamus20XX seamusfit;
  Christy2022 christyfit;

  // looping through the tree ---------------------------------------
  std::cout << std::endl;
  double charge;
  long double totntries = totNtries[0];
  long nevent = 0, nevents = C->GetEntries(), ngoodevs = 0; 
  int treenum = 0, currenttreenum = 0, treeitr = 0;
  while (C->GetEntry(nevent++)) {
   
    // print progress 
    if (nevent%1000 == 0) std::cout << nevent << "/" << nevents << "\r";
    std::cout.flush();

    currenttreenum = C->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum) {
      treenum = currenttreenum;
      // apply global cuts efficiently (AJRP method)
      GlobalCut->UpdateFormulaLeaves();
      
      // getting normalization factors per job
      const char* rftemp = C->GetFile()->GetName();
      SimuJob sjtemp = mnorm[rftemp];
      ebeam = sjtemp.ebeam; charge = sjtemp.charge;
      lumi = sjtemp.lumi; mc_omega = sjtemp.genvol; 
      if (is_simcdeeN) {
	if (sjtemp.process.compare("deep")==0) totntries = totNtries[0]; //p events
	else {totntries = totNtries[1]; maxwtRS_n = sjtemp.maxwtRS;} // n events
      }
      if (sjtemp.usingRS) maxwtRS = sjtemp.maxwtRS;
    } 
    bool passedgCut = GlobalCut->EvalInstance(0) != 0;   
    if (!passedgCut) continue;
    ngoodevs++;

    // cross section weighted normalization factor
    weight = usingRS ? maxwtRS*mc_omega*lumi/totntries : mc_sigma*mc_omega*lumi/totntries;
    weight_norm = usingRS ? maxwtRS*mc_omega*lumi/totntries/charge : mc_sigma*mc_omega*lumi/totntries/charge;

    // kinematic parameters
    double ebeam_corr = ebeam; //- MeanEloss;
    if (generator.compare("simc")==0) ebeam_corr = mc_ebeam;
    double precon = p[0]; //+ MeanEloss_outgoing

    // constructing the 4 vectors
    /* Reaction    : e + e' -> N + N'
       Conservation: Pe + Peprime = PN + PNprime */
    TVector3 vertex(0, 0, vz[0]);
    TLorentzVector Pe(0,0,ebeam_corr,ebeam_corr);   // incoming e-
    TLorentzVector Peprime(px[0] * (precon/p[0]),   // scattered e-
    			   py[0] * (precon/p[0]),
    			   pz[0] * (precon/p[0]),
    			   precon);                 
    TLorentzVector PN;                              // target nucleon [Ntype ??]
    kine::SetPN(Ntype, PN);
    TLorentzVector PNprime;                         // Recoil nucleon 4-vector
    TLorentzVector q = Pe - Peprime;                // 4-momentum of virtual photon

    double etheta = kine::etheta(Peprime);
    double ephi = kine::ephi(Peprime);
    double pelas = kine::pelas(ebeam_corr, etheta, Ntype);
    double thelas = kine::thelas(ebeam_corr, precon, Ntype);

    double nu = 0.;                   // energy of the virtual photon
    double pN_expect = 0.;            // expected recoil nucleon momentum
    double thetaN_expect = 0.;        // expected recoil nucleon theta
    double phiN_expect = ephi + constant::pi; 
    /* Different modes of calculation. Goal is to achieve the best resolution
       model 0 = uses reconstructed p as independent variable
       model 1 = uses reconstructed angles as independent variable 
       model 2 = uses 4-vector calculation */
    TVector3 pNhat;                   // 3-momentum of the recoil nucleon (Unit)
    double Q2recon{0},W2recon{0},epsilon{0},epsilon_p{0},epsilon_n{0};
    if (model == 0) {
      nu = Pe.E() - Peprime.E();
      pN_expect = kine::pN_expect(nu, Ntype);
      thetaN_expect = acos((Pe.E() - Peprime.Pz()) / pN_expect);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);
      PNprime.SetPxPyPzE(pN_expect*pNhat.X(), pN_expect*pNhat.Y(), pN_expect*pNhat.Z(), nu+PN.E());
      Q2recon = kine::Q2(Pe.E(), Peprime.E(), thelas);
      epsilon = kine::epsilon_general(thelas,Q2recon,nu);
      epsilon_p = kine::epsilon(thelas,Q2recon,"p");
      epsilon_n = kine::epsilon(thelas,Q2recon,"n");
      //W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype);
      W2recon = kine::W2_general(Pe.E(), Peprime.E(), etheta, Ntype);
     } else if (model == 1) {
      nu = Pe.E() - pelas;
      pN_expect = kine::pN_expect(nu, Ntype);
      thetaN_expect = acos((Pe.E() - pelas*cos(etheta)) / pN_expect);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);
      PNprime.SetPxPyPzE(pN_expect*pNhat.X(), pN_expect*pNhat.Y(), pN_expect*pNhat.Z(), nu+PN.E());
      Q2recon = kine::Q2(Pe.E(), pelas, etheta);
      epsilon = kine::epsilon_general(etheta,Q2recon,nu);
      epsilon_p = kine::epsilon(etheta,Q2recon,"p");
      epsilon_n = kine::epsilon(etheta,Q2recon,"n");
      //W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype);
      W2recon = kine::W2_general(Pe.E(), Peprime.E(), etheta, Ntype);
    } else if (model == 2) {
      nu = q.E();
      PNprime = q + PN;
      pNhat = PNprime.Vect().Unit();
      Q2recon = -q.M2();
      epsilon = kine::epsilon_general(etheta,Q2recon,nu);
      epsilon_p = kine::epsilon(etheta,Q2recon,"p");
      epsilon_n = kine::epsilon(etheta,Q2recon,"n");
      W2recon = PNprime.M2();
      //--
      pN_expect = kine::pN_expect(nu, Ntype);
    }
    h_Q2->Fill(Q2recon); 
    double Wrecon = sqrt(max(0., W2recon));
    double dpel = Peprime.E()/pelas - 1.0; h_dpel->Fill(dpel);

    T_ebeam = Pe.E();

    // at vertex ("true" values)
    T_veE = mc_veE;
    T_vetheta = mc_vetheta;
    T_vQ2 = kine::Q2(Pe.E(),T_veE,T_vetheta);
    T_vepsilon = kine::epsilon_general(T_vetheta,T_vQ2,T_ebeam-T_veE);

    T_nu = nu;
    T_Q2 = Q2recon;
    T_W2 = W2recon+W2_offset;
    //zT_W2_c = W2recon-W2_offset;
    T_W = Wrecon;
    T_dpel = dpel;
    T_ephi = ephi;
    T_etheta = etheta;
    T_pelas = pelas;
    T_pN_exp = pN_expect;
    T_thN_exp = thetaN_expect;
    T_epsilon = epsilon;
    T_epsilon_p = epsilon_p;
    T_epsilon_n = epsilon_n;

    // defining W cut
    WCut = T_W >= W_cutR[0] && T_W <= W_cutR[1];

    T_vz = vz[0];
    T_trP = p[0];
    T_trX = xTr[0];
    T_trY = yTr[0];
    T_trTh = thTr[0];
    T_trPh = phTr[0];

    T_tgX = xtgt[0];
    T_tgY = ytgt[0];
    T_tgTh = thtgt[0];
    T_tgPh = phtgt[0];

    T_fpX = xfp[0];
    T_fpY = yfp[0];
    T_fpTh = thfp[0];
    T_fpPh = phfp[0];

    // calculating bend angle
    T_ethbend = kine::ethbend(T_tgTh,T_tgPh,T_fpTh,T_fpPh);
    
    // defining BB fiducial cut
    bbfiduCut = abs(T_fpX - 0.9*T_fpTh - bbfidu_cutR[0]) <= bbfidu_cutR[1];

    // EMFF extraction using same parametrization used in MC generators
    double GEp_kelly = kellyfit.GetFF(G_t::kGEp,Q2recon);
    double GMp_kelly = kellyfit.GetFF(G_t::kGMp,Q2recon);
    double GEn_seamus = seamusfit.GetFF(G_t::kGEn,Q2recon);
    double GMn_kelly = kellyfit.GetFF(G_t::kGMn,Q2recon);

    T_sigMott = kine::sigmaMott(Pe.E(),Peprime.E(),etheta);
    T_sigRed_p = kine::sigmaReduced(kine::tau(Q2recon,"p"),epsilon_p,GEp_kelly,GMp_kelly);
    T_sigRed_n = kine::sigmaReduced(kine::tau(Q2recon,"n"),epsilon_n,GEn_seamus,GMn_kelly);
    T_sigBorn_ratio = kine::sigmaBorn_ratio(etheta,Q2recon,GEp_kelly,GMp_kelly,GEn_seamus,GMn_kelly);

    T_ePS = ePS;
    T_rblkPS = rblkPS;
    T_cblkPS = cblkPS;
    T_idblkPS = idblkPS;
    T_atimePS = atimePS;

    T_eSH = eSH;
    T_xSH = xSH;
    T_ySH = ySH;
    T_rblkSH = rblkSH;
    T_cblkSH = cblkSH;
    T_idblkSH = idblkSH;
    T_atimeSH = atimeSH;

    T_EovP = (eSH+ePS)/p[0];

    T_eHCAL = eHCAL;
    T_xHCAL = xHCAL;
    T_yHCAL = yHCAL;
    T_rblkHCAL = rblkHCAL;
    T_cblkHCAL = cblkHCAL;
    T_idblkHCAL = idblkHCAL;
    T_atimeHCAL = atimeHCAL;
    T_tdcHCAL = tdcHCAL[0];

    T_nhitsGEM = nhitsGEM[0];
    T_ngoodhitsGEM = ngoodhitsGEM[0];
    T_trchi2ndf = trchi2ndf[0];

    T_mc_fnucl = int(mc_fnucl);

    // Expected position of the q vector at HCAL
    std::vector<double> xyHCAL_exp; // xyHCAL_exp[0] = xHCAL_exp & xyHCAL_exp[1] = yHCAL_exp
    kine::GetxyHCALexpect(vertex, pNhat, HCAL_origin, HCAL_axes, xyHCAL_exp);
    double dx = xHCAL - xyHCAL_exp[0];  
    double dy = yHCAL - xyHCAL_exp[1]; 

    T_xHCAL_exp = xyHCAL_exp[0];
    T_yHCAL_exp = xyHCAL_exp[1];
    T_dx = dx;
    T_dy = dy+dy_offset;

    /* Calculating thpq (both p & n hypothesis) */
    // n (no deflection)
    TVector3 HCAL_pos = HCAL_origin + xHCAL*HCAL_axes[0] + yHCAL*HCAL_axes[1];
    TVector3 n_dir = (HCAL_pos - vertex);
    T_thpq_n = acos(n_dir.Unit().Dot(pNhat));
    // p 
    double BdL = sbsscalefield * expconst::sbsmaxfield_simu * expconst::sbsdipolegap;
    double proton_thetabend = 0.3 * BdL / PNprime.Vect().Mag();  // p*theta = 0.3*BdL
    double proton_deflection = tan(proton_thetabend)*(sbsconf.GetHCALdist()+hcal_zoffset-(sbsconf.GetSBSdist()+expconst::sbsdipolegap/2.0));
    T_p_def = proton_deflection;
    TVector3 p_dir = (HCAL_pos + proton_deflection*HCAL_axes[0] - vertex);
    T_thpq_p = acos(p_dir.Unit().Dot(pNhat));

    // calculate ToF for neutrons
    double ToF_n = (n_dir.Mag() / constant::c) * sqrt(1. + pow((constant::Mn/PNprime.Vect().Mag()), 2));
    T_ToF_n = ToF_n*1e9; //ns

    /* Calculating Weight Factor for HCAL Efficiency (NDE) */
    // constructing TRUE q-vector
    TVector3 vertex_t(mc_vx, mc_vy, mc_vz);
    T_pN_t = mc_np;
    T_thN_t = acos(mc_npz / mc_np); // final state N's theta
    T_phN_t = TMath::PiOver2() - atan2(mc_npy, mc_npx);  // final state N's phi
    TVector3 pNhat_t = kine::qVect_unit(T_thN_t,T_phN_t);
    // TRUE position of the q vector at HCAL
    vector<double> xyHCAL_exp_t; // xyHCAL_exp[0] = xHCAL_exp & xyHCAL_exp[1] = yHCAL_exp
    kine::GetxyHCALexpect(vertex_t, pNhat_t, HCAL_origin, HCAL_axes, xyHCAL_exp_t);
    T_xHCAL_exp_t = -xyHCAL_exp_t[0];
    T_yHCAL_exp_t = xyHCAL_exp_t[1];
    // Get the efficiency correction value
    int effibin_n = h2_effi_map->FindBin(T_yHCAL_exp_t, T_xHCAL_exp_t);
    int effibin_p = h2_effi_map->FindBin(T_yHCAL_exp_t, T_xHCAL_exp_t-T_p_def);
    double effi_corr = mc_fnucl==0 ? h2_effi_map->GetBinContent(effibin_n)/avg_effi : h2_effi_map->GetBinContent(effibin_p)/avg_effi;
    // updating weight factors accordingly
    weight_effic = weight*effi_corr;
    weight_norm_effic = weight_norm*effi_corr;
    // --    
    
    // HCAL active area and safety margin cuts [Fiducial region]
    ARCut = cut::inHCAL_activeA(xHCAL,yHCAL,hcal_active_area);
    SMCut = cut::inHCAL_safety_margin(target,xyHCAL_exp[0],xyHCAL_exp[1],sbs_kick,hcal_safety_margin);
    fiduCut = ARCut && SMCut; 
    // defining HCAL cuts
    pCut = cut::SpotCut(dx,dx_p_cut[0],dx_p_cut[1],dx_p_cut[2],dy,dy_p_cut[0],dy_p_cut[1],dy_p_cut[2]);
    nCut = cut::SpotCut(dx,dx_n_cut[0],dx_n_cut[1],dx_n_cut[2],dy,dy_n_cut[0],dy_n_cut[1],dy_n_cut[2]);
    // ***
    pdx_nS = fabs(dx-dx_p_cut[0])/dx_p_cut[1];
    ndx_nS = fabs(dx-dx_n_cut[0])/dx_n_cut[1];
    dy_nS = fabs(dy-dy_p_cut[0])/dy_p_cut[1];  // assuming dy is same for n and p 

    // W cut
    if (WCut&&bbfiduCut&&T_eHCAL>0) {
      if (abs(dy)<0.4) h_dxHCAL_nfc->Fill(dx, weight);
      h_dyHCAL_nfc->Fill(dy, weight);
      // fiducial cut
      if (fiduCut) {
    	h_dxHCAL->Fill(dx, weight);
    	h_dyHCAL->Fill(dy, weight);
    	// dx dist. for p & n separately using MC info
    	if (int(mc_fnucl)==0) {
	  if (abs(dy)<0.4) {
	    h_dxHCAL->Fill(dx, weight);
	  }
	  h2_xyHCAL_n->Fill(xyHCAL_exp[1], xyHCAL_exp[0], weight);
    	} else if (int(mc_fnucl)==1) {
	  if (abs(dy)<0.4) {
	    h_dxHCAL_p->Fill(dx, weight);
	    h_dx_w_p_def->Fill(dx+T_p_def, weight);
	  }	  
	  h2_xyHCAL_p->Fill(xyHCAL_exp[1], xyHCAL_exp[0] - sbs_kick, weight);
    	} else {
    	  std::cerr << "*!* Invalid final state nuclei!" << std::endl; 
    	  std::exit(1);
    	}

    	h2_rcHCAL->Fill(cblkHCAL, rblkHCAL, weight);
    	// p & n spots
    	h2_dxdyHCAL->Fill(dy, dx, weight);

    	// hit map to show p & n in fiducial region
    	// if (pCut) h2_xyHCAL_p->Fill(xyHCAL_exp[1], xyHCAL_exp[0] - sbs_kick);
    	// if (nCut) h2_xyHCAL_n->Fill(xyHCAL_exp[1], xyHCAL_exp[0]);
      }
    }

    // fiducial cut but no W cut
    if (Wrecon>0) { 
      if (fiduCut) {
	h_W->Fill(Wrecon,weight);
	if (pCut || nCut) { 
	  h_W_cut->Fill(Wrecon,weight);
	  h_W2_cut->Fill(T_W2,weight);
	  h_W2_cut_noOff->Fill(W2recon,weight);
	} else {
	  h_W_acut->Fill(Wrecon,weight);
	}
      }
    }
      
    // ML data
    if (gen_ML_data) {
      //ml_data << nevent << "," << mc_fnucl << "," << mc_np << "," << dx << "," << dy << ",";
      //ml_data << nevent << "," << mc_fnucl << "," << mc_np << "," << xyHCAL_exp[0] << "," << xyHCAL_exp[1] << ",";
      ml_data << nevent << "," << mc_fnucl << "," << mc_np << "," << dx << "," << dy << "," << xyHCAL_exp[0] << "," << xyHCAL_exp[1] << ",";
      for (int ihit=0;ihit<aHCALndata;ihit++) {ml_arr[int(elHCAL[ihit])] = apHCAL[ihit];}
      for (int i=0; i<288; i++) {ml_data << ml_arr[i] << ","; ml_arr[i] = 0.0; }
      ml_data << std::endl;
    }

    Tout->Fill();
  } // event loop
  std::cout << std::endl << std::endl;

  /////////////////////////////////
  // Generating diagnostic plots //
  /////////////////////////////////
  TString outPlot = outFile; outPlot.ReplaceAll(".root",".pdf");
  /**** Canvas 1 (p&n spots) ****/
  TCanvas *c1 = util_pd::TC("c1",2,2);
  c1->cd(1); //
  h2_dxdyHCAL->Draw("colz");
  TEllipse Ep_p;
  Ep_p.SetFillStyle(0); Ep_p.SetLineColor(2); Ep_p.SetLineWidth(2);
  Ep_p.DrawEllipse(dy_p_cut[0], dx_p_cut[0], dy_p_cut[2]*dy_p_cut[1], dx_p_cut[2]*dx_p_cut[1], 0,360,0);
  TEllipse Ep_n;
  Ep_n.SetFillStyle(0); Ep_n.SetLineColor(3); Ep_n.SetLineWidth(2);
  Ep_n.DrawEllipse(dy_n_cut[0], dx_n_cut[0], dy_n_cut[2]*dy_n_cut[1], dx_n_cut[2]*dx_n_cut[1], 0,360,0);
  c1->cd(2); //
  h_W->Draw(); h_W->SetLineColor(1);
  h_W_cut->Draw("same"); h_W_cut->SetLineColor(2);
  h_W_acut->Draw("same");
  c1->cd(3); //
  h2_xyHCAL_p->Draw("colz");
  util_pd::DrawArea(hcal_active_area,2,4,9);
  util_pd::DrawArea(hcal_safety_margin,4,4,9);
  c1->cd(4); //
  h2_xyHCAL_n->Draw("colz");
  util_pd::DrawArea(hcal_active_area,2,4,9);
  util_pd::DrawArea(hcal_safety_margin,4,4,9);
  c1->SaveAs(Form("%s[",outPlot.Data())); c1->SaveAs(Form("%s",outPlot.Data())); c1->Write();
  //**** -- ***//

  /**** Canvas 2 (dx & dy) ****/
  TCanvas *c2 = util_pd::TC("c2",2,2);
  std::vector<double> hdxp_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_dxHCAL_p_fitR",hdxp_fitR);
  std::vector<double> hdxn_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_dxHCAL_n_fitR",hdxn_fitR);
  std::vector<double> hdy_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_dyHCAL_fitR",hdy_fitR);
  gStyle->SetOptFit(1);
  c2->cd(1); //
  gPad->SetGridx();
  TF1 *fdxp = fit::fit_1gs_nbg(hdxp_fitR,h_dxHCAL);
  TF1 *fdxn = fit::fit_1gs_nbg(hdxn_fitR,h_dxHCAL);
  double dxpM = fdxp->GetParameter(1); double dxpS = fdxp->GetParameter(2);
  double dxnM = fdxn->GetParameter(1); double dxnS = fdxn->GetParameter(2);
  fdxp->Draw("same");
  c2->cd(2); //
  gPad->SetGridx();
  TF1 *fdy = fit::fit_1gs_nbg(hdy_fitR,h_dyHCAL);
  double dyM = fdy->GetParameter(1); double dyS = fdy->GetParameter(2);
  c2->cd(3); //
  gPad->SetGridx();
  TF1 *fdxp_nfc = fit::fit_1gs_nbg(hdxp_fitR,h_dxHCAL_nfc);
  TF1 *fdxn_nfc = fit::fit_1gs_nbg(hdxn_fitR,h_dxHCAL_nfc);
  double dxpM_nfc = fdxp_nfc->GetParameter(1); double dxpS_nfc = fdxp_nfc->GetParameter(2);
  double dxnM_nfc = fdxn_nfc->GetParameter(1); double dxnS_nfc = fdxn_nfc->GetParameter(2);
  fdxp_nfc->Draw("same");
  c2->cd(4); //
  gPad->SetGridx();
  TF1 *fdy_nfc = fit::fit_1gs_nbg(hdy_fitR,h_dyHCAL_nfc);
  double dyM_nfc = fdy->GetParameter(1); double dyS_nfc = fdy->GetParameter(2);
  c2->SaveAs(Form("%s",outPlot.Data())); c2->Write();
  //**** -- ***//

  /**** Canvas 3 (Deflection and W2) ****/
  TCanvas *c3 = util_pd::TC("c3",1,2);
  c3->SetGridx();
  gStyle->SetOptFit(1);
  //
  c3->cd(1);
  gPad->SetGridx();
  std::vector<double> hpdef_fitR{-0.5,0.5,1.2,1.2};
  TF1 *fpdef = fit::fit_1gs_nbg(hpdef_fitR,h_dx_w_p_def);
  double pdefM = fpdef->GetParameter(1); double pdefS = fpdef->GetParameter(2);
  h_dx_w_p_def->Draw("same");
  h_dx_w_p_def->SetLineColor(kBlack);
  //
  c3->cd(2);
  gPad->SetGridx();
  std::vector<double> hW2_fitR{0.6,1.2,1.2,1.2};
  TF1 *fW2 = fit::fit_1gs_nbg(hW2_fitR,h_W2_cut);
  double W2M = fW2->GetParameter(1); double W2S = fW2->GetParameter(2);
  h_W2_cut->Draw("same");
  h_W2_cut->SetLineColor(kBlack);
  h_W2_cut_noOff->Draw("same");
  h_W2_cut_noOff->SetLineColor(kRed);
  c3->SaveAs(Form("%s",outPlot.Data())); c3->Write();
  //**** -- ***//    

  /**** Summary Canvas ****/
  TCanvas *cSummary = new TCanvas("cSummary","Summary");
  cSummary->cd();
  TPaveText *pt = new TPaveText(.05,.1,.95,.8);
  pt->AddText(Form(" Date of creation: %s", util_pd::getDate().c_str()));
  pt->AddText(Form("Configfile: %s",configfilename));
  pt->AddText(Form(" Analyzing %s generated QE events for SBS%d-SBS%dp settings",generator.c_str(),conf,sbsmag));
  pt->AddText(Form(" Analysis model: %d",model));
  pt->AddText(Form(" HCAL offsets: v = %.4f, h = %.4f, z = %.4f",hcal_voffset,hcal_hoffset,hcal_zoffset));
  pt->AddText(Form(" W2 offset = %.4f, dy offset = %.4f",W2_offset,dy_offset));
  pt->AddText(Form(" Total # events analyzed: %ld",nevents));
  pt->AddText(Form(" Global cuts: "));
  std::string tmpstr = "";
  for (std::size_t i=0; i<gCutList.size(); i++) {
    if (i>0 && i%3==0) {pt->AddText(Form(" %s",tmpstr.c_str())); tmpstr="";}
    tmpstr += gCutList[i] + ", "; 
  }
  if (!tmpstr.empty()) pt->AddText(Form(" %s",tmpstr.c_str()));
  pt->AddText(Form(" # events passed global cuts: %ld", ngoodevs));
  pt->AddText(" Elastic cuts: ");
  pt->AddText(Form(" Inbuilt W cut: %.2f #leq W #leq %.2f GeV/c",W_cutR[0],W_cutR[1]));
  pt->AddText(Form(" Inbuilt p cut (#Deltax): Mean = %.4f, %.1f#sigma = %.4f",dx_p_cut[0],dx_p_cut[2],dx_p_cut[1]));
  pt->AddText(Form(" Inbuilt p cut (#Deltay): Mean = %.4f, %.1f#sigma = %.4f",dy_p_cut[0],dy_p_cut[2],dy_p_cut[1]));
  pt->AddText(Form(" Inbuilt n cut (#Deltax): Mean = %.4f, %.1f#sigma = %.4f",dx_n_cut[0],dx_n_cut[2],dx_n_cut[1]));
  pt->AddText(Form(" Inbuilt n cut (#Deltay): Mean = %.4f, %.1f#sigma = %.4f",dy_n_cut[0],dy_n_cut[2],dy_n_cut[1]));
  pt->AddText(Form(" Inbuilt BB fiducial cut: |fpX-0.9*fpTh-%.2f| #leq %.2f",bbfidu_cutR[0],bbfidu_cutR[1]));
  pt->AddText(" Fit info: ");
  pt->AddText(" p & n peaks, w/ fiducial cut: dxpM,dxpS,dxnM,dxnS,dyM,dyS ");
  pt->AddText(Form(" %.5f,%.5f,%.5f,%.5f,%.5f,%.5f",dxpM,dxpS,dxnM,dxnS,dyM,dyS));
  pt->AddText(" p & n peaks, w/o fiducial cut: dxpM_nfc,dxpS_nfc,dxnM_nfc,dxnS_nfc,dyM_nfc,dyS_nfc ");
  pt->AddText(Form(" %.5f,%.5f,%.5f,%.5f,%.5f,%.5f",dxpM_nfc,dxpS_nfc,dxnM_nfc,dxnS_nfc,dyM_nfc,dyS_nfc));
  if (usingRS) {
    pt->AddText(Form(" Rejection Sampling (RS) Summary:")); 
    if (is_simcdeeN) { 
      pt->AddText(Form(" Chosen maximum weight (p): %f",maxwtRS)); 
      pt->AddText(Form(" Total # tries (p): %.0Lf",totNtries[0])); 
      pt->AddText(Form(" Chosen maximum weight (n): %f",maxwtRS_n)); 
      pt->AddText(Form(" Total # tries (n): %.0Lf",totNtries[1])); 
    } else {
      pt->AddText(Form(" Chosen maximum weight: %f",maxwtRS)); 
      pt->AddText(Form(" Total # tries: %.0Lf",totNtries[0])); 
    }
    TText *t5 = pt->GetLineWith(" Rejection"); t5->SetTextColor(kMagenta+2);
  }
  sw->Stop();
  pt->AddText(Form("Macro processing time: CPU %.1fs | Real %.1fs",sw->CpuTime(),sw->RealTime()));
  TText *t1 = pt->GetLineWith("Configfile"); t1->SetTextColor(kRed);
  TText *t2 = pt->GetLineWith(" Global"); t2->SetTextColor(kBlue);
  TText *t3 = pt->GetLineWith(" Elastic"); t3->SetTextColor(kBlue);
  TText *t4 = pt->GetLineWith(" Fit info"); t4->SetTextColor(kBlue);
  TText *t6 = pt->GetLineWith("Macro"); t6->SetTextColor(kGreen+3);
  pt->Draw(); 
  cSummary->SaveAs(Form("%s",outPlot.Data())); cSummary->SaveAs(Form("%s]",outPlot.Data())); cSummary->Write();  
  //**** -- ***//

  std::cout << "\n----Fit info----" << "\n";
  std::cout << "dxpM,dxpS,dxnM,dxnS,dyM,dyS,dxpM_nfc,dxpS_nfc,dxnM_nfc,dxnS_nfc,dyM_nfc,dyS_nfc" << "\n";
  std::cout << dxpM<<","<<dxpS<<","<<dxnM<<","<<dxnS<<","<<dyM<<","<<dyS<<","<<dxpM_nfc<<","<<dxpS_nfc<<","<<dxnM_nfc<<","<<dxnS_nfc<<","<<dyM_nfc<<","<<dyS_nfc<<"\n";
  std::cout << "----------------" << "\n\n"; 

  std::cout << "------" << "\n";
  std::cout << " Summary plots  : " << outPlot << std::endl;
  std::cout << " Output ROOT file  : " << outFile << std::endl;
  std::cout << "------" << "\n\n";

  std::cout << "CPU time = " << sw->CpuTime() << " s. Real time = " << sw->RealTime() << " s.\n\n";

  Tout->Write("",TObject::kOverwrite);
  h_Q2->Write();
  h_dpel->Write(); h_W->Write();
  h_W_cut->Write(); h_W_acut->Write();
  h_W2_cut->Write(); h_W2_cut_noOff->Write();  
  h_dxHCAL->Write(); h_dxHCAL_nfc->Write();
  h_dyHCAL->Write(); h_dyHCAL_nfc->Write();
  h_dxHCAL_p->Write(); h_dxHCAL_n->Write();
  h2_rcHCAL->Write(); h2_dxdyHCAL->Write();
  h2_xyHCAL_p->Write(); h2_xyHCAL_n->Write();
  h_dx_w_p_def->Write();
  h2_effi_map->Write();
  sw->Delete();
  delete jmgr;
  return 0;
}


/*
  ////////////////////////////////////////////////////////
  // Brief description of configuration file parameters //
  ////////////////////////////////////////////////////////
  ** Ntype : Struck nucleon type. (Valid options: n, p, np)
     - n(p) => neutrons(protons) 
     - np => Avg. of n and p masses. Use for LD2 data since the struck nucleon is not known apriori.
  ** SBS_config : SBS configuration. (Valid options: 4,7,11,14,8,9)
  ** SBS_magnet_percent : SBS magnet field current as a percentage of 2100A
  ** model : Model of analysis. (Valid options: 0, 1, 2)
     - 0 => uses reconstructed p as independent variable
     - 1 => uses reconstructed angles as independent variable
     - 2 => uses 4-vector calculation
  ** generator : MC enevt generator (Valid options: g4sbs, simc)
  ** prefix_to_filebase : <prefix>_sbs<sbsconfig>_sbs<sbsmagfield>p_<generator>_<process>.root
  ** Njobs_to_ana : # MC jobs to analyze
  ** rootfile_dir : Directory name w/ path containing the MC ROOT files to analyze
  ** global_cut : set of global cuts to apply at the start of event processing
  ** SBS_field : 
  ** hcal_v(h)offset :
  ** dx_p(n)_cut : deltax p(n) peak cut definitions. Needed to constitute p(n) spot cuts. 
     - dx_p(n)_cut[0] => p(n) peak position, 
     - dx_p(n)_cut[1] => p(n) peak RMS, 
     - dx_p(n)_cut[2] => # sigma to include in the cut
  ** dy_p(n)_cut : deltay p(n) peak cut definitions. Needed to constitute p(n) spot cuts.
     - dy_p(n)_cut[0] => p(n) peak position, 
     - dy_p(n)_cut[1] => p(n) peak RMS, 
     - dy_p(n)_cut[2] => # sigma to include in the cut
  ** h_dx(dy)HCAL_lims : h_dx(dy)HCAL histogram limits (Can be found in the output ROOT file)
     - h_dx(dy)HCAL_lims[0] : No. of bins of h_dx(dy)HCAL histograms
     - h_dx(dy)HCAL_lims[1] : xmin
     - h_dx(dy)HCAL_lims[2] : xmax
  ** h_dxHCAL_p(n)_fitR : Fit ranges for the proton(neutron) signal peak in deltax dist. (h_dxHCAL histogram)
                          Algorithm fits the distribution twice for optimization.
     - h_dxHCAL_p(n)_fitR[0] : xmin for 1st fit (crude). Try to avoid any secondary peak.
     - h_dxHCAL_p(n)_fitR[1] : xmax for 1st fit (crude). Try to avoid any secondary peak.
     - h_dxHCAL_p(n)_fitR[2] : # sigma below the peak for 2nd fit (fine)
     - h_dxHCAL_p(n)_fitR[3] : # sigma above the peak for 2nd fit (fine)
  ** h_dyHCAL_fitR : Fit ranges for the signal peak in deltay dist. (h_dxHCAL histogram)
                     Algorithm fits the distribution twice for optimization.
     - h_dyHCAL_fitR[0] : xmin for 1st fit (crude). Try to avoid any secondary peak.
     - h_dyHCAL_fitR[1] : xmax for 1st fit (crude). Try to avoid any secondary peak.
     - h_dyHCAL_fitR[2] : # sigma below the peak for 2nd fit (fine)
     - h_dyHCAL_fitR[3] : # sigma above the peak for 2nd fit (fine)
  ** W_cutR : W cut range.
     - W_cutR[0] : lower limit
     - W_cutR[1] : upper limit
*/
